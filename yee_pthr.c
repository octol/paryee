/*
 *    Copyright (C) 2012 Jon Haggblad
 *
 *    This file is part of ParYee.
 *
 *    ParYee is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ParYee is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ParYee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>

#include "yee_common.h"

/* 
 * NOTE: Barriers are defined in the optional part of the POSIX standard,
 * hence with gcc one cannot compile with -ansi / -std=cXX 
 */
pthread_barrier_t barrier;

/* 
 * Data structures
 */
struct thread_param {
    double Nt;
    double dt;
    struct field_variable *p;
    struct field_variable *u;
    struct partition part;
};

void *thread_main(void *arg)
{
    struct thread_param *param = (struct thread_param *) arg;
    struct field_variable *p = param->p;;
    struct field_variable *u = param->u;
    long n;
    long i;
    double dt = param->dt;
    double dx = p->dx;

    for (n = 0; n < param->Nt; ++n) {
        /* update the pressure (p) */
        for (i = param->part.p[0]; i <= param->part.p[1]; ++i)
            p->value[i] += dt / dx * (u->value[i + 1] - u->value[i]);
        pthread_barrier_wait(&barrier);
        /* update the velocity (u) */
        for (i = param->part.u[0]; i <= param->part.u[1]; ++i)
            u->value[i] += dt / dx * (p->value[i] - p->value[i - 1]);
        pthread_barrier_wait(&barrier);
    }

    return NULL;
}

/* 
 * Main program 
 */
int main(int argc, char *argv[])
{
    /* Parameters */
    double length = 1;
    double cfl = 1;
    double T = 0.3;
    double c = 1;
    long nx = 2048;
    double tic, toc;
    struct field f;
    long threads = 4;
    long cells_per_thread;
    /*struct Partition part; */
    pthread_t *thr;
    pthread_attr_t attr;
    struct thread_param *param;
    char outfile_p[STR_SIZE] = "yee_pthr_p.tsv";
    char outfile_u[STR_SIZE] = "yee_pthr_u.tsv";

    /* Parse parameters from commandline */
    parse_cmdline(&nx, &threads, outfile_p, outfile_u, argc, argv);
    printf("Running with: N=%ld, threads=%ld\n", nx, threads);

    /* Initialize */
    alloc_field(&f.p, nx);
    alloc_field(&f.u, nx + 1);
    set_grid(&f.p, 0.5 * length / nx, length - 0.5 * length / nx);
    set_grid(&f.u, 0, length);
    apply_func(&f.p, gauss);    /* initial data */
    apply_func(&f.u, zero);     /* initial data */

    /* Depends on the numerical variables initialized above */
    f.dt = cfl * f.p.dx / c;    /* CFL condition is: c*dt/dx = cfl <= 1 */
    f.Nt = T / f.dt;

    /* Setup threads */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_barrier_init(&barrier, NULL, threads);
    thr = malloc(sizeof(pthread_t) * (threads - 1));
    if (!thr) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    /* Partion the grid and assemble structs to be used as arguments to the
     * threads */
    param = malloc(sizeof(struct thread_param) * threads);
    if (!param) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    cells_per_thread = nx / threads;
    assert(cells_per_thread * threads == nx);
    unsigned int i;
    for (i = 0; i < threads; ++i) {
        param[i].part = partition_grid(i, cells_per_thread);
        param[i].Nt = f.Nt;
        param[i].dt = f.dt;
        param[i].p = &f.p;
        param[i].u = &f.u;
    }

    /* timestep */
    tic = gettime();

    /* Spawn additional [threads-1] threads */
    for (i = 0; i < threads - 1; ++i)
        pthread_create(&thr[i], &attr, thread_main, &param[i]);

    /* Current thread */
    thread_main(&param[threads - 1]);

    /* Collect threads when finished */
    for (i = 0; i < threads - 1; ++i)
        pthread_join(thr[i], NULL);

    toc = gettime();
    printf("Elapsed: %f seconds\n", toc - tic);

    /* write data to disk and free data */
    write_to_disk(f.p, outfile_p);
    write_to_disk(f.u, outfile_u);
    free(param);
    free(thr);
    free_field(f.p);
    free_field(f.u);

    pthread_attr_destroy(&attr);
    pthread_barrier_destroy(&barrier);
    return EXIT_SUCCESS;
}
