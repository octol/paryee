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

/*
 * Pthread shared memory parallelization
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

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
    struct field *f;
    struct cell_partition part;
    long tid;
};

void *thread_main(void *arg)
{
    struct thread_param *param = (struct thread_param *) arg;
    struct field *f = param->f;
    long tid = param->tid;

    long p0, p1, u0, u1, v0, v1;
    cellindex_to_nodeindex(tid, param->part, &p0, &p1, &u0, &u1, &v0, &v1);
#ifdef DEBUG
    printf("tid=%lu  p0=%lu  p1=%lu  u0=%lu  u1=%lu  v0=%lu  v1=%lu\n",
           tid, p0, p1, u0, u1, v0, v1);
#endif

    /* for the index macro */
    double *p = f->p.value;
    double *u = f->u.value;
    double *v = f->v.value;
    long nx = f->p.size_x;
    long ny = f->p.size_y;

    double dt = f->dt;
    long n, i, j;
    for (n = 0; n < f->Nt; ++n) {

        /* update the pressure (p) */
        for (i = p0; i < p1; ++i) {
            for (j = 0; j < ny; ++j) {
                P(i, j) +=
                    dt / f->u.dx * (U(i + 1, j) - U(i, j)) +
                    dt / f->v.dy * (V(i, j + 1) - V(i, j));
            }
        }
        pthread_barrier_wait(&barrier);

        /* update the velocity (u,v) */
        for (i = u0; i < u1; ++i) {
            for (j = 0; j < ny; ++j) {
                U(i, j) += dt / f->p.dx * (P(i, j) - P(i - 1, j));
            }
        }
        for (i = v0; i < v1; ++i) {
            /*for (j = 1; j < ny - 1; ++j) { */
            for (j = 1; j < ny; ++j) {
                V(i, j) += dt / f->p.dy * (P(i, j) - P(i, j - 1));
            }
        }
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
    double x[] = { 0, 1 };
    double y[] = { 0, 1 };
    double cfl = 0.99 / sqrt(2);        /* CFL condition: c*dt/dx = cfl <= 1/sqrt(2) */
    double T = 0.3;
    double c = 1;
    long nx = 32;
    long ny = 32;
    struct field f;
    char outfile[STR_SIZE] = "yee_pthr.tsv";
    int write = 0;
    long threads = 4;
    struct cell_partition *part;
    pthread_t *thr;
    pthread_attr_t attr;
    struct thread_param *param;

    /* Parse parameters from commandline */
    parse_cmdline(&nx, &threads, outfile, &write, argc, argv);
    ny = nx;                    /* square domain */
    printf("Domain: %li x %li\n", nx, ny);
    printf("Pthreads: %li\n", threads);

    /* Initialize */
    f = init_acoustic_field(nx, ny, x, y);
    apply_func(&f.p, gauss2d);  /* initial data */
    set_boundary(&f);

    /* Depends on the numerical variables initialized above */
    f.dt = cfl * f.p.dx / c;
    f.Nt = T / f.dt;

    /* Setup pthreads */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_barrier_init(&barrier, NULL, threads);
    /* Note that we also compute in the main threads, hence we only need to
     * allocate an additinal #(threads-1) threads. */
    thr = malloc(sizeof(pthread_t) * (threads - 1));
    if (!thr) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    /* Partition the grid */
    part = partition_grid(threads, nx);

    /* Assemble structs to be used as arguments to the threads */
    param = malloc(sizeof(struct thread_param) * threads);
    if (!param) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    long i;
    for (i = 0; i < threads; ++i) {
        param[i].tid = i;
        param[i].f = &f;
        param[i].part = part[i];
    }

    /* timestep */
    double tic, toc;
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

    /* write to disk and free data */
    if (write)
        write_to_disk(f.p, outfile);
    free(param);
    free(thr);
    free(part);
    free_acoustic_field(f);

    pthread_attr_destroy(&attr);
    pthread_barrier_destroy(&barrier);
    return EXIT_SUCCESS;
}
