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
struct ThreadParam {
    struct UpdateParam p;       /* to update p */
    struct UpdateParam u;       /* to update u */
    double Nt;
};

void *thread_main(void *arg)
{
    struct ThreadParam *param = (struct ThreadParam *) arg;
    struct UpdateParam p = param->p;
    struct UpdateParam u = param->u;
    int n;

    for (n = 0; n < param->Nt; ++n) {
        /* update the pressure (p) */
        update_field_i(p.dst, p.dst1, p.dst2, p.src, p.src1, p.dt);
        pthread_barrier_wait(&barrier);
        /* update the velocity (u) */
        update_field_i(u.dst, u.dst1, u.dst2, u.src, u.src1, u.dt);
        pthread_barrier_wait(&barrier);
    }

    return NULL;
}

/* 
 * Main program 
 */
int main(int argc, char *argv[])
{
    double length = 1, cfl = 1, T = 0.3, c = 1;
    unsigned long nx = 2048, nodes = 4;
    unsigned int i, cells_per_node;
    double tic, toc;
    struct Field f;
    struct Partition part;
    pthread_t *thr;
    pthread_attr_t attr;
    struct ThreadParam *param;
    char outfile_p[STR_SIZE] = "yee_pthr_p.tsv";
    char outfile_u[STR_SIZE] = "yee_pthr_u.tsv";

    /* Parse parameters from commandline */
    parse_cmdline(&nx, &nodes, outfile_p, outfile_u, argc, argv);
    printf("Running with: N=%ld, threads=%ld\n", nx, nodes);

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
    pthread_barrier_init(&barrier, NULL, nodes);
    thr = malloc(sizeof(pthread_t) * (nodes - 1));
    if (!thr) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    /* Partion the grid and assemble structs to be used as arguments to the
     * threads */
    param = malloc(sizeof(struct ThreadParam) * nodes);
    if (!param) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    cells_per_node = nx / nodes;
    assert(cells_per_node * nodes == nx);
    for (i = 0; i < nodes; ++i) {
        part = partition_grid(i, cells_per_node);
        param[i].p = collect_param(&f.p, part.p[0], part.p[1],
                                   &f.u, part.u[0] - 1, f.dt);
        param[i].u = collect_param(&f.u, part.u[0], part.u[1],
                                   &f.p, part.p[0], f.dt);
        param[i].Nt = f.Nt;
    }

    /* timing */
    tic = gettime();

    /* Spawn additional [nodes-1] threads */
    for (i = 0; i < nodes - 1; ++i)
        pthread_create(&thr[i], &attr, thread_main, &param[i]);
    /* Current thread */
    thread_main(&param[nodes - 1]);

    /* Collect threads when finished */
    for (i = 0; i < nodes - 1; ++i)
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
