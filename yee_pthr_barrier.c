/* 
 * Basic implementation of 1D wave equation on system form
 *
 *      p_t = a u_x
 *      u_t = b p_x
 *
 * This formulation is sometimes used in acoustics, where then p is the
 * pressure and u is the velocity field in the direction of the x-axis.
 * 
 * We use a staggared grid with u defined on the endpoints, and thus on full
 * integer indices. N denotes the number of interval lengths, i.e., we have
 * N+1 number of nodes for u
 *
 * The outer (x=0, x=L) boundary is set to u=0. This means we start the
 * enumeration in the second u point when we partition the domain. 
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
typedef struct ThreadParam ThreadParam;

struct ThreadParam {
    UpdateParam p; /* to update p */
    UpdateParam u; /* to update u */
    int Nt;
};

void* thread_main (void* arg)
{
    ThreadParam* param = (ThreadParam*) arg;
    UpdateParam p = param->p;
    UpdateParam u = param->u;
    int n;
    
    for (n=0; n<param->Nt; ++n) {
        /* update the pressure (p) */
        update_field (p.dst, p.dst1, p.dst2, p.src, p.src1, p.dt);
        pthread_barrier_wait (&barrier);
        /* update the velocity (u) */
        update_field (u.dst, u.dst1, u.dst2, u.src, u.src1, u.dt);
        pthread_barrier_wait (&barrier);
    }

    return NULL;
}

/* 
 * Main program 
 */
int main (int argc, char* argv[])
{
    double length=1, cfl=1, T=0.3, c=1;
    unsigned long nx=2048, nodes=2;
    int i, cells_per_node;
    double tic, toc;
    Field f;
    FieldPartition part;
    pthread_t* thr;
    pthread_attr_t attr;
    ThreadParam* param;

    /* Parse parameters from commandline */
    parse_cmdline (&nx, &nodes, argc, argv);
    printf("Running with: N=%ld, threads=%ld\n", nx, nodes);

    /* Initialize */
    alloc_field(&f.p, nx);
    alloc_field(&f.u, nx+1);
    set_grid (&f.p, 0.5*length/nx, length-0.5*length/nx);
    set_grid (&f.u, 0, length);
    field_func (&f.p, gauss); /* initial data */
    field_func (&f.u, zero);  /* initial data */

    /* Depends on the numerical variables initialized above */
    f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
    f.Nt = T/f.dt;

    /* Setup threads */
    pthread_attr_init (&attr);
    pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
    pthread_barrier_init (&barrier, NULL, nodes);
    thr = malloc (sizeof(pthread_t)*(nodes-1));
    if (!thr) {
        fprintf(stderr,"Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    /* Partion the grid and assemble structs to be used as arguments to the
     * threads */
    param = malloc (sizeof(ThreadParam)*nodes);
    if (!param) {
        fprintf(stderr,"Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    cells_per_node = nx/nodes;
    assert(cells_per_node*nodes == nx);
    for (i=0; i<nodes; ++i) {
        part = partition_grid (i, nodes, cells_per_node);
        param[i].p = collect_param (&f.p, part.start_p, part.end_p, 
                                    &f.u, part.start_u-1, f.dt);
        param[i].u = collect_param (&f.u, part.start_u, part.end_u, 
                                    &f.p, part.start_p, f.dt);
        param[i].Nt = f.Nt;
    }
    
    /* timing */
    tic = gettime ();

    /* Spawn additional [nodes-1] threads */
    for (i=0; i<nodes-1; ++i)
        pthread_create (&thr[i], &attr, thread_main, &param[i]);
    /* Current thread */
    thread_main (&param[nodes-1]);

    /* Collect threads when finished */
    for (i=0; i<nodes-1; ++i)
        pthread_join (thr[i], NULL);

    toc = gettime ();
    printf ("Elapsed: %f seconds\n", toc-tic);

    /* write data to disk and free data */
    write_to_disk(f.p, "yee_pthr_barrier_p"); 
    write_to_disk(f.u, "yee_pthr_barrier_u"); 
    free (param);
    free (thr);
    free_field (f.p);
    free_field (f.u);

    pthread_attr_destroy (&attr);
    pthread_barrier_destroy (&barrier);
    return EXIT_SUCCESS;
}
