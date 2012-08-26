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
 * Leapfrog time update adapted to use return value and argument suitable
 * for pthreads.
 */
void* thread_update (void* arg)
{
    UpdateParam* p = (UpdateParam*) arg;
    update_field (p->dst, p->dst1, p->dst2, p->src, p->src1, p->dt);
    pthread_exit(NULL);
}

/* 
 * Main program 
 */
int main(int argc, char* argv[])
{
    double length=1, cfl=1, T=0.3, c=1;
    unsigned long nx=2048, nodes=2;
    int cells_per_node;
    int n, i;
    double tic, toc;
    Field f;
    FieldPartition part;
    UpdateParam* param_p;
    UpdateParam* param_u;
    pthread_t* thr;
    pthread_attr_t attr;

    /* Parse parameters */
    parse_cmdline (&nx, &nodes, argc, argv);
    printf("Running with: N=%ld, threads=%ld\n", nx, nodes);

    /* Initialize */
    alloc_field(&f.p, nx);
    alloc_field(&f.u, nx+1);
    set_grid (&f.p, 0.5*length/nx, length-0.5*length/nx);
    set_grid (&f.u, 0, length);
    field_func (&f.p, gauss);
    field_func (&f.u, zero);

    /* Depends on the numerical variables initialzed above */
    f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
    f.Nt = T/f.dt;

    /* Setup thread */
    pthread_attr_init (&attr);
    pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
    thr = malloc (sizeof(pthread_t)*nodes);
    if (!thr) {
        fprintf(stderr,"Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    /* Partition the grid */
    param_p = malloc (sizeof(UpdateParam)*nodes);
    param_u = malloc (sizeof(UpdateParam)*nodes);
    if (!param_p || !param_u) {
        fprintf(stderr,"Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    cells_per_node = nx/nodes;
    assert(cells_per_node*nodes == nx);
    for (i=0; i<nodes; ++i) {
        part = partition_grid (i, nodes, cells_per_node);
        param_p[i] = collect_param (&f.p, part.start_p, part.end_p, 
                                    &f.u, part.start_u-1, f.dt);
        param_u[i] = collect_param (&f.u, part.start_u, part.end_u, 
                                    &f.p, part.start_p, f.dt);
    }

    /* Timestep */
    tic = gettime ();
    for (n=0; n<f.Nt; ++n) {
        /* update the pressure (p) */
        for (i=0; i<nodes; ++i) 
            pthread_create (&thr[i], &attr, thread_update, &param_p[i]);
        for (i=0; i<nodes; ++i)
            pthread_join (thr[i], NULL);

        /* update the velocity (u) */
        for (i=0; i<nodes; ++i)
            pthread_create (&thr[i], &attr, thread_update, &param_u[i]);
        for (i=0; i<nodes; ++i)
            pthread_join (thr[i], NULL);
    }
    toc = gettime ();
    printf ("Elapsed: %f seconds\n", toc-tic);

    /* write data to disk and free data */
    write_to_disk (f.p, "output_p"); 
    write_to_disk (f.u, "output_u"); 
    free (thr);
    free (param_p);
    free (param_u);
    free_field (f.p);
    free_field (f.u);

    pthread_attr_destroy (&attr);
    return EXIT_SUCCESS;
}
