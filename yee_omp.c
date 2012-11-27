#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include <unistd.h>

#include "yee_common.h"

/* 
 * Main program 
 */
int main (int argc, char* argv[])
{
    /* Paramters */
    double length=1, cfl=1, T=0.3, c=1;
    unsigned long nx=2048, nodes=2;
    unsigned int n, i, cells_per_node;
    int tid;
    double tic, toc;
    struct Field f;
    struct Partition* part;

    /* Parse parameters from commandline */
    parse_cmdline (&nx, &nodes, argc, argv);
    omp_set_num_threads (nodes);
    printf("Running with: N=%ld, threads=%ld\n", nx, nodes);

    /* Initialize */
#pragma omp parallel sections 
{
#pragma omp section 
{
    alloc_field(&f.p, nx);
    set_grid (&f.p, 0.5*length/nx, length-0.5*length/nx);
    apply_func (&f.p, gauss); /* initial data */
}
#pragma omp section 
{
    alloc_field(&f.u, nx+1);
    set_grid (&f.u, 0, length);
    apply_func (&f.u, zero);  /* initial data */
}
}

    /* Depends on the numerical variables initialized above */
    f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
    f.Nt = T/f.dt;

    /* Setup parallellization */
    cells_per_node = nx/nodes;
    assert(cells_per_node*nodes == nx);
    part = malloc (sizeof(struct Partition)*nodes);
    if (!part) {
        fprintf(stderr,"Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (i=0; i<nodes; ++i)
        part[i] = partition_grid (i, nodes, cells_per_node);

    tic = gettime();

#pragma omp parallel default(shared) private(tid,n)
{
    tid = omp_get_thread_num ();
    for (n=0; n<f.Nt; ++n) {

        /* update the pressure (p) */
        /*update_field3 (&f.p, part[tid].start_p, part[tid].end_p,*/
                       /*&f.u, part[tid].start_u-1, f.dt);*/
        update_field_i (&f.p, part[tid].p[0], part[tid].p[1],
                        &f.u, part[tid].u[0]-1, f.dt);
#pragma omp barrier

        /* update the velocity (u) */
        /*update_field (&f.u, part[tid].start_u, part[tid].end_u,*/
                      /*&f.p, part[tid].start_p, f.dt);*/
        update_field_i (&f.u, part[tid].u[0], part[tid].u[1],
                        &f.p, part[tid].p[0], f.dt);
#pragma omp barrier
    }
}

    toc = gettime();
    printf ("Elapsed: %f seconds\n", toc-tic);

    /* write data to disk and free data */
#pragma omp parallel sections 
{
#pragma omp section
    write_to_disk(f.p, "yee_omp_p"); 
#pragma omp section
    write_to_disk(f.u, "yee_omp_u"); 
}
    free (part);
    free_field (f.p);
    free_field (f.u);

    return EXIT_SUCCESS;
}
