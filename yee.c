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
 * The outer (x=0, x=L) boundary is set to u=0. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "yee_common.h"

/* 
 * Main program 
 */
int main (int argc, char* argv[])
{
    /* Paramters */
    double length=1, cfl=1, T=0.3, c=1;
    unsigned long nx=2048, n;
    double tic, toc;
    Field f;

    /* Parse parameters from commandline */
    parse_cmdline (&nx, NULL, argc, argv);
    printf("Running with: N=%ld\n", nx);

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

    tic = gettime();
    for (n=0; n<f.Nt; ++n) {

        /* update the pressure (p) */
        update_field (&f.p, 0, nx, &f.u, 0, f.dt);

        /* update the velocity (u) */
        update_field (&f.u, 1, nx, &f.p, 0, f.dt);
    }
    toc = gettime();
    printf ("Elapsed: %f seconds\n", toc-tic);

    /* write data to disk and free data */
    write_to_disk(f.p, "output_p"); 
    write_to_disk(f.u, "output_u"); 
    free_field(f.p);
    free_field(f.u);

    return EXIT_SUCCESS;
}
