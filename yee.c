#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "yee_common.h"

/* 
 * Main program 
 */
int main (int argc, char* argv[])
{
    /* Parameters */
    double length = 1; 
    double cfl = 1; 
    double T = 0.3;
    double c = 1;
    unsigned long nx = 2048;
    unsigned long n;
    double tic, toc;
    struct Field f;
    char outfile_p[STR_SIZE] = "yee_p.tsv";
    char outfile_u[STR_SIZE] = "yee_u.tsv";

    /* Parse parameters from commandline */
    parse_cmdline (&nx, NULL, outfile_p, outfile_u, argc, argv);
    printf("Running with: N=%ld\n", nx);

    /* Initialize */
    alloc_field(&f.p, nx);
    alloc_field(&f.u, nx+1);
    set_grid (&f.p, 0.5*length/nx, length-0.5*length/nx);
    set_grid (&f.u, 0, length);
    apply_func (&f.p, gauss); /* initial data */
    apply_func (&f.u, zero);  /* initial data */

    /* Depends on the numerical variables initialized above */
    f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
    f.Nt = T/f.dt;

    tic = gettime();
    for (n=0; n<f.Nt; ++n) {

        /* update the pressure (p) */
        update_field_i (&f.p, 0, nx, &f.u, 0, f.dt);

        /* update the velocity (u) */
        update_field_i (&f.u, 1, nx, &f.p, 0, f.dt);
    }
    toc = gettime();
    printf ("Elapsed: %f seconds\n", toc-tic);

    /* write data to disk and free data */
    write_to_disk(f.p, outfile_p); 
    write_to_disk(f.u, outfile_u); 
    free_field(f.p);
    free_field(f.u);

    return EXIT_SUCCESS;
}
