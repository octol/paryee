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
 * Basic single threaded implementation to compare parallelization overhead.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "yee_common.h"

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
    unsigned long nx = 32;
    unsigned long ny = 32;
    double tic, toc;
    struct field f;
    char outfile[STR_SIZE] = "yee.tsv";

    /* Parse parameters from commandline */
    parse_cmdline(&nx, NULL, outfile, argc, argv);
    ny = nx;                    /* square domain */
    printf("Domain: %li x %li\n", nx, ny);

    /* Initialize */
    f = init_acoustic_field(nx, ny, x, y);
    apply_func(&f.p, gauss2d);  /* initial data */

    /* Depends on the numerical variables initialized above */
    f.dt = cfl * f.p.dx / c;
    f.Nt = T / f.dt;

    /* timestep */
    tic = gettime();
    timestep_leapfrog(&f, f.Nt);
    toc = gettime();
    printf("Elapsed: %f seconds\n", toc - tic);

    /* write to disk and free data */
    write_to_disk(f.p, outfile);
    free_acoustic_field(f);

    return EXIT_SUCCESS;
}
