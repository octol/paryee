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
    double length = 1;
    double cfl = 0.99 / sqrt(2);
    double T = 0.3;
    double c = 1;
    long nx = 8;
    long n;                     /* time step index */
    double tic, toc;
    struct field f;
    char outfile[STR_SIZE] = "yee.tsv";

    /* Parse parameters from commandline */
    /*parse_cmdline(&nx, NULL, outfile, argc, argv); */
    /*printf("Running with: N=%li\n", nx); */

    /* Initialize */
    /*f = init_acoustic_field(nx, ny, 0, length);                           */
    /*apply_func(&f.p, gauss);    [> initial data <]                        */
    /*apply_func(&f.u, zero);     [> initial data <]                        */

    /*[> Depends on the numerical variables initialized above <]            */
    /*double dx = f.p.dx;                                                   */
    /*double dt = cfl * dx / c;   [> CFL condition is: c*dt/dx = cfl <= 1 <] */
    /*double Nt = T / dt;                                                   */
    /*long i;                                                               */

    /*[> timestep <]                                                        */
    /*tic = gettime();                                                      */
    /*for (n = 0; n < Nt; ++n) {                                            */
    /*    [> update p <]                                                    */
    /*    for (i = 0; i < nx; ++i)                                          */
    /*        f.p.value[i] += dt / dx * (f.u.value[i + 1] - f.u.value[i]);  */

    /*    [> update u <]                                                    */
    /*    for (i = 1; i < nx - 1; ++i)                                      */
    /*        f.u.value[i] += dt / dx * (f.p.value[i] - f.p.value[i - 1]);  */
    /*}                                                                     */
    /*toc = gettime();                                                      */
    /*printf("Elapsed: %f seconds\n", toc - tic);                           */

    /*[> write data to disk and free data <]                                */
    /*write_field_to_disk(f, outfile_p, outfile_u);                         */
    /*free_acoustic_field(f);                                               */

    return EXIT_SUCCESS;
}
