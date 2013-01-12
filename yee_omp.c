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
 * OpenMP shared memory parallelization
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <omp.h>

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
    long nx = 32;
    long ny = 32;
    struct field f;
    char outfile[STR_SIZE] = "yee_omp.tsv";
    int write = 1;
    long threads = 4;
    struct cell_partition *part;

    /* Parse parameters from commandline */
    parse_cmdline(&nx, &threads, outfile, &write, argc, argv);
    ny = nx;                    /* square domain */
    omp_set_num_threads(threads);
    printf("Domain: %li x %li\n", nx, ny);
    printf("OpenMP threads: %li\n", threads);

    /* Initialize */
    f = init_acoustic_field(nx, ny, x, y);
    apply_func(&f.p, gauss2d);  /* initial data */
    set_boundary(&f);

    /* Depends on the numerical variables initialized above */
    f.dt = cfl * f.p.dx / c;
    f.Nt = T / f.dt;

    /* Partition the grid */
    part = partition_grid(threads, nx);

    /* Maybe we can make optimization of the inner loop a bit easier for the
     * compiler? 
     * This is probably not needed. */
    double *p = f.p.value;
    double *u = f.u.value;
    double *v = f.v.value;

    /* timestep */
    long n, i, j;
    double tic, toc;
    tic = gettime();
#pragma omp parallel default(shared) private(i,j,n)
    {
        /* private variables used in the time stepping */
        double dt = f.dt;
        double Nt = f.Nt;
        long tid = omp_get_thread_num();
        long p0, p1, u0, u1, v0, v1;
        cellindex_to_nodeindex(tid, part[tid], &p0, &p1, &u0, &u1, &v0,
                               &v1);
#ifdef DEBUG
        printf("tid=%lu  p0=%lu  p1=%lu  u0=%lu  u1=%lu  v0=%lu  v1=%lu\n",
               tid, p0, p1, u0, u1, v0, v1);
#endif

        for (n = 0; n < Nt; ++n) {

            /* update the pressure (p) */
            for (i = p0; i < p1; ++i) {
                for (j = 0; j < ny; ++j) {
                    P(i, j) +=
                        dt / f.u.dx * (U(i + 1, j) - U(i, j)) +
                        dt / f.v.dy * (V(i, j + 1) - V(i, j));
                }
            }
#pragma omp barrier

            /* update the velocity (u,v) */
            for (i = u0; i < u1; ++i) {
                for (j = 0; j < ny; ++j) {
                    U(i, j) += dt / f.p.dx * (P(i, j) - P(i - 1, j));
                }
            }

            for (i = v0; i < v1; ++i)
                /*for (j = 1; j < ny - 1; ++j) */
                for (j = 1; j < ny; ++j)
                    V(i, j) += dt / f.p.dy * (P(i, j) - P(i, j - 1));
#pragma omp barrier
        }
    }

    toc = gettime();
    printf("Elapsed: %f seconds\n", toc - tic);

    /* write to disk and free data */
    if (write)
        write_to_disk(f.p, outfile);
    free(part);
    free_acoustic_field(f);

    return EXIT_SUCCESS;
}
