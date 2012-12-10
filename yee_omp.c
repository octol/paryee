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
#include <unistd.h>
#include <math.h>

#include <omp.h>

#include "yee_common.h"

/* 
 * Main program 
 */
int main(int argc, char *argv[])
{
    /* Parameters */
    double length = 1;
    double cfl = 1;
    double T = 0.3;
    double c = 1;
    long nx = 2048;
    long n;                     /* time step index */
    double tic, toc;
    struct field f;
    long threads = 4;
    long cells_per_thread;
    struct partition *part;
    char outfile_p[STR_SIZE] = "yee_omp_p.tsv";
    char outfile_u[STR_SIZE] = "yee_omp_u.tsv";

    /* Parse parameters from commandline */
    parse_cmdline(&nx, &threads, outfile_p, outfile_u, argc, argv);
    omp_set_num_threads(threads);
    printf("Running with: N=%ld, threads=%ld\n", nx, threads);

    /* Initialize */
    f = init_acoustic_field(nx, 0, length);
#pragma omp parallel sections
    {
#pragma omp section
        {
            /*alloc_field(&f.p, nx); */
            /*set_grid(&f.p, 0.5 * length / nx, length - 0.5 * length / nx); */
            apply_func(&f.p, gauss);    /* initial data */
        }
#pragma omp section
        {
            /*alloc_field(&f.u, nx + 1); */
            /*set_grid(&f.u, 0, length); */
            apply_func(&f.u, zero);     /* initial data */
        }
    }

    /* check integrity of the generated structures */
    assert(fabs(f.p.dx - f.u.dx) < 1e-14);

    /* Depends on the numerical variables initialized above */
    f.dt = cfl * f.p.dx / c;    /* CFL condition is: c*dt/dx = cfl <= 1 */
    f.Nt = T / f.dt;

    /* Setup parallellization */
    long i;
    cells_per_thread = nx / threads;
    assert(cells_per_thread * threads == nx);
    part = malloc((long) sizeof(struct partition) * threads);
    if (!part) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < threads; ++i)
        part[i] = partition_grid(i, cells_per_thread);

    /* Maybe we can make optimization of the inner loop a bit easier for the
     * compiler? 
     * This is probably not needed. */
    double *p = f.p.value;
    double *u = f.u.value;

    /* timestep */
    tic = gettime();
#pragma omp parallel default(shared) private(i,n)
    {
        /* private variables used in the time stepping */
        double dx = f.p.dx;
        double dt = f.dt;
        double Nt = f.Nt;
        long tid = omp_get_thread_num();
        long p0 = part[tid].p[0];       /* start index */
        long p1 = part[tid].p[1];       /* end index */
        long u0 = part[tid].u[0];       /* start index */
        long u1 = part[tid].u[1];       /* end index */

        for (n = 0; n < Nt; ++n) {

            /* update the pressure (p) */
            for (i = p0; i <= p1; ++i)
                p[i] += dt / dx * (u[i + 1] - u[i]);
#pragma omp barrier

            /* update the velocity (u) */
            for (i = u0; i <= u1; ++i)
                u[i] += dt / dx * (p[i] - p[i - 1]);
#pragma omp barrier
        }
    }

    toc = gettime();
    printf("Elapsed: %f seconds\n", toc - tic);

    /* write data to disk and free data */
#pragma omp parallel sections
    {
#pragma omp section
        write_to_disk(f.p, outfile_p);
#pragma omp section
        write_to_disk(f.u, outfile_u);
    }
    free(part);
    /*free_field(f.p); */
    /*free_field(f.u); */
    free_acoustic_field(f);

    return EXIT_SUCCESS;
}
