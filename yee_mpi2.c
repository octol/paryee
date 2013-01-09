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
 * Basic MPI distributed memory parallelizaton.
 * This version tries to be non-blocking.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "mpi.h"

#include "yee_common.h"
#include "yee_common_mpi.h"

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
    char outfile[STR_SIZE] = "yee_mpi.tsv";

    /* Parse parameters from commandline */
    parse_cmdline(&nx, NULL, outfile, argc, argv);
    ny = nx;                    /* square domain */

    /* Setup MPI and initialize variables */
    MPI_Status status;
    int rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    int numtasks;
    int taskid;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    /* Field variables - data structures for the simulation */
    struct field f;             /* Contains pressure and velocity */
    struct cell_partition *part;        /* Array of domain partitions */

    if (taskid == MASTER) {
        /********************* Master code *********************/
        /* MASTER deals with the setting up, sending out, receiving the result,
         * writing to disk. 
         * It also deals with the bottom partition (first partition). */

        /* Initialize */
        printf("Domain: %li x %li\n", nx, ny);
        printf("MPI processes: %d\n", numtasks);
        f = init_acoustic_field(nx, ny, x, y);
        apply_func(&f.p, gauss2d);      /* initial data */
        set_boundary(&f);

        /* Depends on the numerical variables initialized above */
        /* CFL condition is: c*dt/dx = cfl <= 1 */
        f.dt = cfl * f.p.dx / c;
        f.Nt = T / f.dt;

        /* Compute partition. NOTE: partitioning across ny! (not nx) */
        part = partition_grid(numtasks, ny);

        /* Timing */
        double tic, toc;
        tic = gettime();

        /* Send out data to the workers */
        for (long taskid = 1; taskid < numtasks; ++taskid) {
            /* Compute neighbours */
            long left = taskid - 1;
            long right = (taskid == numtasks - 1) ? NONE : taskid + 1;

            /* Compute (start/end) coordinates of the partition */
            double y_part[2];
            get_partition_coords(part[taskid], &f, y_part);

            /* Send out data to the other workers / computational nodes */
            /* We send the number of inner cells (size) */
            send_grid_data(taskid, left, right, part[taskid].size, y_part);
            send_field_data(taskid, &f, part[taskid]);
        }

        /* The neighbours of the master node */
        long left = NONE;
        long right = (numtasks > 1) ? 1 : NONE;

        /* Time step on the leftmost partition.  */
        timestep_mpi2(&f, left, right, part[0].size);

        /* Wait for returning data from the workers */
        for (long taskid = 1; taskid < numtasks; ++taskid)
            collect_data(taskid, part[taskid], &f, &status);

        toc = gettime();
        printf("Elapsed: %f seconds\n", toc - tic);

        write_to_disk(f.p, outfile);
        free(part);
        free_acoustic_field(f);
        MPI_Finalize();

    } else {
        /********************* Worker code *********************/
        /* Receives grid and node parameters from master */

        /* Partition data */
        long left, right, size;
        double y_part[2];       /* coord interval on y-axis */
        receive_grid_data(&left, &right, &size, y_part, &status);

        /* Allocate and receive the local copy of the field */
        f = init_local_acoustic_field(nx, size, x, y_part);
        receive_field_data(&f, size, &status);

        /* Depends on the numerical variables initialized above */
        /* CFL condition is: c*dt/dx = cfl <= 1 */
        f.dt = cfl * f.p.dx / c;
        f.Nt = T / f.dt;

        /* Set outer v to zero for the rightmost (top) partition. */
        if (right == NONE) {
            long j = size;
            for (long i = 0; i < nx; ++i)
                assign_to(f.v, i, j, 0);
        }

        timestep_mpi2(&f, left, right, size);
        return_data(&f, size);  /* Send back data to master */
        free_acoustic_field(f);
        MPI_Finalize();
    }

    return 0;
}
