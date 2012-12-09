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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mpi.h"

#include "yee_common.h"

#define BEGIN 1                 /* message tag */
#define COLLECT 2               /* message tag */
#define UTAG 3                  /* communication tag */
#define PTAG 4                  /* communication tag */

int main(int argc, char *argv[])
{
    /* Simulation parameters */
    double length = 1;
    double cfl = 1;
    double T = 0.3;
    double c = 1;
    char outfile_p[STR_SIZE] = "yee_mpi_p.tsv";
    char outfile_u[STR_SIZE] = "yee_mpi_u.tsv";

    /* Numerical paramters */
    long nx = 2048;

    long cells_per_worker = 0;
    long left = 0;
    long right = 0;

    double tic, toc;

    /* Parse parameters from commandline */
    parse_cmdline(&nx, NULL, outfile_p, outfile_u, argc, argv);

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
    int numworkers = numtasks - 1;

    /* Field variables - data structures for the simulation */
    /* Contains field data (pressure and velocity) */
    struct field f;
    /* Array of domain partitions */
    struct partition *part;

    if (taskid == MASTER) {
        /********************* Master code *********************/
        if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)) {
            printf("ERROR: number of tasks must be between %d and %d.\n",
                   MINWORKER + 1, MAXWORKER + 1);
            printf("Quitting...\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(EXIT_FAILURE);
        }

        /*printf ("Starting yee_mpi with %i workers\n",numworkers); */
        printf("Running with: N=%ld, workers=%i\n", nx, numworkers);

        /* Initialize */
        alloc_field(&f.p, nx);
        alloc_field(&f.u, nx + 1);
        set_grid(&f.p, 0.5 * length / nx, length - 0.5 * length / nx);
        set_grid(&f.u, 0, length);
        apply_func(&f.p, gauss);        /* initial data */
        apply_func(&f.u, zero); /* initial data */

        /* Compute partition, neighbours and then send out data to the
         * workers */
        cells_per_worker = nx / numworkers;
        if (cells_per_worker * numworkers != nx) {
            fprintf(stderr, "Only 2^n+1, n=1,2,..., threads supported\n");
            exit(EXIT_FAILURE);
        }
        part = malloc(sizeof(struct partition) * numworkers);
        if (!part) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }

        /* remember that tid=0 is the master */
        for (int tid = 1; tid <= numworkers; ++tid) {
            /* partition_grid() requires the nodes to be enumerated starting
             * from 0: 0,...,total_nodes.
             * Note that partitition_grid2 returns the inner nodes only */
            part[tid - 1] = partition_grid(tid - 1, cells_per_worker);
            /* compute neighbours */
            left = (tid == 1) ? NONE : tid - 1;
            right = (tid == numworkers) ? NONE : tid + 1;

            /* Check that everything is ok before we send out */
            verify_grid_integrity(part[tid - 1], tid, nx, numworkers,
                                  left);

            /* Make field data more compact, so we can more easily use as 
             * array indices. For convenience only.
             * Note, these are the inner nodes. */
            long bp, ep, sp, bu, eu, su;
            expand_indices(part[tid - 1], &bp, &ep, &sp, &bu, &eu, &su);

            /* timing */
            tic = gettime();

            /* Send out neighbour information */
            MPI_Send(&left, 1, MPI_LONG, tid, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&right, 1, MPI_LONG, tid, BEGIN, MPI_COMM_WORLD);

            /* Send out partition information */
            MPI_Send(&bp, 1, MPI_LONG, tid, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&ep, 1, MPI_LONG, tid, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&sp, 1, MPI_LONG, tid, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&bu, 1, MPI_LONG, tid, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&eu, 1, MPI_LONG, tid, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&su, 1, MPI_LONG, tid, BEGIN, MPI_COMM_WORLD);

            /* Send out field data */
            MPI_Send(&f.p.value[bp], sp, MPI_DOUBLE, tid, BEGIN,
                     MPI_COMM_WORLD);
            MPI_Send(&f.u.value[bu], su, MPI_DOUBLE, tid, BEGIN,
                     MPI_COMM_WORLD);

            /* Send out grid data */
            MPI_Send(&f.p.x[bp], sp, MPI_DOUBLE, tid, BEGIN,
                     MPI_COMM_WORLD);
            MPI_Send(&f.u.x[bu], su, MPI_DOUBLE, tid, BEGIN,
                     MPI_COMM_WORLD);
            MPI_Send(&f.p.dx, 1, MPI_DOUBLE, tid, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&f.u.dx, 1, MPI_DOUBLE, tid, BEGIN, MPI_COMM_WORLD);

#ifdef DEBUG
            printf("Sent to task %d: left=%ld, right=%ld\n", tid, left,
                   right);
#endif
        }

        /* Wait for returning data from the workers */
        for (int tid = 1; tid <= numworkers; ++tid) {
            /* Make field data more compact, so we can more easily use as 
             * array indices. For convenience only. */
            long bp, ep, sp, bu, eu, su;
            expand_indices(part[tid - 1], &bp, &ep, &sp, &bu, &eu, &su);

            MPI_Recv(&f.p.value[bp], sp, MPI_DOUBLE, tid, COLLECT,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&f.u.value[bu], su, MPI_DOUBLE, tid, COLLECT,
                     MPI_COMM_WORLD, &status);

#ifdef DEBUG
            printf("Results collected from task %d\n", tid);
#endif
        }

        toc = gettime();
        printf("Elapsed: %f seconds\n", toc - tic);

        write_to_disk(f.p, outfile_p);
        write_to_disk(f.u, outfile_u);

        free(part);
        free_field(f.p);
        free_field(f.u);

        MPI_Finalize();

    } else {
        /********************* Worker code *********************/
        /* array indices and sizes */
        long bp, bu, ep, eu, sp, su;
        /* local array indices and sizes */
        long lbp, lbu, lep, leu, lsp, lsu;

        part = malloc(sizeof(struct partition));
        if (!part) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }

        /* Receive grid and node parameters from master */

        /* Data on neighbours */
        MPI_Recv(&left, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD,
                 &status);
        MPI_Recv(&right, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD,
                 &status);
        /* Partition data */
        MPI_Recv(&bp, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&ep, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&sp, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&bu, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&eu, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv(&su, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, &status);

        /* Given the partition data we compute the corresponding local array
         * indices in the array that is expanded to also contain the
         * neighbouring points needed to update. */
        set_local_index(sp, su, left, &lbp, &lep, &lsp, &lbu, &leu, &lsu);

#ifdef DEBUG
        printf("taskid=%i:", taskid);
        printf("bp=%li  ep=%li  bu=%li  eu=%li  ", bp, ep, bu, eu);
        printf("sp=%li  su=%li  ", sp, su);
        printf("lbp=%li  lep=%li  lbu=%li  leu=%li  ", lbp, lep, lbu, leu);
        printf("lsp=%li  lsu=%li\n", lsp, lsu);
#endif

        /* Now that we have the grid info we can setup the data structures
         * containing the field. */
        alloc_field(&f.p, lsp);
        alloc_field(&f.u, lsu);
        apply_func(&f.p, zero); /* set everything to zero */
        apply_func(&f.u, zero); /* set everything to zero */

        /* Field data */
        /* NOTE: we only recieve sp/su number of values, not lsp/lsu. */
        MPI_Recv(&f.p.value[lbp], sp, MPI_DOUBLE, MASTER, BEGIN,
                 MPI_COMM_WORLD, &status);
        MPI_Recv(&f.u.value[lbu], su, MPI_DOUBLE, MASTER, BEGIN,
                 MPI_COMM_WORLD, &status);
        /* Grid data */
        MPI_Recv(&f.p.x[lbp], sp, MPI_DOUBLE, MASTER, BEGIN,
                 MPI_COMM_WORLD, &status);
        MPI_Recv(&f.u.x[lbu], su, MPI_DOUBLE, MASTER, BEGIN,
                 MPI_COMM_WORLD, &status);
        MPI_Recv(&f.p.dx, 1, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD,
                 &status);
        MPI_Recv(&f.u.dx, 1, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD,
                 &status);
#ifdef DEBUG
        printf("Task=%d received:  left=%ld  right=%ld", taskid, left,
               right);
        printf("  bp=%ld  ep=%ld", bp, ep);
        printf("  bu=%ld  eu=%ld\n", bu, eu);
#endif

        /* Depends on the numerical variables initialized above */
        /* CFL condition is: c*dt/dx = cfl <= 1 */
        f.dt = cfl * f.p.dx / c;
        f.Nt = T / f.dt;

        long i;
        for (long n = 0; n < f.Nt; ++n) {
            /* Communicate */
            /* send u to the left */
            /* receive u from the right */
            if (left != NONE) {
                MPI_Send(&f.u.value[lbu], 1, MPI_DOUBLE, left, UTAG,
                         MPI_COMM_WORLD);
            }
            if (right != NONE) {
                MPI_Recv(&f.u.value[leu + 1], 1, MPI_DOUBLE, right, UTAG,
                         MPI_COMM_WORLD, &status);
            }

            /* update the pressure (p) */
            /*update_field_s(&f.p, lbp, sp, &f.u, 0, f.dt); */
            for (i = lbp; i <= lbp + sp - 1; ++i)
                f.p.value[i] += f.dt / f.p.dx
                    * (f.u.value[i + 1 - lbp] - f.u.value[i - lbp]);

            /* Communicate */
            /* receive p from the left */
            /* send p to the right */
            if (left != NONE) {
                MPI_Recv(&f.p.value[lbp - 1], 1, MPI_DOUBLE, left, PTAG,
                         MPI_COMM_WORLD, &status);
            }
            if (right != NONE) {
                MPI_Send(&f.p.value[lep], 1, MPI_DOUBLE, right, PTAG,
                         MPI_COMM_WORLD);
            }

            /* update the velocity (u) */
            /*update_field_s(&f.u, lbu, su, &f.p, 0, f.dt); */
            for (i = lbu; i <= lbu + su - 1; ++i)
                f.u.value[i] += f.dt / f.p.dx
                    * (f.p.value[i - lbu + 1] - f.p.value[i - lbu]);
        }

        /* Send back data to master */
        MPI_Send(&f.p.value[lbp], sp, MPI_DOUBLE, MASTER, COLLECT,
                 MPI_COMM_WORLD);
        MPI_Send(&f.u.value[lbu], su, MPI_DOUBLE, MASTER, COLLECT,
                 MPI_COMM_WORLD);

        /* Remember to deallocate */
        free(part);
        free_field(f.p);
        free_field(f.u);

        MPI_Finalize();
    }

    return 0;
}
