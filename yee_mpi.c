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
 * Basic MPI distributed memory parallelizaton
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>             // sleep

#include "mpi.h"

#include "yee_common.h"

#define BEGIN 1                 /* message tag */
#define COLLECT 2               /* message tag */
#define PTAG 3                  /* communication tag */
#define UTAG 4                  /* communication tag */
#define VTAG 5                  /* communication tag */

void send_grid_data(int taskid, long left, long right, long size,
                    double y[2])
{
    /* Send out neighbour information */
    MPI_Send(&left, 1, MPI_LONG, taskid, BEGIN, MPI_COMM_WORLD);
    MPI_Send(&right, 1, MPI_LONG, taskid, BEGIN, MPI_COMM_WORLD);

    /* Send out partition information */
    MPI_Send(&size, 1, MPI_LONG, taskid, BEGIN, MPI_COMM_WORLD);

    /* Send the coordinates of the local partition (domain) */
    MPI_Send(y, 2, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);
}

void receive_grid_data(long *left, long *right, long *size, double *y,
                       MPI_Status * status)
{
    /* Data on neighbours */
    MPI_Recv(left, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, status);
    MPI_Recv(right, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, status);

    /* Partition information */
    MPI_Recv(size, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, status);

    /* Domain coordinates */
    MPI_Recv(y, 2, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, status);
}

void send_field_data(long taskid, struct field *f,
                     struct cell_partition part)
{
    /* Send out field data.
     * Note that MASTER deals with the taskid=0 special case, i.e., we do
     * not have to worry about it here. */
    long begin = part.begin;    /* for convenience */
    long size = part.size;

    long nx = f->p.size_x;
    long begin_p = 0 + begin * nx;
    long begin_u = 0 + begin * (nx + 1);
    long begin_v = 0 + begin * nx;
    long size_p = size * nx;
    long size_u = size * (nx + 1);
    long size_v = size * nx;

    MPI_Send(&f->p.value[begin_p], size_p, MPI_DOUBLE, taskid, BEGIN,
             MPI_COMM_WORLD);
    MPI_Send(&f->u.value[begin_u], size_u, MPI_DOUBLE, taskid, BEGIN,
             MPI_COMM_WORLD);
    MPI_Send(&f->v.value[begin_v], size_v, MPI_DOUBLE, taskid, BEGIN,
             MPI_COMM_WORLD);
}

void receive_field_data(struct field *f, long size, MPI_Status * status)
{
    long nx = f->p.size_x;
    long begin_p = 0 + 1 * nx;
    long begin_u = 0 + 0 * (nx + 1);
    long begin_v = 0 + 0 * nx;
    long size_p = size * nx;
    long size_u = size * (nx + 1);
    long size_v = size * nx;

    MPI_Recv(&f->p.value[begin_p], size_p, MPI_DOUBLE, MASTER, BEGIN,
             MPI_COMM_WORLD, status);
    MPI_Recv(&f->u.value[begin_u], size_u, MPI_DOUBLE, MASTER, BEGIN,
             MPI_COMM_WORLD, status);
    MPI_Recv(&f->v.value[begin_v], size_v, MPI_DOUBLE, MASTER, BEGIN,
             MPI_COMM_WORLD, status);
}

void return_data(struct field *f, long size)
{
    long nx = f->p.size_x;
    long begin_p = 0 + 1 * nx;
    long begin_u = 0 + 0 * (nx + 1);
    long begin_v = 0 + 0 * nx;
    long size_p = size * nx;
    long size_u = size * (nx + 1);
    long size_v = size * nx;

    MPI_Send(&f->p.value[begin_p], size_p, MPI_DOUBLE, MASTER, COLLECT,
             MPI_COMM_WORLD);
    MPI_Send(&f->u.value[begin_u], size_u, MPI_DOUBLE, MASTER, COLLECT,
             MPI_COMM_WORLD);
    MPI_Send(&f->v.value[begin_v], size_v, MPI_DOUBLE, MASTER, COLLECT,
             MPI_COMM_WORLD);
}


void collect_data(long taskid, struct cell_partition part, struct field *f,
                  MPI_Status * status)
{
    long begin = part.begin;    /* for convenience */
    long size = part.size;

    long nx = f->p.size_x;
    long begin_p = 0 + begin * nx;
    long begin_u = 0 + begin * (nx + 1);
    long begin_v = 0 + begin * nx;
    long size_p = size * nx;
    long size_u = size * (nx + 1);
    long size_v = size * nx;

    MPI_Recv(&f->p.value[begin_p], size_p, MPI_DOUBLE, taskid, COLLECT,
             MPI_COMM_WORLD, status);
    MPI_Recv(&f->u.value[begin_u], size_u, MPI_DOUBLE, taskid, COLLECT,
             MPI_COMM_WORLD, status);
    MPI_Recv(&f->v.value[begin_v], size_v, MPI_DOUBLE, taskid, COLLECT,
             MPI_COMM_WORLD, status);
}

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
    /* Contains field data (pressure and velocity) */
    struct field f;
    /* Array of domain partitions */
    struct cell_partition *part;

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

        /* Compute partition, neighbours and then send out data to the
         * workers */
        /* NOTE: partitioning across ny! (not nx) */
        part = partition_grid(numtasks, ny);

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

        /*
         * Time stepping section
         */

        /* Time step on the leftmost partition.  */
        /*long left = NONE; */
        long right = (numtasks > 1) ? 1 : NONE;
        long size = part[0].size;

        /* Depends on the numerical variables initialized above */
        /* CFL condition is: c*dt/dx = cfl <= 1 */
        f.dt = cfl * f.p.dx / c;
        f.Nt = T / f.dt;

        /* used by index macro */
        double *p = f.p.value;
        double *u = f.u.value;
        double *v = f.v.value;

        /* timing */
        double tic, toc;
        tic = gettime();

        long i, j;
        for (long n = 0; n < f.Nt; ++n) {

            /* Communicate */
            /* receive v from the right (top) */
            if (right != NONE) {
                long begin_v = 0 + size * nx;
                long size_v = 1 * nx;
                MPI_Recv(&f.v.value[begin_v], size_v, MPI_DOUBLE, right,
                         VTAG, MPI_COMM_WORLD, &status);
            }

            /* update the pressure (p) */
            /*printf("Master: updating p\n"); */
            for (i = 0; i < nx; ++i) {
                for (j = part[taskid].begin; j < part[taskid].end + 1; ++j) {
                    P(i, j) +=
                        f.dt / f.u.dx * (U(i + 1, j) - U(i, j)) +
                        f.dt / f.v.dy * (V(i, j + 1) - V(i, j));
                }
            }

            /* Communicate */
            /* send p to the right (top) */
            if (right != NONE) {
                long begin_p = 0 + (size - 1) * nx;     /* since we don't have extra ghost p */
                long size_p = 1 * nx;
                MPI_Send(&f.p.value[begin_p], size_p, MPI_DOUBLE, right,
                         PTAG, MPI_COMM_WORLD);
            }

            /* update the velocity (u,v) */
            for (i = 1; i < nx; ++i)
                for (j = part[taskid].begin; j < part[taskid].end + 1; ++j)
                    U(i, j) += f.dt / f.p.dx * (P(i, j) - P(i - 1, j));

            for (i = 0; i < nx; ++i)
                for (j = 1; j < part[taskid].end + 1; ++j)
                    V(i, j) += f.dt / f.p.dy * (P(i, j) - P(i, j - 1));
        }

        /*
         * End time stepping section
         */

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
        double y_part[2];
        receive_grid_data(&left, &right, &size, y_part, &status);

#ifdef DEBUG
        printf("Task %i received:  ", taskid);
        printf("left=%li  right=%li  ", left, right);
        printf("size=%li  ", size);
        printf("y[0]=%.2f  y[1]=%.2f  ", y_part[0], y_part[1]);
        printf("\n");
#endif

        /* Allocate and receive the local copy of the field */
        f = init_local_acoustic_field(nx, size, x, y_part);
        receive_field_data(&f, size, &status);

        /* Set out v to zero for the rightmost (top) partition. */
        if (right == NONE) {
            long j = size;
            long nx = f.p.size_x;
            double *v = f.v.value;
            for (long i = 0; i < nx; ++i) {
                V(i, j) = 0;
            }
        }

        /* Depends on the numerical variables initialized above */
        /* CFL condition is: c*dt/dx = cfl <= 1 */
        f.dt = cfl * f.p.dx / c;
        f.Nt = T / f.dt;

        /* used by index macro */
        long nx = f.p.size_x;
        double *p = f.p.value;
        double *u = f.u.value;
        double *v = f.v.value;

        /* Time step */
        long i, j;
        for (long n = 0; n < f.Nt; ++n) {
            /*printf("Worker %i: %li\n", taskid, n); */
            /* Communicate */
            /* send v to the left (bottom) */
            /* receive v from the right (top) */
            if (left != NONE) {
                long begin_v = 0 + 0 * nx;
                long size_v = 1 * nx;
                MPI_Send(&f.v.value[begin_v], size_v, MPI_DOUBLE, left,
                         VTAG, MPI_COMM_WORLD);
            }
            if (right != NONE) {
                long begin_v = 0 + size * nx;
                long size_v = 1 * nx;
                MPI_Recv(&f.v.value[begin_v], size_v, MPI_DOUBLE, right,
                         VTAG, MPI_COMM_WORLD, &status);
            }

            /* update the pressure (p) */
            for (i = 0; i < nx; ++i) {
                for (j = 0; j < size; ++j) {
                    P(i, j + 1) +=
                        f.dt / f.u.dx * (U(i + 1, j) - U(i, j)) +
                        f.dt / f.v.dy * (V(i, j + 1) - V(i, j));
                }
            }

            /* Communicate */
            /* receive p from the left */
            /* send p to the right */
            if (left != NONE) {
                long begin_p = 0 + 0 * nx;
                long size_p = 1 * nx;
                MPI_Recv(&f.p.value[begin_p], size_p, MPI_DOUBLE, left,
                         PTAG, MPI_COMM_WORLD, &status);
            }
            if (right != NONE) {
                long begin_p = 0 + size * nx;
                long size_p = 1 * nx;
                MPI_Send(&f.p.value[begin_p], size_p, MPI_DOUBLE, right,
                         PTAG, MPI_COMM_WORLD);
            }

            /* update the velocity (u) */
            for (i = 1; i < nx; ++i)
                for (j = 0; j < size; ++j)
                    U(i, j) +=
                        f.dt / f.p.dx * (P(i, j + 1) - P(i - 1, j + 1));

            /*printf("worker %i: updating v\n",taskid); */
            for (i = 0; i < nx; ++i)
                for (j = 0; j < size; ++j)
                    V(i, j) += f.dt / f.p.dy * (P(i, j + 1) - P(i, j));
        }

        /* Send back data to master */
        return_data(&f, size);

        /* Remember to deallocate */
        free_acoustic_field(f);

        MPI_Finalize();
    }

    return 0;
}
