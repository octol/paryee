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

#include "yee_common_mpi.h"

/*
 * Sends the grid data from the master node to the worker nodes.
 */
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

/*
 * Each worker node receives its local grid data from the master node.
 */
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

/*
 * The actual field (pressure and velocity) is sent from master to worker node. 
 */
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

/*
 * Each worker node receives its field (pressure and velocity) from the master
 * node.
 */
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

/* 
 * After simulation/time stepping is done, each worker node sends its field
 * data back to the master node.
 */
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

/*
 * The master collects all data into one unit.
 */
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
 * Compute one time step for the pressure. This is the same expression as in
 * the shared memory implementations.
 */
void leapfrog_master_p(double *restrict p, double *restrict u,
                       double *restrict v, long nx, double dt, double dx,
                       double dy, long size)
{
    for (long i = 0; i < nx; ++i)
        for (long j = 0; j < size; ++j)
            P(i, j) +=
                dt / dx * (U(i + 1, j) - U(i, j)) +
                dt / dy * (V(i, j + 1) - V(i, j));
}

/*
 * Compute one time step for the pressure on the worker nodes. This differs
 * since we now have extra ghost cells to take into consideration. In
 * particular, grid indices for pressure needs to be shifter up by one along
 * the y-axis.
 */
void leapfrog_worker_p(double *restrict p, double *restrict u,
                       double *restrict v, long nx, double dt, double dx,
                       double dy, long size)
{
    for (long i = 0; i < nx; ++i)
        for (long j = 0; j < size; ++j)
            P(i, j + 1) +=
                dt / dx * (U(i + 1, j) - U(i, j)) +
                dt / dy * (V(i, j + 1) - V(i, j));
}

/* 
 * Compute one time step for the velocity field. This is the same expression as
 * in the shared memory implementations.
 */
void leapfrog_master_uv(double *restrict p, double *restrict u,
                        double *restrict v, long nx, double dt, double dx,
                        double dy, long size)
{
    for (long i = 1; i < nx; ++i)
        for (long j = 0; j < size; ++j)
            U(i, j) += dt / dx * (P(i, j) - P(i - 1, j));

    for (long i = 0; i < nx; ++i)
        for (long j = 1; j < size; ++j)
            V(i, j) += dt / dy * (P(i, j) - P(i, j - 1));
}

/*
 * Compute one time step for the velocity field on the worker nodes. This
 * differs since we now have extra ghost cells to take into consideration. In
 * particular, grid indices for pressure needs to be shifter up by one along
 * the y-axis.
 */
void leapfrog_worker_uv(double *restrict p, double *restrict u,
                        double *restrict v, long nx, double dt, double dx,
                        double dy, long size)
{
    for (long i = 1; i < nx; ++i)
        for (long j = 0; j < size; ++j)
            U(i, j) += dt / dx * (P(i, j + 1) - P(i - 1, j + 1));

    for (long i = 0; i < nx; ++i)
        for (long j = 0; j < size; ++j)
            V(i, j) += dt / dy * (P(i, j + 1) - P(i, j));
}

/*
 * Before/after each time step we need to send the boundary data between each
 * computational node.
 *  Send v to the left (bottom),
 *  Receive v from the right (top). 
 */
void communicate_v(double *v, long left, long right, long nx, long size,
                   MPI_Status * status)
{
    if (left != NONE) {
        long begin_v = 0 + 0 * nx;
        long size_v = 1 * nx;
        MPI_Send(&v[begin_v], size_v, MPI_DOUBLE, left, VTAG,
                 MPI_COMM_WORLD);
    }
    if (right != NONE) {
        long begin_v = 0 + size * nx;
        long size_v = 1 * nx;
        MPI_Recv(&v[begin_v], size_v, MPI_DOUBLE, right, VTAG,
                 MPI_COMM_WORLD, status);
    }
}

/*
 * Before/after each time step we need to send the boundary data between each
 * computational node.
 *  Receive p from the left (bottom)
 *  Send p to the right (top)
 */
void communicate_p(double *p, long left, long right, long nx, long size,
                   MPI_Status * status)
{
    if (left != NONE) {
        long begin_p = 0 + 0 * nx;
        long size_p = 1 * nx;
        MPI_Recv(&p[begin_p], size_p, MPI_DOUBLE, left, PTAG,
                 MPI_COMM_WORLD, status);
    }
    if (right != NONE) {
        /* For the master node we don't have an extra ghost point p to the 
         * left (below) */
        long begin_p =
            (left == NONE) ? 0 + (size - 1) * nx : 0 + size * nx;
        long size_p = 1 * nx;
        MPI_Send(&p[begin_p], size_p, MPI_DOUBLE, right, PTAG,
                 MPI_COMM_WORLD);
    }
}

/*
 * The full update done in each time step. This includes sending/receiving data between the computational nodes
 * as well as performing the leapfrog update of the field.
 */
void leapfrog_mpi(struct field *f, long left, long right, long size)
{
    MPI_Status status;

    /* Used by index macro */
    double *p = f->p.value;
    double *u = f->u.value;
    double *v = f->v.value;
    long nx = f->p.size_x;

    /* Communicate: */
    /* send v to the left (bottom), receive v from the right (top) */
    communicate_v(v, left, right, nx, size, &status);

    /* Update the pressure (p) */
    if (left == NONE)
        leapfrog_master_p(p, u, v, nx, f->dt, f->u.dx, f->v.dy, size);
    else
        leapfrog_worker_p(p, u, v, nx, f->dt, f->u.dx, f->v.dy, size);

    /* Communicate: */
    /* send p to the right (top), receive p from the left (bottom) */
    communicate_p(p, left, right, nx, size, &status);

    /* Update the velocity (u,v) */
    if (left == NONE)
        leapfrog_master_uv(p, u, v, nx, f->dt, f->p.dx, f->p.dy, size);
    else
        leapfrog_worker_uv(p, u, v, nx, f->dt, f->p.dx, f->p.dy, size);
}

/*
 * Update the field data from the start time to the end time.
 */
void timestep(struct field *f, long left, long right, long size)
{
    for (long n = 0; n < f->Nt; ++n)
        leapfrog_mpi(f, left, right, size);
}
