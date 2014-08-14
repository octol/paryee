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

void send_field_data(long taskid, struct py_field *f,
                     struct cell_partition part)
{
    /* Send out py_field data.
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

void receive_field_data(struct py_field *f, long size, MPI_Status * status)
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

void return_data(struct py_field *f, long size)
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

void collect_data(long taskid, struct cell_partition part, struct py_field *f,
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

void leapfrog_master_p2(double *restrict p, double *restrict u,
                        double *restrict v, long nx, double dt, double dx,
                        double dy, long size, long right,
                        MPI_Request * received)
{
    for (long i = 0; i < nx; ++i)
        for (long j = 0; j < size - 1; ++j)
            P(i, j) +=
                dt / dx * (U(i + 1, j) - U(i, j)) +
                dt / dy * (V(i, j + 1) - V(i, j));

    /* Bordering py_field points */
    if (right != NONE)
        MPI_Wait(received, MPI_STATUS_IGNORE);

    long j = size - 1;
    for (long i = 0; i < nx; ++i)
        P(i, j) +=
            dt / dx * (U(i + 1, j) - U(i, j)) +
            dt / dy * (V(i, j + 1) - V(i, j));
}

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

void leapfrog_worker_p2(double *restrict p, double *restrict u,
                        double *restrict v, long nx, double dt, double dx,
                        double dy, long size, long right,
                        MPI_Request * received)
{
    for (long i = 0; i < nx; ++i)
        for (long j = 0; j < size - 1; ++j)
            P(i, j + 1) +=
                dt / dx * (U(i + 1, j) - U(i, j)) +
                dt / dy * (V(i, j + 1) - V(i, j));

    /* Bordering py_field points */
    if (right != NONE)
        MPI_Wait(received, MPI_STATUS_IGNORE);

    long j = size - 1;
    for (long i = 0; i < nx; ++i)
        P(i, j + 1) +=
            dt / dx * (U(i + 1, j) - U(i, j)) +
            dt / dy * (V(i, j + 1) - V(i, j));
}

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

void leapfrog_worker_uv2(double *restrict p, double *restrict u,
                         double *restrict v, long nx, double dt, double dx,
                         double dy, long size, long left,
                         MPI_Request * received)
{
    for (long i = 1; i < nx; ++i)
        for (long j = 0; j < size; ++j)
            U(i, j) += dt / dx * (P(i, j + 1) - P(i - 1, j + 1));
    for (long i = 0; i < nx; ++i)
        for (long j = 1; j < size; ++j)
            V(i, j) += dt / dy * (P(i, j + 1) - P(i, j));

    /* Bordering py_field points */
    if (left != NONE)
        MPI_Wait(received, MPI_STATUS_IGNORE);

    long j = 0;
    for (long i = 0; i < nx; ++i)
        V(i, j) += dt / dy * (P(i, j + 1) - P(i, j));
}

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

void communicate_v2(double *v, long left, long right, long nx, long size,
                    MPI_Request * sent, MPI_Request * received)
{
    if (left != NONE) {
        long begin_v = 0 + 0 * nx;
        long size_v = 1 * nx;
        MPI_Isend(&v[begin_v], size_v, MPI_DOUBLE, left, VTAG,
                  MPI_COMM_WORLD, sent);
    }
    if (right != NONE) {
        long begin_v = 0 + size * nx;
        long size_v = 1 * nx;
        MPI_Irecv(&v[begin_v], size_v, MPI_DOUBLE, right, VTAG,
                  MPI_COMM_WORLD, received);
    }
}

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

void communicate_p2(double *p, long left, long right, long nx, long size,
                    MPI_Request * sent, MPI_Request * received)
{
    if (left != NONE) {
        long begin_p = 0 + 0 * nx;
        long size_p = 1 * nx;
        MPI_Irecv(&p[begin_p], size_p, MPI_DOUBLE, left, PTAG,
                  MPI_COMM_WORLD, received);
    }
    if (right != NONE) {
        /* For the master node we don't have an extra ghost point p to the 
         * left (below) */
        long begin_p =
            (left == NONE) ? 0 + (size - 1) * nx : 0 + size * nx;
        long size_p = 1 * nx;
        MPI_Isend(&p[begin_p], size_p, MPI_DOUBLE, right, PTAG,
                  MPI_COMM_WORLD, sent);
    }
}

void leapfrog_mpi(struct py_field *f, long left, long right, long size)
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

void leapfrog_mpi2(struct py_field *f, long left, long right, long size)
{
    MPI_Request sent;
    MPI_Request received;

    /* Used by index macro */
    double *p = f->p.value;
    double *u = f->u.value;
    double *v = f->v.value;
    long nx = f->p.size_x;

    /* Communicate: */
    /* send v to the left (bottom), receive v from the right (top) */
    communicate_v2(v, left, right, nx, size, &sent, &received);

    /* Update the pressure (p) */
    if (left == NONE)
        leapfrog_master_p2(p, u, v, nx, f->dt, f->u.dx, f->v.dy, size,
                           right, &received);
    else
        leapfrog_worker_p2(p, u, v, nx, f->dt, f->u.dx, f->v.dy, size,
                           right, &received);
    if (left != NONE)
        MPI_Wait(&sent, MPI_STATUS_IGNORE);

    /* Communicate: */
    /* send p to the right (top), receive p from the left (bottom) */
    communicate_p2(p, left, right, nx, size, &sent, &received);

    /* Update the velocity (u,v) */
    if (left == NONE)
        leapfrog_master_uv(p, u, v, nx, f->dt, f->p.dx, f->p.dy, size);
    else
        leapfrog_worker_uv2(p, u, v, nx, f->dt, f->p.dx, f->p.dy, size,
                            left, &received);
    if (right != NONE)
        MPI_Wait(&sent, MPI_STATUS_IGNORE);
}

void timestep_mpi(struct py_field *f, long left, long right, long size)
{
    for (long n = 0; n < f->Nt; ++n)
        leapfrog_mpi(f, left, right, size);
}

void timestep_mpi2(struct py_field *f, long left, long right, long size)
{
    for (long n = 0; n < f->Nt; ++n)
        leapfrog_mpi2(f, left, right, size);
}
