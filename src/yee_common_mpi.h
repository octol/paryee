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
 * Data structures and functions common to the various implementations.
 */
#ifndef YEE_COMMON_MPI_H
#define YEE_COMMON_MPI_H

#include "mpi.h"
#include "yee_common.h"

#define BEGIN 1                 /* message tag */
#define COLLECT 2               /* message tag */
#define PTAG 3                  /* communication tag */
#define UTAG 4                  /* communication tag */
#define VTAG 5                  /* communication tag */

#define MASTER 0                /* id of master process */
#define NONE -1                 /* no neighbour */

/*
 * Sends the grid data from the master node to the worker nodes.
 */
void send_grid_data(int taskid, long left, long right, long size,
                    double y[2]);

/*
 * Each worker node receives its local grid data from the master node.
 */
void receive_grid_data(long *left, long *right, long *size, double *y,
                       MPI_Status * status);

/*
 * The actual py_field (pressure and velocity) is sent from master to worker node. 
 */
void send_field_data(long taskid, struct py_field *f,
                     struct py_cell_partition part);

/*
 * Each worker node receives its py_field (pressure and velocity) from the master
 * node.
 */
void receive_field_data(struct py_field *f, long size,
                        MPI_Status * status);

/* 
 * After simulation/time stepping is done, each worker node sends its field
 * data back to the master node.
 */
void return_data(struct py_field *f, long size);

/*
 * The master collects all data into one unit.
 */
void collect_data(long taskid, struct py_cell_partition part,
                  struct py_field *f, MPI_Status * status);

/* 
 * Compute one time step for the pressure. This is the same expression as in
 * the shared memory implementations.
 */
void leapfrog_master_p(double *restrict p, double *restrict u,
                       double *restrict v, long nx, double dt, double dx,
                       double dy, long size);
/* Non-blocking version */
void leapfrog_master_p2(double *restrict p, double *restrict u,
                        double *restrict v, long nx, double dt, double dx,
                        double dy, long size, long right,
                        MPI_Request * recived);

/*
 * Compute one time step for the pressure on the worker nodes. This differs
 * since we now have extra ghost cells to take into consideration. In
 * particular, grid indices for pressure needs to be shifter up by one along
 * the y-axis.
 */
void leapfrog_worker_p(double *restrict p, double *restrict u,
                       double *restrict v, long nx, double dt, double dx,
                       double dy, long size);
/* Non-blocking version */
void leapfrog_worker_p2(double *restrict p, double *restrict u,
                        double *restrict v, long nx, double dt, double dx,
                        double dy, long size, long right,
                        MPI_Request * received);

/* 
 * Compute one time step for the velocity field. This is the same expression as
 * in the shared memory implementations.
 */
void leapfrog_master_uv(double *restrict p, double *restrict u,
                        double *restrict v, long nx, double dt, double dx,
                        double dy, long size);

/*
 * Compute one time step for the velocity py_field on the worker nodes. This
 * differs since we now have extra ghost cells to take into consideration. In
 * particular, grid indices for pressure needs to be shifter up by one along
 * the y-axis.
 */
void leapfrog_worker_uv(double *restrict p, double *restrict u,
                        double *restrict v, long nx, double dt, double dx,
                        double dy, long size);
/* Non-blocking version */
void leapfrog_worker_uv2(double *restrict p, double *restrict u,
                         double *restrict v, long nx, double dt, double dx,
                         double dy, long size, long left,
                         MPI_Request * recieved);

/*
 * Before/after each time step we need to send the boundary data between each
 * computational node.
 *  Send v to the left (bottom),
 *  Receive v from the right (top). 
 */
void communicate_v(double *v, long left, long right, long nx, long size,
                   MPI_Status * status);
/* Non-blocking version */
void communicate_v2(double *v, long left, long right, long nx, long size,
                    MPI_Request * sent, MPI_Request * received);

/*
 * Before/after each time step we need to send the boundary data between each
 * computational node.
 *  Receive p from the left (bottom)
 *  Send p to the right (top)
 */
void communicate_p(double *p, long left, long right, long nx, long size,
                   MPI_Status * status);
/* Non-blocking version */
void communicate_p2(double *p, long left, long right, long nx, long size,
                    MPI_Request * sent, MPI_Request * received);

/*
 * The full update done in each time step. This includes sending/receiving
 * data between the computational nodes as well as performing the leapfrog
 * update of the field.
 */
void leapfrog_mpi(struct py_field *f, long left, long right, long size);
/* Non-blocking version */
void leapfrog_mpi2(struct py_field *f, long left, long right, long size);

/*
 * Update the py_field data from the start time to the end time.
 */
void timestep_mpi(struct py_field *f, long left, long right, long size);
/* Non-blocking version */
void timestep_mpi2(struct py_field *f, long left, long right, long size);

#endif
