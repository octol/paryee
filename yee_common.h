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
#ifndef YEE_COMMON_H
#define YEE_COMMON_H

#define MASTER 0                /* id of master process */
#define NONE -1                 /* no neighbour */
#define STR_SIZE 256

#define P(x, y) (p[(x) + (y)*nx])
#define U(x, y) (u[(x) + (y)*(nx + 1)])
#define V(x, y) (v[(x) + (y)*nx])

/***************************************************************************
 * Data structures
 **************************************************************************/

/*
 * A field variable is scalar valued and defined on a 2D grid.
 */
struct field_variable {
    double *value;
    double *x;
    double *y;
    double dx;
    double dy;
    long size_x;       /* number of points */
    long size_y;       /* number of points */
};

/*
 * Contains the entire field
 */
struct field {
    struct field_variable p;
    struct field_variable u;
    struct field_variable v;
    double dt;
    double Nt;
};

/*
 * Specifies the start and end cell index of the internal domain (along the
 * x-axis).
 */
struct cell_partition {
    long begin;
    long end;
    long size;
};

/***************************************************************************
 * Functions
 **************************************************************************/

/* 
 * Allocate and deallocate the memory used by the grid (field) 
 */
void alloc_field(struct field_variable *, const long size_x,
                 const long size_y);
void free_field(struct field_variable);

/* 
 * Generate grid from start to end
 */
void set_grid(struct field_variable *, const double x[2],
              const double y[2]);

/*
 * Vectorization of function.
 * Apply the function func on the array dst with argument arg.
 */
void vec_func(double *dst, double (*func) (double), const double *arg,
              const long size);

void vec_func2d(double *dst, double (*func) (double, double),
                const double *arg_x, const long size_x,
                const double *arg_y, const long size_y);

/*
 * Vectorization of function applied to field_variable.
 */
void apply_func(struct field_variable *f, double (*func) (double, double));

/* 
 * Allocates field and sets initial value to zero
 */
struct field init_acoustic_field(long cells_x,
                                 long cells_y, double x[2],
                                 double y[2]);

/* 
 * Allocates field and sets initial value to zero. This variant is for
 * distributed memory local domains. The difference is we get one more p
 * node to the left (below). This the generated grid will actually start
 * dy/2 before y[0].
 * NOTE: Function not valid for the left most (bottom) partition.
 */
struct field init_local_acoustic_field(long cells_x, long cells_y, double x[2], double y[2]);

/*
 * Deallocate field
 */
void free_acoustic_field(struct field);

/*
 * Assign and retrive value from fv(i,j)
 */
double assign_to(struct field_variable fv, long i,
                 long j, double value);
double get_from(struct field_variable fv, long i,
                long j);

/*
 * Sets the outer boundary (u,v) to zero.
 */
void set_boundary(struct field *f);

/*
 * Leapfrog update of the field 
 */
void leapfrog(struct field *f);

/*
 * Time step n Leapfrog updates of the field 
 */
void timestep_leapfrog(struct field *f, double n);

/*
 * Divide the grid for the different threads. The partition is done as
 * strips. We divide the grid according to cells.
 * Note: only the inner nodes are returned.
 */
struct cell_partition *partition_grid(long total_threads,
                                      long cells);

/*
 * Given partition information, get the corresponding domain coordinates.
 * NOTE: this is for when we divide on the y-axis for the MPI version.
 */
void get_partition_coords(struct cell_partition part, struct field *f, double *y);

/*
 * Convert the cell indices that specifies the partition, to node indices
 * that specifies which nodes to update (leapfrog update / time step).
 * Note that the end index is such that the indices <p1, <u1, <v1 are
 * update.
 */
void cellindex_to_nodeindex(long tid, struct cell_partition part,
                            long *p0, long *p1,
                            long *u0, long *u1,
                            long *v0, long *v1);

/* 
 * Expand the partition struct into indices 
 */
/*void expand_indices(struct cell_partition partition,
                    long *begin_p, long *end_p,
                    long *size_p, long *begin_u,
                    long *end_u, long *size_u);*/

/*
 * Checks that the grid is according to the specs.
 */
/*void verify_grid_integrity(struct cell_partition partition, long tid,
                           long nx, long numworkers,
                           long left);*/

/*
 * Compute the local start/end indices and local array sizes from the
 * partition sizes as well as if the we are on the left boundary.
 */
/*void set_local_index(long size_p, long size_u,
                     long left, long *local_begin_p,
                     long *local_end_p,
                     long *local_size_p,
                     long *local_begin_u,
                     long *local_end_u,
                     long *local_size_u);*/

/*
 * Parse commandline argument and set the number of grid points, * the
 * number of threads as well as the output filenames.
 */
void parse_cmdline(long *nx, long *threads,
                   char *outfile, int argc, char *argv[]);

/*
 * Write field to dist
 */
int write_to_disk(struct field_variable f, char *fstr);

/* 
 * Some grid functions 
 */
double gauss(double);
double gauss2d(double x, double y);
double zero(double x);
double zero2d(double x, double y);
double one2d(double x, double y);
double identity(double);
double identity2d(double x, double y);

/*
 * Round up integer division
 */
int round_up_divide(int x, int y);

/*
 * Get time in a way that is compatible with threading
 */
double gettime(void);

#endif
