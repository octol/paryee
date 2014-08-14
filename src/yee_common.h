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

#define STR_SIZE 256

#define P(x, y) (p[(x) + (y)*nx])
#define U(x, y) (u[(x) + (y)*(nx + 1)])
#define V(x, y) (v[(x) + (y)*nx])

/***************************************************************************
 * Data structures
 **************************************************************************/

/*
 * A py_field variable is scalar valued and defined on a 2D grid.
 */
struct py_field_variable {
    double *value;
    double *x;
    double *y;
    double dx;
    double dy;
    long size_x;                /* number of points */
    long size_y;                /* number of points */
};

/*
 * Contains the entire field
 */
struct py_field {
    struct py_field_variable p;
    struct py_field_variable u;
    struct py_field_variable v;
    double dt;
    double Nt;
};

/*
 * Specifies the start and end cell index of the internal domain (along the
 * x-axis).
 */
struct py_cell_partition {
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
void py_alloc_field(struct py_field_variable *, const long size_x,
                 const long size_y);
void py_free_field(struct py_field_variable);

/* 
 * Generate grid from start to end
 */
void py_set_grid(struct py_field_variable *, const double x[2],
              const double y[2]);

/*
 * Vectorization of function.
 * Apply the function func on the array dst with argument arg.
 */
void py_vec_func(double *dst, double (*func) (double), const double *arg,
              const long size);

void py_vec_func2d(double *dst, double (*func) (double, double),
                const double *arg_x, const long size_x,
                const double *arg_y, const long size_y);

/*
 * Vectorization of function applied to py_field_variable.
 */
void py_apply_func(struct py_field_variable *f, double (*func) (double, double));

/* 
 * Allocates py_field and sets initial value to zero
 */
struct py_field py_init_acoustic_field(long cells_x,
                                 long cells_y, double x[2], double y[2]);

/* 
 * Allocates py_field and sets initial value to zero. This variant is for
 * distributed memory local domains. The difference is we get one more p
 * node to the left (below). This the generated grid will actually start
 * dy/2 before y[0].
 * NOTE: Function not valid for the left most (bottom) partition.
 */
struct py_field py_init_local_acoustic_field(long cells_x, long cells_y,
                                       double x[2], double y[2]);

/*
 * Deallocate field
 */
void py_free_acoustic_field(struct py_field);

/*
 * Assign and retrive value from fv(i,j)
 */
double py_assign_to(struct py_field_variable fv, long i, long j, double value);
double get_from(struct py_field_variable fv, long i, long j);

/*
 * Sets the outer boundary (u,v) to zero.
 */
void set_boundary(struct py_field *f);

/*
 * Leapfrog update of the py_field 
 */
void leapfrog(struct py_field *f);

/*
 * Time step n Leapfrog updates of the py_field 
 */
void timestep_leapfrog(struct py_field *f, double n);

/*
 * Divide the grid for the different threads. The partition is done as
 * strips. We divide the grid according to cells.
 * Note: only the inner nodes are returned.
 */
struct py_cell_partition *partition_grid(long total_threads, long cells);

/*
 * Given partition information, get the corresponding domain coordinates.
 * NOTE: this is for when we divide on the y-axis for the MPI version.
 */
void get_partition_coords(struct py_cell_partition part, struct py_field *f,
                          double *y);

/*
 * Convert the cell indices that specifies the partition, to node indices
 * that specifies which nodes to update (leapfrog update / time step).
 * Note that the end index is such that the indices <p1, <u1, <v1 are
 * update.
 */
void cellindex_to_nodeindex(long tid, struct py_cell_partition part,
                            long *p0, long *p1,
                            long *u0, long *u1, long *v0, long *v1);

/*
 * Parse commandline argument and set the number of grid points, * the
 * number of threads as well as the output filenames.
 */
void parse_cmdline(long *nx, long *threads, char *outfile, int *write,
                   int argc, char *argv[]);

/*
 * Write py_field to dist
 */
int write_to_disk(struct py_field_variable f, char *fstr);

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
