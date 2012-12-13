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
#define MINWORKER 1
#define MAXWORKER 8
#define NONE -1                 /* no neighbour */
#define STR_SIZE 256

#define p(x, y) (p[(x) + (y)*nx])
#define u(x, y) (u[(x) + (y)*(nx + 1)])
#define v(x, y) (v[(x) + (y)*nx])

/***************************************************************************
 * Data structures
 **************************************************************************/
struct field_variable {
    double *value;
    double *x;
    double *y;
    double dx;
    double dy;
    unsigned long size_x;                  /* number of points */
    unsigned long size_y;                  /* number of points */
};

struct field {
    struct field_variable p;
    struct field_variable u;
    struct field_variable v;
    double dt;
    double Nt;
};

struct partition {
    unsigned long p[2];                  /* start end indices of internal domain */
    unsigned long u[2];
    unsigned long v[2];
};

/***************************************************************************
 * Functions
 **************************************************************************/

/* 
 * Allocate and deallocate the memory used by the grid (field) 
 */
void alloc_field(struct field_variable *, const unsigned long size_x, const unsigned long size_y);
void free_field(struct field_variable);

/* 
 * Generate grid from start to end
 */
void set_grid(struct field_variable *, const double x[2], const double y[2]);

/*
 * Vectorization of function.
 * Apply the function func on the array dst with argument arg.
 */
void vec_func(double *dst, double (*func) (double), const double *arg, 
              const unsigned long size);

void vec_func2d(double *dst, double (*func) (double, double), 
                const double *arg_x, const unsigned long size_x, 
                const double *arg_y, const unsigned long size_y);

/*
 * Vectorization of function applied to field_variable.
 */
void apply_func(struct field_variable *f, double (*func) (double, double));

/* 
 * Allocates field and sets initial value to zero
 */
struct field init_acoustic_field(unsigned long cells_x, unsigned long cells_y, double x[2], double y[2]);

/*
 * Deallocate field
 */
void free_acoustic_field(struct field);

/*
 * Assign and retrive value from fv(i,j)
 */
double assign_to(struct field_variable fv, unsigned long i, unsigned long j, double value);
double get_from(struct field_variable fv, unsigned long i, unsigned long j);

/*
 * Leapfrog update of the field 
 */
void leapfrog(struct field *f);

/*
 * Time step n Leapfrog updates of the field 
 */
void timestep_leapfrog(struct field *f, unsigned long n);

/*
 * Divide the grid for the different threads.
 * We divide the grid according to cells.
 * Note: only the inner nodes are returned.
 */
struct partition partition_grid(int current_thread, int cells_per_thread);

/* 
 * Expand the partition struct into indices 
 */
void expand_indices(struct partition partition,
                    long *begin_p, long *end_p,
                    long *size_p, long *begin_u,
                    long *end_u, long *size_u);

/*
 * Checks that the grid is according to the specs.
 */
void verify_grid_integrity(struct partition partition, int tid,
                           long nx, int numworkers, int left);

/*
 * Compute the local start/end indices and local array sizes from the
 * partition sizes as well as if the we are on the left boundary.
 */
void set_local_index(long size_p, long size_u,
                     long left, long *local_begin_p,
                     long *local_end_p,
                     long *local_size_p,
                     long *local_begin_u,
                     long *local_end_u, long *local_size_u);

/*
 * Parse commandline argument and set the number of grid points, * the
 * number of threads as well as the output filenames.
 */
void parse_cmdline(unsigned long *nx, unsigned long *threads,
                   char *outfile, int argc, char *argv[]);

/*
 * Write field to dist
 */
int write_to_disk(struct field_variable f, char *fstr);
int write_field_to_disk(struct field f, char *p_str, char *u_str);

/* 
 * Some grid functions 
 */
double gauss(double);
double zero(double);
double zero2d(double x, double y);
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
