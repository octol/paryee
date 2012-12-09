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
#define NONE 0                  /* no neighbour */
#define STR_SIZE 256

/***************************************************************************
 * Data structures
 **************************************************************************/
struct field_variable {
    double *value;
    double *x;
    double dx;
    long size;                  /* number of points */
};

struct field {
    struct field_variable p;
    struct field_variable u;
    double dt;
    double Nt;
};

struct partition {
    long p[2];                  /* start end indices of internal domain */
    long u[2];
};

/***************************************************************************
 * Functions
 **************************************************************************/

/* 
 * Allocate and deallocate the memory used by the grid (field) 
 */
void alloc_field(struct field_variable *, long size);
void free_field(struct field_variable);

/* 
 * Generate grid from start to end with N number of points  
 */
void set_grid(struct field_variable *, double start, double end);

/*
 * Vectorization of function.
 * Apply the function func on the array dst with argument arg.
 */
void vec_func(double *dst, double *arg, double (*func) (double),
              long i1, long i2);

/*
 * Vectorization of function applied to field_variable.
 */
void apply_func(struct field_variable *f, double (*func) (double));

/* 
 * Leapfrog time update.
 * Update the size number of field points starting at the index idst, using
 * the field points starting at isrc. Note: _s in the name indicates the
 * input (size).
 */
void update_field_s(struct field_variable *restrict dst, int idst,
                    int size, struct field_variable *restrict src,
                    int isrc, double dt);

/* 
 * Leapfrog time update.
 * Update the field points starting at the index dst1 and ending at dst2,
 * using the field points starting at src1. Note: _i in the name indicates
 * the input (index).
 */
void update_field_i(struct field_variable *restrict dst, int dst1,
                    int dst2, struct field_variable *restrict src,
                    int src1, double dt);

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
void parse_cmdline(long *nx, long *threads,
                   char *outfile_p, char *outfile_u,
                   int argc, char *argv[]);

/*
 * Write field to dist
 */
int write_to_disk(struct field_variable f, char *str);

/* 
 * Some grid functions 
 */
double gauss(double);
double zero(double);

/*
 * Round up integer division
 */
int round_up_divide(int x, int y);

/*
 * Get time in a way that is compatible with threading
 */
double gettime(void);

#endif
