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
#include "yee_common.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <getopt.h>
#include <assert.h>

void alloc_field(struct field_variable *f, const unsigned long size_x,
                 const unsigned long size_y)
{
    unsigned long block = sizeof(double) * size_x * size_y;
    f->size_x = size_x;
    f->size_y = size_y;
    f->value = (double *) malloc(block);
    f->x = (double *) malloc(sizeof(double) * size_x);
    f->y = (double *) malloc(sizeof(double) * size_y);

    if (!f->value || !f->x || !f->y) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

void free_field(struct field_variable f)
{
    free(f.value);
    free(f.x);
    free(f.y);
}

void set_grid(struct field_variable *f, const double x[2],
              const double y[2])
{
    unsigned long i;

    /* NOTE: f->size is the number of grid points, not intervals */
    f->dx = (x[1] - x[0]) / ((double) f->size_x - 1.0);
    f->dy = (y[1] - y[0]) / ((double) f->size_y - 1.0);
    for (i = 0; i < f->size_x; ++i)
        f->x[i] = x[0] + (double) i *f->dx;
    for (i = 0; i < f->size_y; ++i)
        f->y[i] = y[0] + (double) i *f->dy;
}

void vec_func(double *dst, double (*func) (double), const double *arg,
              const unsigned long size)
{
    for (unsigned long i = 0; i < size; ++i)
        dst[i] = func(arg[i]);
}

void vec_func2d(double *dst, double (*func) (double, double),
                const double *arg_x, const unsigned long size_x,
                const double *arg_y, const unsigned long size_y)
{
    for (unsigned long ii = 0; ii < size_x; ++ii)
        for (unsigned long jj = 0; jj < size_y; ++jj)
            dst[ii + jj * size_x] = func(arg_x[ii], arg_y[jj]);
}

void apply_func(struct field_variable *f, double (*func) (double, double))
{
    vec_func2d(f->value, func, f->x, f->size_x, f->y, f->size_y);
}

struct field init_acoustic_field(unsigned long cells_x,
                                 unsigned long cells_y, double x[2],
                                 double y[2])
{
    struct field f;

    alloc_field(&f.p, cells_x, cells_y);
    alloc_field(&f.u, cells_x + 1, cells_y);
    alloc_field(&f.v, cells_x, cells_y + 1);

    double dx = (x[1] - x[0]) / (double) cells_x;
    double dy = (y[1] - y[0]) / (double) cells_y;
    double p_x[2] = { x[0] + dx / 2.0, x[1] - dx / 2.0 };
    double p_y[2] = { y[0] + dy / 2.0, y[1] - dy / 2.0 };
    double u_x[2] = { x[0], x[1] };
    double u_y[2] = { y[0] + dy / 2.0, y[1] - dy / 2.0 };
    double v_x[2] = { x[0] + dx / 2.0, x[1] - dx / 2.0 };
    double v_y[2] = { y[0], y[1] };

    assert(fabs(p_x[0] - (x[0] + dx / 2.0)) < 1e-14);
    assert(fabs(p_y[0] - (y[0] + dy / 2.0)) < 1e-14);
    assert(fabs(u_x[0] - x[0]) < 1e-14);
    assert(fabs(u_y[0] - (y[0] + dy / 2.0)) < 1e-14);
    assert(fabs(v_x[0] - (x[0] + dx / 2.0)) < 1e-14);
    assert(fabs(v_y[0] - y[0]) < 1e-14);

    set_grid(&f.p, p_x, p_y);
    set_grid(&f.u, u_x, u_y);
    set_grid(&f.v, v_x, v_y);
    apply_func(&f.p, zero2d);
    apply_func(&f.u, zero2d);
    apply_func(&f.v, zero2d);

    return f;
}

void free_acoustic_field(struct field f)
{
    free_field(f.p);
    free_field(f.u);
    free_field(f.v);
}

double assign_to(struct field_variable fv, unsigned long i,
                 unsigned long j, double value)
{
    fv.value[i + j * fv.size_x] = value;
    return value;
}

double get_from(struct field_variable fv, unsigned long i, unsigned long j)
{
    return fv.value[i + j * fv.size_x];
}

void leapfrog(struct field *f)
{
    unsigned long i, j;
    unsigned long nx = f->p.size_x;
    unsigned long ny = f->p.size_y;
    double dt = f->dt;
    double *p = f->p.value;
    double *u = f->u.value;
    double *v = f->v.value;

    for (i = 0; i < nx; ++i) {
        for (j = 0; j < ny; ++j) {
            p(i, j) +=
#ifdef EXTRA_WORK
                /* NOTE: the pow(exp(-pow(sin(dt),2)),3.2) factor is to make
                 * more operations per memory access. This is to test the
                 * parallel performance.  */
                pow(exp(-pow(sin(dt), 2)), 3.2) *
                pow(exp(-pow(sin(dt), 4)), 1.2) *
                pow(exp(-pow(sin(dt), 2)), 4.2) *
#endif
                dt / f->u.dx * (u(i + 1, j) - u(i, j)) +
                dt / f->v.dy * (v(i, j + 1) - v(i, j));
        }
    }

    for (i = 1; i < nx - 1; ++i)
        for (j = 0; j < ny; ++j) 
            u(i, j) += dt / f->p.dx * (p(i, j) - p(i - 1, j));

    for (i = 0; i < nx; ++i) 
        for (j = 1; j < ny - 1; ++j) 
            v(i, j) += dt / f->p.dy * (p(i, j) - p(i, j - 1));
}

void timestep_leapfrog(struct field *f, unsigned long Nt)
{
    for (unsigned long n = 0; n < Nt; ++n) 
        leapfrog(f);
}

/*struct partition partition_grid(int current_thread, int cells_per_thread)  */
/*{                                                                          */
   /*struct partition partition;                                            */

   /*partition.p[0] = current_thread * cells_per_thread;                    */
   /*partition.p[1] = (current_thread + 1) * cells_per_thread - 1;          */
   /*partition.u[0] = current_thread * cells_per_thread;                    */
   /*partition.u[1] = (current_thread + 1) * cells_per_thread - 1;          */

   /*[> First thread skips the first u, as it sits on the boundary <]       */
   /*if (current_thread == 0)                                               */
       /*++partition.u[0];                                                  */

   /*return partition;                                                      */
/*}                                                                          */

/*void expand_indices(struct partition partition,                             */
/*                    long *begin_p, long *end_p,                             */
/*                    long *size_p, long *begin_u, long *end_u, long *size_u) */
/*{                                                                           */
/*    *begin_p = partition.p[0];  [> begin index <]                           */
/*    *begin_u = partition.u[0];                                              */
/*    *end_p = partition.p[1];    [> end index <]                             */
/*    *end_u = partition.u[1];                                                */

/*    [> sizes <]                                                             */
/*    *size_p = *end_p - *begin_p + 1;                                        */
/*    *size_u = *end_u - *begin_u + 1;                                        */
/*}                                                                           */

/*void verify_grid_integrity(struct partition partition, int tid,             */
/*                           long nx, int numworkers, int left)               */
/*{                                                                           */
/*    long bp, ep, sp, bu, eu, su;                                            */
/*    expand_indices(partition, &bp, &ep, &sp, &bu, &eu, &su);                */

/*    [> Sanity checks <]                                                     */
/*    assert(bp == (tid) * nx / numworkers);                                  */
/*    assert(ep == (tid + 1) * nx / numworkers - 1);                          */
/*    assert(sp == nx / numworkers);                                          */
/*    if (left == NONE) {                                                     */
/*        assert(bu == 1);                                                    */
/*        assert(su == nx / numworkers - 1);                                  */
/*    } else {                                                                */
/*        assert(bu == (tid) * nx / numworkers);                              */
/*        assert(su == nx / numworkers);                                      */
/*    }                                                                       */
/*    assert(eu == (tid + 1) * nx / numworkers - 1);                          */

/*    [> Specific sanity checks <]                                            */
/*    if (nx == 8 && numworkers == 4) {                                       */
/*        switch (tid) {                                                      */
/*        case 0:                                                             */
/*            assert(bp == 0 && ep == 1 && sp == 2);                          */
/*            assert(bu == 1 && eu == 1 && su == 1);                          */
/*            break;                                                          */
/*        case 1:                                                             */
/*            assert(bp == 2 && ep == 3 && sp == 2);                          */
/*            assert(bu == 2 && eu == 3 && su == 2);                          */
/*            break;                                                          */
/*        case 2:                                                             */
/*            assert(bp == 4 && ep == 5 && sp == 2);                          */
/*            assert(bu == 4 && eu == 5 && su == 2);                          */
/*            break;                                                          */
/*        case 3:                                                             */
/*            assert(bp == 6 && ep == 7 && sp == 2);                          */
/*            assert(bu == 6 && eu == 7 && su == 2);                          */
/*            break;                                                          */
/*        }                                                                   */
/*    }                                                                       */
/*}                                                                           */

/*void set_local_index(long size_p, long size_u,                              */
/*                     long left, long *local_begin_p,                        */
/*                     long *local_end_p,                                     */
/*                     long *local_size_p,                                    */
/*                     long *local_begin_u,                                   */
/*                     long *local_end_u, long *local_size_u)                 */
/*{                                                                           */
/*    [> since p is padded by one point to the left <]                        */
/*    *local_begin_p = 1;                                                     */
/*    [> since u is not padded on the left <]                                 */
/*    *local_begin_u = 0;                                                     */
/*    if (left == NONE) {                                                     */
/*        *local_begin_p = 0;                                                 */
/*        *local_begin_u = 1;                                                 */
/*    }                                                                       */
/*    *local_end_p = *local_begin_p + size_p - 1;                             */
/*    *local_end_u = *local_begin_u + size_u - 1;                             */

/*    *local_size_p = *local_end_p + 1;                                       */
/*    *local_size_u = *local_end_u + 2;                                       */
/*}                                                                           */

void parse_cmdline(unsigned long *nx, unsigned long *threads,
                   char *outfile, int argc, char *argv[])
{
    int opt;
    while ((opt = getopt(argc, argv, "t:n:o:")) != -1) {
        switch (opt) {
        case 't':
            *threads = atoi(optarg);
            break;
        case 'n':
            *nx = atoi(optarg);
            break;
        case 'o':
            strncpy(outfile, optarg, STR_SIZE);
            outfile[STR_SIZE - 1] = '\0';       /* force null termination */
            break;
        default:               /* '?' */
            fprintf(stderr, "Usage: %s ", argv[0]);
            if (threads != NULL)
                fprintf(stderr, "[-t threads] ");
            fprintf(stderr, "[-n intervals] ");
            fprintf(stderr, "[-o outfile]\n");
            exit(EXIT_FAILURE);
        }
    }
}

int write_to_disk(struct field_variable f, char *fstr)
{
    FILE *fp;

    printf("Writing to: %s\n", fstr);
    fp = fopen(fstr, "w");
    if (fp == NULL) {
        perror("Error: can not write to disk");
        return EXIT_FAILURE;
    }
    for (unsigned long i = 0; i < f.size_x; ++i)
        for (unsigned long j = 0; j < f.size_y; ++j)
            fprintf(fp, "%e\t%e\t%e\n", f.x[i], f.y[j], get_from(f, i, j));
    fclose(fp);
    return EXIT_SUCCESS;
}

double gauss(double x)
{
    return exp(-pow(x - 0.5, 2) / pow(0.05, 2));
}

double gauss2d(double x, double y)
{
    double r_squared = pow(x-0.5,2) + pow(y-0.5,2);
    return exp(-r_squared / pow(0.05, 2));
}

double zero(double x)
{
    return 0.0;
}

double zero2d(double x, double y)
{
    return 0.0;
}

double identity(double x)
{
    return x;
}

double identity2d(double x, double y)
{
    return x;
}

int round_up_divide(int x, int y)
{
    return (x - 1) / y + 1;
}

double gettime(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1.0e-6 * t.tv_usec + t.tv_sec);
}
