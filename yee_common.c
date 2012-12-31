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

void alloc_field(struct field_variable *f, const long size_x,
                 const long size_y)
{
    long block = sizeof(double) * size_x * size_y;
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
    long i;

    /* NOTE: f->size is the number of grid points, not intervals */
    f->dx = (x[1] - x[0]) / ((double) f->size_x - 1.0);
    f->dy = (y[1] - y[0]) / ((double) f->size_y - 1.0);
    for (i = 0; i < f->size_x; ++i)
        f->x[i] = x[0] + (double) i *f->dx;
    for (i = 0; i < f->size_y; ++i)
        f->y[i] = y[0] + (double) i *f->dy;
}

void vec_func(double *dst, double (*func) (double), const double *arg,
              const long size)
{
    for (long i = 0; i < size; ++i)
        dst[i] = func(arg[i]);
}

void vec_func2d(double *dst, double (*func) (double, double),
                const double *arg_x, const long size_x,
                const double *arg_y, const long size_y)
{
    for (long ii = 0; ii < size_x; ++ii)
        for (long jj = 0; jj < size_y; ++jj)
            dst[ii + jj * size_x] = func(arg_x[ii], arg_y[jj]);
}

void apply_func(struct field_variable *f, double (*func) (double, double))
{
    vec_func2d(f->value, func, f->x, f->size_x, f->y, f->size_y);
}

struct field init_acoustic_field(long cells_x,
                                 long cells_y, double x[2],
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

    assert(fabs(p_x[1] - (x[1] - dx / 2.0)) < 1e-14);
    assert(fabs(p_y[1] - (y[1] - dy / 2.0)) < 1e-14);
    assert(fabs(u_x[1] - x[1]) < 1e-14);
    assert(fabs(u_y[1] - (y[1] - dy / 2.0)) < 1e-14);
    assert(fabs(v_x[1] - (x[1] - dx / 2.0)) < 1e-14);
    assert(fabs(v_y[1] - y[1]) < 1e-14);

    set_grid(&f.p, p_x, p_y);
    set_grid(&f.u, u_x, u_y);
    set_grid(&f.v, v_x, v_y);

    assert(fabs(f.p.dx - f.u.dx) < 1e-14);
    assert(fabs(f.p.dx - f.v.dx) < 1e-14);
    assert(fabs(f.p.dy - f.u.dy) < 1e-14);
    assert(fabs(f.p.dy - f.v.dy) < 1e-14);

    apply_func(&f.p, zero2d);
    apply_func(&f.u, zero2d);
    apply_func(&f.v, zero2d);

    return f;
}

struct field init_local_acoustic_field(long cells_x, long cells_y, double x[2], double y[2])
{
    struct field f;

    alloc_field(&f.p, cells_x, cells_y + 1); /* differs */
    alloc_field(&f.u, cells_x + 1, cells_y); 
    alloc_field(&f.v, cells_x, cells_y + 1);

    double dx = (x[1] - x[0]) / (double) cells_x;
    double dy = (y[1] - y[0]) / (double) cells_y;
    double p_x[2] = { x[0] + dx / 2.0, x[1] - dx / 2.0 };
    double p_y[2] = { y[0] - dy / 2.0, y[1] - dy / 2.0 }; /* differs */
    double u_x[2] = { x[0], x[1] };
    double u_y[2] = { y[0] + dy / 2.0, y[1] - dy / 2.0 }; 
    double v_x[2] = { x[0] + dx / 2.0, x[1] - dx / 2.0 };
    double v_y[2] = { y[0], y[1] };

    assert(fabs(p_x[0] - (x[0] + dx / 2.0)) < 1e-14);
    assert(fabs(p_y[0] - (y[0] - dy / 2.0)) < 1e-14);
    assert(fabs(u_x[0] - x[0]) < 1e-14);
    assert(fabs(u_y[0] - (y[0] + dy / 2.0)) < 1e-14);
    assert(fabs(v_x[0] - (x[0] + dx / 2.0)) < 1e-14);
    assert(fabs(v_y[0] - y[0]) < 1e-14);

    assert(fabs(p_x[1] - (x[1] - dx / 2.0)) < 1e-14);
    assert(fabs(p_y[1] - (y[1] - dy / 2.0)) < 1e-14);
    assert(fabs(u_x[1] - x[1]) < 1e-14);
    assert(fabs(u_y[1] - (y[1] - dy / 2.0)) < 1e-14);
    assert(fabs(v_x[1] - (x[1] - dx / 2.0)) < 1e-14);
    assert(fabs(v_y[1] - y[1]) < 1e-14);

    set_grid(&f.p, p_x, p_y);
    set_grid(&f.u, u_x, u_y);
    set_grid(&f.v, v_x, v_y);

    assert(fabs(f.p.dx - dx) < 1e-14);
    assert(fabs(f.p.dy - dy) < 1e-14);
    assert(fabs(f.p.dx - f.u.dx) < 1e-14);
    assert(fabs(f.p.dx - f.v.dx) < 1e-14);
    assert(fabs(f.p.dy - f.u.dy) < 1e-14);
    assert(fabs(f.p.dy - f.v.dy) < 1e-14);

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

double assign_to(struct field_variable fv, long i,
                 long j, double value)
{
    fv.value[i + j * fv.size_x] = value;
    return value;
}

double get_from(struct field_variable fv, long i, long j)
{
    return fv.value[i + j * fv.size_x];
}

void set_boundary(struct field *f)
{
    long nx, i, j;

    nx = f->v.size_x;
    j = 0;
    for (i = 0; i < f->v.size_x; ++i)
        f->v.value[i + j*nx] = 0;
    j = f->v.size_y - 1;
    for (i = 0; i < f->v.size_x; ++i)
        f->v.value[i + j*nx] = 0;

    nx = f->u.size_x;
    i = 0;
    for (j = 0; j < f->u.size_y; ++j)
        f->u.value[i + j*nx] = 0;
    i = f->u.size_x - 1;
    for (j = 0; j < f->u.size_y; ++j)
        f->u.value[i + j*nx] = 0;
}

void leapfrog(struct field *f)
{
    long i, j;
    long nx = f->p.size_x;
    long ny = f->p.size_y;
    double dt = f->dt;

    /* used by index macro */
    double *p = f->p.value;
    double *u = f->u.value;
    double *v = f->v.value;

    for (i = 0; i < nx; ++i) {
        for (j = 0; j < ny; ++j) {

            /*if (j == 1) {*/
                /*printf("i=%li:  dt=%e  f.u.dx=%e  f.v.dy=%e  ",i,dt,f->u.dx, f->v.dy);*/
                /*[>printf("U(i+1,j)=%e  U(i,j)=%e  \n",U(i + 1, j), U(i, j));<]*/
                /*printf("V(i,j+1)=%e  V(i,j)=%e  \n",V(i, j + 1), V(i, j));*/
            /*}*/

            P(i, j) +=
#ifdef EXTRA_WORK
                /* NOTE: the pow(exp(-pow(sin(dt),2)),3.2) factor is to make
                 * more operations per memory access. This is to test the
                 * parallel performance.  */
                /*pow(exp(-pow(sin(dt), 2)), 3.2) * */
                /*pow(exp(-pow(sin(dt), 4)), 1.2) * */
                /*pow(exp(-pow(sin(dt), 2)), 4.2) * */
#endif
                dt / f->u.dx * (U(i + 1, j) - U(i, j)) +
                dt / f->v.dy * (V(i, j + 1) - V(i, j));
        }
    }

    /*for (i = 1; i < nx - 1; ++i) {*/
    for (i = 1; i < nx; ++i) {
        for (j = 0; j < ny; ++j) {
            /*if (i == 3) {*/
                /*printf("i=%li:  dt=%e  f.p.dy=%e  ",i,dt,f->p.dx);*/
                /*printf("P(i,j)=%e  P(i - 1,j)=%e  \n",P(i, j), P(i - 1, j));*/
            /*}*/
            U(i, j) += dt / f->p.dx * (P(i, j) - P(i - 1, j));
        }
    }

    for (i = 0; i < nx; ++i) {
        for (j = 1; j < ny; ++j) {

            /*printf("j=%li\n",j);*/
            /*if (j == 3) {*/
                /*printf("i=%li:  dt=%e  f.p.dy=%e  ",i,dt,f->p.dy);*/
                /*printf("P(i,j)=%e  P(i,j - 1)=%e  \n",P(i, j), P(i, j - 1));*/
            /*}*/

            V(i, j) += dt / f->p.dy * (P(i, j) - P(i, j - 1));
        }
    }
}

void timestep_leapfrog(struct field *f, double Nt)
{
    for (long n = 0; n < Nt; ++n)
        leapfrog(f);
}

struct cell_partition *partition_grid(long total_threads,
                                      long cells)
{
    assert(total_threads <= cells);
    long i;
    long cells_per_thread;
    struct cell_partition *partition;

    partition = malloc(sizeof(struct cell_partition) * total_threads);
    cells_per_thread = ceil((double) cells / (double) total_threads);

    for (i = 0; i < total_threads - 1; ++i) {
        partition[i].begin = i * cells_per_thread;
        partition[i].end = (i + 1) * cells_per_thread - 1;
        partition[i].size = partition[i].end - partition[i].begin + 1;
        assert(partition[i].size != 0);
    }

    /* last cell gets is sometimes smaller to make up for the fact that 
     * the number of cells is not evenly divided among the threads. */
    i = total_threads - 1;
    partition[i].begin = (i == 0) ? 0 : partition[i - 1].end + 1;
    partition[i].end = cells - 1;
    partition[i].size = partition[i].end - partition[i].begin + 1;
    assert(partition[i].size != 0);

    return partition;
}

void get_partition_coords(struct cell_partition part, long left, long right, struct field *f, double *y)
{
    long begin = part.begin;
    long end = part.end + 1;
    y[0] = f->v.y[begin];
    y[1] = f->v.y[end];
    /* Due to the staggered grid we need one half cell to the left (bottom) */
    /*if (left != NONE)*/
        /*y[0] -= f->v.dy/2.0;*/
    /*if (right != NONE)*/
        /*y[1] += f->v.dy/2.0;*/
}

void cellindex_to_nodeindex(long tid, struct cell_partition part,
                            long *p0, long *p1,
                            long *u0, long *u1,
                            long *v0, long *v1)
{
    *p0 = part.begin;
    *p1 = part.end + 1;
    *u0 = part.begin;
    *u1 = part.end + 1;
    *v0 = part.begin;
    *v1 = part.end + 1;
    if (tid == 0) {
        ++(*u0);
    }
}

/*void expand_indices(struct cell_partition partition,*/
                    /*long *begin_p, long *end_p,*/
                    /*long *size_p, long *begin_u, long *end_u, long *size_u)*/
/*{*/
    /**begin_p = partition.begin;*/
    /**end_p = partition.end;    [> end index <]*/
    /**begin_u = partition.begin;*/
    /**end_u = partition.end;*/

    /*[> sizes <]*/
    /**size_p = *end_p - *begin_p + 1;*/
    /**size_u = *end_u - *begin_u + 1;*/
    /**size_v = *end_v - *begin_v + 1;*/
/*}*/

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

void parse_cmdline(long *nx,
                   long *threads, char *outfile,
                   int argc, char *argv[])
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
    for (long i = 0; i < f.size_x; ++i) {
        for (long j = 0; j < f.size_y; ++j) {
            fprintf(fp, "%.16e\t%.16e\t%.16e\n", f.x[i], f.y[j], get_from(f, i, j));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return EXIT_SUCCESS;
}

double gauss(double x)
{
    return exp(-pow(x - 0.5, 2) / pow(0.05, 2));
}

double gauss2d(double x, double y)
{
    /*double r_squared = pow(x - 0.2, 2) + pow(y - 0.5, 2);*/
    double r_squared = pow(x - 0.5, 2) + pow(y - 0.5, 2);
    return exp(-r_squared / pow(0.1, 2));
}

double zero(double x)
{
    (void) x;
    return 0.0;
}

double zero2d(double x, double y)
{
    (void) x;
    (void) y;
    return 0.0;
}

double one2d(double x, double y)
{
    (void) x;
    (void) y;
    return 1.0;
}

double identity(double x)
{
    return x;
}

double identity2d(double x, double y)
{
    (void) y;
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
