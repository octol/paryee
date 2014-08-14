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

void py_alloc_field(struct py_field_variable *f, const long size_x,
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

void py_free_field(struct py_field_variable f)
{
    free(f.value);
    free(f.x);
    free(f.y);
}

void py_set_grid(struct py_field_variable *f, const double x[2],
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

void py_vec_func(double *dst, double (*func) (double), const double *arg,
              const long size)
{
    for (long i = 0; i < size; ++i)
        dst[i] = func(arg[i]);
}

void py_vec_func2d(double *dst, double (*func) (double, double),
                const double *arg_x, const long size_x,
                const double *arg_y, const long size_y)
{
    for (long ii = 0; ii < size_x; ++ii)
        for (long jj = 0; jj < size_y; ++jj)
            dst[ii + jj * size_x] = func(arg_x[ii], arg_y[jj]);
}

void py_apply_func(struct py_field_variable *f, double (*func) (double, double))
{
    py_vec_func2d(f->value, func, f->x, f->size_x, f->y, f->size_y);
}

struct py_field py_init_acoustic_field(long cells_x,
                                 long cells_y, double x[2], double y[2])
{
    struct py_field f;

    py_alloc_field(&f.p, cells_x, cells_y);
    py_alloc_field(&f.u, cells_x + 1, cells_y);
    py_alloc_field(&f.v, cells_x, cells_y + 1);

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

    py_set_grid(&f.p, p_x, p_y);
    py_set_grid(&f.u, u_x, u_y);
    py_set_grid(&f.v, v_x, v_y);

    assert(fabs(f.p.dx - f.u.dx) < 1e-14);
    assert(fabs(f.p.dx - f.v.dx) < 1e-14);
    assert(fabs(f.p.dy - f.u.dy) < 1e-14);
    assert(fabs(f.p.dy - f.v.dy) < 1e-14);

    py_apply_func(&f.p, py_zero2d);
    py_apply_func(&f.u, py_zero2d);
    py_apply_func(&f.v, py_zero2d);

    return f;
}

struct py_field py_init_local_acoustic_field(long cells_x, long cells_y,
                                       double x[2], double y[2])
{
    struct py_field f;

    py_alloc_field(&f.p, cells_x, cells_y + 1);    /* differs */
    py_alloc_field(&f.u, cells_x + 1, cells_y);
    py_alloc_field(&f.v, cells_x, cells_y + 1);

    double dx = (x[1] - x[0]) / (double) cells_x;
    double dy = (y[1] - y[0]) / (double) cells_y;
    double p_x[2] = { x[0] + dx / 2.0, x[1] - dx / 2.0 };
    double p_y[2] = { y[0] - dy / 2.0, y[1] - dy / 2.0 };       /* differs */
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

    py_set_grid(&f.p, p_x, p_y);
    py_set_grid(&f.u, u_x, u_y);
    py_set_grid(&f.v, v_x, v_y);

    assert(fabs(f.p.dx - dx) < 1e-14);
    assert(fabs(f.p.dy - dy) < 1e-14);
    assert(fabs(f.p.dx - f.u.dx) < 1e-14);
    assert(fabs(f.p.dx - f.v.dx) < 1e-14);
    assert(fabs(f.p.dy - f.u.dy) < 1e-14);
    assert(fabs(f.p.dy - f.v.dy) < 1e-14);

    py_apply_func(&f.p, py_zero2d);
    py_apply_func(&f.u, py_zero2d);
    py_apply_func(&f.v, py_zero2d);

    return f;
}

void py_free_acoustic_field(struct py_field f)
{
    py_free_field(f.p);
    py_free_field(f.u);
    py_free_field(f.v);
}

double py_assign_to(struct py_field_variable fv, long i, long j, double value)
{
    fv.value[i + j * fv.size_x] = value;
    return value;
}

double py_get_from(struct py_field_variable fv, long i, long j)
{
    return fv.value[i + j * fv.size_x];
}

void py_set_boundary(struct py_field *f)
{
    long nx, i, j;

    nx = f->v.size_x;
    j = 0;
    for (i = 0; i < f->v.size_x; ++i)
        f->v.value[i + j * nx] = 0;
    j = f->v.size_y - 1;
    for (i = 0; i < f->v.size_x; ++i)
        f->v.value[i + j * nx] = 0;

    nx = f->u.size_x;
    i = 0;
    for (j = 0; j < f->u.size_y; ++j)
        f->u.value[i + j * nx] = 0;
    i = f->u.size_x - 1;
    for (j = 0; j < f->u.size_y; ++j)
        f->u.value[i + j * nx] = 0;
}

void py_leapfrog(struct py_field *f)
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

            P(i, j) +=
#ifdef EXTRA_WORK
                /* NOTE: the pow(exp(-pow(sin(dt),2)),3.2) factor is to make
                 * more operations per memory access. This is to test the
                 * parallel performance.  */
                pow(exp(-pow(sin(dt), 2)), 3.2) *
                pow(exp(-pow(sin(dt), 4)), 1.2) *
                pow(exp(-pow(sin(dt), 2)), 4.2) *
#endif
                dt / f->u.dx * (U(i + 1, j) - U(i, j)) +
                dt / f->v.dy * (V(i, j + 1) - V(i, j));
        }
    }

    for (i = 1; i < nx; ++i)
        for (j = 0; j < ny; ++j)
            U(i, j) += dt / f->p.dx * (P(i, j) - P(i - 1, j));

    for (i = 0; i < nx; ++i)
        for (j = 1; j < ny; ++j)
            V(i, j) += dt / f->p.dy * (P(i, j) - P(i, j - 1));
}

void py_timestep_leapfrog(struct py_field *f, double Nt)
{
    for (long n = 0; n < Nt; ++n)
        py_leapfrog(f);
}

struct py_cell_partition *py_partition_grid(long total_threads, long cells)
{
    assert(total_threads <= cells);
    long i;
    long cells_per_thread;
    struct py_cell_partition *partition;

    partition = malloc(sizeof(struct py_cell_partition) * total_threads);
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

void py_get_partition_coords(struct py_cell_partition part, struct py_field *f,
                          double *y)
{
    long begin = part.begin;
    long end = part.end + 1;
    y[0] = f->v.y[begin];
    y[1] = f->v.y[end];
}

void py_cellindex_to_nodeindex(long tid, struct py_cell_partition part,
                            long *p0, long *p1,
                            long *u0, long *u1, long *v0, long *v1)
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

void py_parse_cmdline(long *nx, long *threads, char *outfile, int *write,
                   int argc, char *argv[])
{
    int opt;
    while ((opt = getopt(argc, argv, "t:n:o:wq")) != -1) {
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
        case 'q':
            *write = 0;
            break;
        default:               /* '?' */
            fprintf(stderr, "Usage: %s ", argv[0]);
            if (threads != NULL)
                fprintf(stderr, "[-t threads] ");
            fprintf(stderr, "[-n intervals] ");
            fprintf(stderr, "[-o outfile]");
            fprintf(stderr, "[-q]\n");
            exit(EXIT_FAILURE);
        }
    }
}

int py_write_to_disk(struct py_field_variable f, char *fstr)
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
            fprintf(fp, "%.16e\t%.16e\t%.16e\n", f.x[i], f.y[j],
                    py_get_from(f, i, j));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return EXIT_SUCCESS;
}

double py_gauss(double x)
{
    return exp(-pow(x - 0.5, 2) / pow(0.05, 2));
}

double py_gauss2d(double x, double y)
{
    double r_squared = pow(x - 0.2, 2) + pow(y - 0.5, 2);
    return exp(-r_squared / pow(0.1, 2));
}

double py_zero(double x)
{
    (void) x;
    return 0.0;
}

double py_zero2d(double x, double y)
{
    (void) x;
    (void) y;
    return 0.0;
}

double py_one2d(double x, double y)
{
    (void) x;
    (void) y;
    return 1.0;
}

double py_identity(double x)
{
    return x;
}

double py_identity2d(double x, double y)
{
    (void) y;
    return x;
}

int py_round_up_divide(int x, int y)
{
    return (x - 1) / y + 1;
}

double py_gettime(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1.0e-6 * t.tv_usec + t.tv_sec);
}
