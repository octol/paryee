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

void alloc_field(struct field_variable *f, long size)
{
    long block = sizeof(double) * size;
    f->size = size;
    f->value = (double *) malloc(block);
    f->x = (double *) malloc(block);

    if (!f->value || !f->x) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

void free_field(struct field_variable f)
{
    free(f.value);
    free(f.x);
}

void set_grid(struct field_variable *f, double start, double end)
{
    long i;

    /* NOTE: f->size is the number of grid points, not intervals */
    f->dx = (end - start) / (f->size - 1);
    for (i = 0; i < f->size; ++i)
        f->x[i] = start + i * f->dx;
}

void vec_func(double *dst, double *arg, double (*func) (double),
              long i1, long i2)
{
    long i;
    for (i = i1; i < i2; ++i)
        dst[i] = func(arg[i]);
}

void apply_func(struct field_variable *f, double (*func) (double))
{
    vec_func(f->value, f->x, func, 0, f->size);
}

void parse_cmdline(long *nx, long *threads,
                   char *outfile_p, char *outfile_u,
                   int argc, char *argv[])
{
    int opt;
    while ((opt = getopt(argc, argv, "t:n:p:u:")) != -1) {
        switch (opt) {
        case 't':
            *threads = atoi(optarg);
            break;
        case 'n':
            *nx = atoi(optarg);
            break;
        case 'p':
            strncpy(outfile_p, optarg, STR_SIZE);
            outfile_p[STR_SIZE - 1] = '\0';     /* force null termination */
            break;
        case 'u':
            strncpy(outfile_u, optarg, STR_SIZE);
            outfile_u[STR_SIZE - 1] = '\0';     /* force null termination */
            break;
        default:               /* '?' */
            fprintf(stderr, "Usage: %s ", argv[0]);
            if (threads != NULL)
                fprintf(stderr, "[-t threads] ");
            fprintf(stderr, "[-n intervals] ");
            fprintf(stderr, "[-p outfile_p] [-u outfile_u]\n");
            exit(EXIT_FAILURE);
        }
    }
}

void update_field_s(struct field_variable *restrict dst, int idst,
                    int size, struct field_variable *restrict src,
                    int isrc, double dt)
{
    int i;
    for (i = idst; i < idst + size; ++i, ++isrc)
        dst->value[i] += dt *
#ifdef EXTRA_WORK
            /* 
             * NOTE: the pow(exp(-pow(sin(dt),2)),3.2) factor is to make
             * more operations per memory access. This is to test the
             * parallel performance. 
             */
            pow(exp(-pow(sin(dt), 2)), 3.2) *
            pow(exp(-pow(sin(dt), 4)), 1.2) *
            pow(exp(-pow(sin(dt), 2)), 4.2) *
#endif
            (src->value[isrc + 1] - src->value[isrc]) / src->dx;
}

void update_field_i(struct field_variable *restrict dst, int dst1,
                    int dst2, struct field_variable *restrict src,
                    int src1, double dt)
{
    int i;
    for (i = dst1; i <= dst2; ++i, ++src1)
        dst->value[i] += dt *
#ifdef EXTRA_WORK
            /* 
             * NOTE: the pow(exp(-pow(sin(dt),2)),3.2) factor is to make
             * more operations per memory access. This is to test the
             * parallel performance. 
             */
            pow(exp(-pow(sin(dt), 2)), 3.2) *
            pow(exp(-pow(sin(dt), 4)), 1.2) *
            pow(exp(-pow(sin(dt), 2)), 4.2) *
#endif
            (src->value[src1 + 1] - src->value[src1]) / src->dx;
}

struct partition partition_grid(int current_thread, int cells_per_thread)
{
    struct partition partition;

    partition.p[0] = current_thread * cells_per_thread;
    partition.p[1] = (current_thread + 1) * cells_per_thread - 1;
    partition.u[0] = current_thread * cells_per_thread;
    partition.u[1] = (current_thread + 1) * cells_per_thread - 1;

    /* First thread skips the first u, as it sits on the boundary */
    if (current_thread == 0)
        ++partition.u[0];

    return partition;
}

void expand_indices(struct partition partition,
                    long *begin_p, long *end_p,
                    long *size_p, long *begin_u,
                    long *end_u, long *size_u)
{
    *begin_p = partition.p[0];  /* begin index */
    *begin_u = partition.u[0];
    *end_p = partition.p[1];    /* end index */
    *end_u = partition.u[1];

    /* sizes */
    *size_p = *end_p - *begin_p + 1;
    *size_u = *end_u - *begin_u + 1;
}

void verify_grid_integrity(struct partition partition, int tid,
                           long nx, int numworkers, int left)
{
    long bp, ep, sp, bu, eu, su;
    expand_indices(partition, &bp, &ep, &sp, &bu, &eu, &su);

    /* Sanity checks */
    assert(bp == (tid - 1) * nx / numworkers);
    assert(ep == (tid) * nx / numworkers - 1);
    assert(sp == nx / numworkers);
    if (left == NONE) {
        assert(bu == 1);
        assert(su == nx / numworkers - 1);
    } else {
        assert(bu == (tid - 1) * nx / numworkers);
        assert(su == nx / numworkers);
    }
    assert(eu == (tid) * nx / numworkers - 1);

    /* Specific sanity checks */
    if (nx == 8 && numworkers == 4) {
        switch (tid) {
        case 1:
            assert(bp == 0 && ep == 1 && sp == 2);
            assert(bu == 1 && eu == 1 && su == 1);
            break;
        case 2:
            assert(bp == 2 && ep == 3 && sp == 2);
            assert(bu == 2 && eu == 3 && su == 2);
            break;
        case 3:
            assert(bp == 4 && ep == 5 && sp == 2);
            assert(bu == 4 && eu == 5 && su == 2);
            break;
        case 4:
            assert(bp == 6 && ep == 7 && sp == 2);
            assert(bu == 6 && eu == 7 && su == 2);
            break;
        }
    }
}

void set_local_index(long size_p, long size_u,
                     long left, long *local_begin_p,
                     long *local_end_p,
                     long *local_size_p,
                     long *local_begin_u,
                     long *local_end_u,
                     long *local_size_u)
{
    /* since p is padded by one point to the left */
    *local_begin_p = 1;
    /* since u is not padded on the left */
    *local_begin_u = 0;
    if (left == NONE) {
        *local_begin_p = 0;
        *local_begin_u = 1;
    }
    *local_end_p = *local_begin_p + size_p - 1;
    *local_end_u = *local_begin_u + size_u - 1;

    *local_size_p = *local_end_p + 1;
    *local_size_u = *local_end_u + 2;
}

int write_to_disk(struct field_variable f, char *str)
{
    /* append file suffix */
    char fstr[256];
    FILE *fp;
    unsigned int i;
    strcpy(fstr, str);
    /*strcat(fstr,".tsv"); */

    printf("Writing to: %s\n", fstr);

    fp = fopen(fstr, "w");
    if (fp == NULL) {
        perror("Error: can not write to disk");
        return EXIT_FAILURE;
    }
    for (i = 0; i < f.size; ++i)
        fprintf(fp, "%e\t%e\n", f.x[i], f.value[i]);
    fclose(fp);
    return EXIT_SUCCESS;
}

double gauss(double x)
{
    return exp(-pow(x - 0.5, 2) / pow(0.05, 2));
}

double zero(double x)
{
    return 0.0 * x;
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
