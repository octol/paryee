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

void alloc_field (struct FieldVariable* f, unsigned long size)
{
    f->size = size;
    f->value = (double*) malloc (sizeof(double)*size);
    f->x = (double*) malloc (sizeof(double)*size);

    if (!f->value || !f->x) {
        fprintf(stderr,"Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

void free_field (struct FieldVariable f)
{
    free(f.value);
    free(f.x);
}

void set_grid (struct FieldVariable* f, double start, double end)
{
    unsigned int i;

    /* NOTE: f->size is the number of grid points, not intervals */
    f->dx = (end-start)/(f->size - 1);
    for (i=0; i<f->size; ++i) 
        f->x[i] = start + i*f->dx;
}

void vec_func (double* dst, double* arg, double (*func) (double), 
               unsigned long i1, unsigned long i2)
{
    unsigned long i;
    for (i=i1; i<i2; ++i)
        dst[i] = func (arg[i]);
}

void apply_func (struct FieldVariable* f, double (*func) (double))
{
    vec_func (f->value, f->x, func, 0, f->size);
}

struct UpdateParam collect_param (
                            struct FieldVariable* dst, int dst1, int dst2, 
                            struct FieldVariable* src, int src1, double dt)
{
    struct UpdateParam param;
    param.dst = dst; param.dst1 = dst1; param.dst2 = dst2; 
    param.src = src; param.src1 = src1; 
    param.dt = dt;
    return param;
}

void parse_cmdline (unsigned long* nx, unsigned long* nodes, 
                    int argc, char* argv[])
{
    int opt;
    while ((opt = getopt(argc, argv, "t:n:")) != -1) {
        switch (opt) {
            case 't':
                *nodes = atoi (optarg);
                break;
            case 'n':
                *nx = atoi(optarg);
                break;
            default: /* '?' */
                fprintf(stderr, 
                        "Usage: %s [-t threads] [-n intervals]\n",
                        argv[0]);
                exit(EXIT_FAILURE);
        }
    }
}

void update_field_s (struct FieldVariable* dst, int idst, int size, 
                     struct FieldVariable* src, int isrc, double dt)
{
    int i;
    for (i=idst; i<idst+size; ++i, ++isrc)
        dst->value[i] += dt*
#ifdef EXTRA_WORK
            /* 
             * NOTE: the pow(exp(-pow(sin(dt),2)),3.2) factor is to make
             * more operations per memory access. This is to test the
             * parallel performance. 
             */
            pow(exp(-pow(sin(dt),2)),3.2)*
            pow(exp(-pow(sin(dt),4)),1.2)*
            pow(exp(-pow(sin(dt),2)),4.2)*
#endif
            (src->value[isrc+1] - src->value[isrc])/src->dx;
}

void update_field_i (struct FieldVariable* dst, int dst1, int dst2, 
                     struct FieldVariable* src, int src1, double dt)
{
    int i;
    for (i=dst1; i<=dst2; ++i, ++src1)
        dst->value[i] += dt*
#ifdef EXTRA_WORK
            /* 
             * NOTE: the pow(exp(-pow(sin(dt),2)),3.2) factor is to make
             * more operations per memory access. This is to test the
             * parallel performance. 
             */
            pow(exp(-pow(sin(dt),2)),3.2)*
            pow(exp(-pow(sin(dt),4)),1.2)*
            pow(exp(-pow(sin(dt),2)),4.2)*
#endif
            (src->value[src1+1] - src->value[src1])/src->dx;
}

struct Partition partition_grid (int current_node,int nodes, 
                                 int cells_per_node)
{
    struct Partition partition;

    partition.p[0] = current_node*cells_per_node;
    partition.p[1] = (current_node+1)*cells_per_node - 1;
    partition.u[0] = current_node*cells_per_node;
    partition.u[1] = (current_node+1)*cells_per_node - 1;

    /* First node skips the first u, as it sits on the boundary */
    if (current_node==0) 
        ++partition.u[0];

    return partition;
}

void expand_indices (struct Partition partition, 
                     int* begin_p, int* end_p, int* size_p, 
                     int* begin_u, int* end_u, int* size_u)
{
    *begin_p = partition.p[0]; /* begin index */
    *begin_u = partition.u[0]; 
    *end_p = partition.p[1]; /* end index */
    *end_u = partition.u[1]; 

    /* sizes */
    *size_p = *end_p - *begin_p + 1;
    *size_u = *end_u - *begin_u + 1;  
}

void verify_grid_integrity (struct Partition partition, int tid, 
                            unsigned long nx, int numworkers, int left)
{
    int bp, ep, sp, bu, eu, su;
    expand_indices (partition, &bp, &ep, &sp, &bu, &eu, &su);

    /* Sanity checks */
    assert (bp == (tid-1)*nx/numworkers);
    assert (ep == (tid)*nx/numworkers-1);
    assert (sp == nx/numworkers);
    if (left == NONE) {
        assert (bu == 1);
        assert (su == nx/numworkers-1);
    } else {
        assert (bu == (tid-1)*nx/numworkers);
        assert (su == nx/numworkers);
    }
    assert (eu == (tid)*nx/numworkers-1);
    
    /* Specific sanity checks */
    if (nx==8 && numworkers==4) {
        switch (tid) {
        case 1:
            assert (bp==0 && ep==1 && sp==2);
            assert (bu==1 && eu==1 && su==1);
            break; 
        case 2:
            assert (bp==2 && ep==3 && sp==2);
            assert (bu==2 && eu==3 && su==2);
            break; 
        case 3:
            assert (bp==4 && ep==5 && sp==2);
            assert (bu==4 && eu==5 && su==2);
            break;
        case 4:
            assert (bp==6 && ep==7 && sp==2);
            assert (bu==6 && eu==7 && su==2);
            break;
        }
    }
}

void set_local_index (int size_p, int size_u, int left, 
                int* local_begin_p, int* local_end_p, int* local_size_p, 
                int* local_begin_u, int* local_end_u, int* local_size_u)
{
        *local_begin_p = 1; /* since p is padded by one point to the left */
        *local_begin_u = 0; /* since u is not padded on the left */
        if (left == NONE) {
            *local_begin_p = 0;
            *local_begin_u = 1;
        }
        *local_end_p = *local_begin_p + size_p - 1;
        *local_end_u = *local_begin_u + size_u - 1;

        *local_size_p = *local_end_p + 1;
        *local_size_u = *local_end_u + 2; 
}

int write_to_disk (struct FieldVariable f, char* str)
{
    /* append file suffix */
    char fstr[256];
    FILE* fp;
    unsigned int i;
    strcpy(fstr,str);
    strcat(fstr,".tsv");

    fp = fopen(fstr,"w");
    if (fp == NULL) {
        perror ("Error: can not write to disk");
        return EXIT_FAILURE;
    } 
    for (i=0; i<f.size; ++i)
        fprintf (fp, "%e\t%e\n", f.x[i], f.value[i]);
    fclose (fp);
    return EXIT_SUCCESS;
}

double gauss (double x) 
{ 
    return exp(-pow(x-0.5,2)/pow(0.05,2)); 
}

double zero (double x) 
{ 
    return 0.0; 
}

int round_up_divide (int x, int y)
{
    return (x-1)/y + 1;
}

double gettime (void) 
{
    struct timeval t;
    gettimeofday (&t, NULL);
    return (1.0e-6*t.tv_usec + t.tv_sec);
}

