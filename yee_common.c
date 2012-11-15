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

void alloc_field (FieldVariable* f, unsigned long size)
{
    f->size = size;
    f->value = (double*) malloc (sizeof(double)*size);
    f->x = (double*) malloc (sizeof(double)*size);

    if (!f->value || !f->x) {
        fprintf(stderr,"Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

void free_field (FieldVariable f)
{
    free(f.value);
    free(f.x);
}

void set_grid (FieldVariable* f, double start, double end)
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

void field_func (FieldVariable* f, double (*func) (double))
{
    vec_func (f->value, f->x, func, 0, f->size);
}

UpdateParam collect_param (FieldVariable* dst, int dst1, int dst2, 
                           FieldVariable* src, int src1, double dt)
{
    UpdateParam param;
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

void update_field (FieldVariable* dst, int dst1, int dst2, 
                   FieldVariable* src, int src1, double dt)
{
    int i;
    for (i=dst1; i<dst2; ++i, ++src1)
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

void update_field2 (FieldVariable* dst, int dst1, int dst2, 
                    FieldVariable* src, int src1, double dt)
{
    int i;
    for (i=dst1; i<dst2; ++i, ++src1)
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

FieldPartition partition_grid(int current_node,int nodes,int cells_per_node)
{
    FieldPartition part;

    part.start_p = current_node*cells_per_node;
    part.end_p = (current_node+1)*cells_per_node;
    part.start_u = 1 + current_node*cells_per_node;
    part.end_u = 1 + (current_node+1)*cells_per_node;
    if (current_node==nodes-1) /* last node gets one less */
        --part.end_u;

    return part;
}

FieldPartition partition_grid_into_cells(int current_node,int nodes,
                                         int cells_per_node)
{
    FieldPartition part;

    part.start_p = current_node*cells_per_node;
    part.end_p = (current_node+1)*cells_per_node - 1;
    part.start_u = current_node*cells_per_node;
    part.end_u = (current_node+1)*cells_per_node - 1;

    /*int i = current_node;                                   */
    /*int bp = part.start_p;                                  */
    /*int ep = part.end_p;                                    */
    /*int bu = part.start_u;                                  */
    /*int eu = part.end_u;                                    */
    /*printf ("Partition %i: cells_per_node=%d  bp=%d  ep=%d",*/
    /*        i,cells_per_node,bp,ep);                        */
    /*printf ("  bu=%d  eu=%d\n",bu,eu);                      */

    /* First node skips the first u, as it sits on the boundary */
    if (current_node==0) 
        ++part.start_u;

    return part;
}

int write_to_disk (FieldVariable f, char* str)
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

