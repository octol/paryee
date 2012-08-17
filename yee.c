/* 
 * Basic implementation of 1D wave equation on system form
 *
 *      p_t = a u_x
 *      u_t = b p_x
 *
 * This formulation is sometimes used in acoustics, where then p is the
 * pressure and u is the velocity field in the direction of the x-axis.
 * 
 * We use a staggared grid with u defined on the endpoints, and thus on full
 * integer indices. N denotes the number of interval lengths, i.e., we have
 * N+1 number of nodes for u
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <pthread.h>
#include <time.h>

#ifndef NX
#define NX 1024*2*2*2
#endif

#ifndef NODES
#define NODES 4
#endif

/* 
 * Data structures
 */
typedef struct FieldVariable FieldVariable;
typedef struct Field Field;
typedef struct FieldUpdateParam FieldUpdateParam;

struct FieldVariable {
    double* value;
    double* x;
    double dx;
    unsigned long N;    /* number of points */
};

struct Field {
    FieldVariable p;
    FieldVariable u;
    double dt;
    double Nt;
};

struct FieldUpdateParam {
    FieldVariable* dst;
    int dst1;
    int dst2;
    FieldVariable* src;
    int src1;
    double dt;
};

void alloc_field (FieldVariable* f, unsigned long N)
{
    f->N = N;
    f->value = (double*) malloc (sizeof(double)*f->N);
    f->x = (double*) malloc (sizeof(double)*f->N);
}

void free_field (FieldVariable f)
{
    free(f.value);
    free(f.x);
}

/* 
 * Generate grid from start to end with N number of points  
 */
void set_grid (FieldVariable* f, double start, double end)
{
    /* NOTE: f->N is the number of grid points, not intervals */
    f->dx = (end-start)/(f->N - 1);
    for (unsigned int i=0; i<f->N; ++i) 
        f->x[i] = start + i*f->dx;
}

/*
 * Vectorization of function.
 * Apply the function func on the array dst with argument arg.
 */
void vec_func (double* dst, double* arg, double (*func) (double), 
               unsigned long i1, unsigned long i2)
{
    for (unsigned long i=i1; i<i2; ++i)
        dst[i] = func (arg[i]);
}

/*
 * Vectorization of function applied to FieldVariable.
 */
void field_func (FieldVariable* f, double (*func) (double))
{
    vec_func (f->value, f->x, func, 0, f->N);
}

/* 
 * Some grid functions 
 */
double gauss (double x) { return exp(-pow(x-0.5,2)/pow(0.05,2)); }
double zero (double x) { return 0.0; }

int write_to_disk (FieldVariable f, char* str)
{
    /* append file suffix */
    int n = strlen(str) + 4;
    char fstr[n];
    strcpy(fstr,str);
    strcat(fstr,".tsv");

    FILE* fp;
    fp = fopen(fstr,"w");
    if (fp == NULL) {
        perror ("Error: can not write to disk");
        return EXIT_FAILURE;
    } 
    for (unsigned int i=0; i<f.N; ++i)
        fprintf (fp, "%e\t%e\n", f.x[i], f.value[i]);
    fclose (fp);
    return EXIT_SUCCESS;
}

/*
 * Leapfrog time update
 * Update the field points dst for the indices dst1,...,dst2, using the
 * field points src using indices src1,...,src1 + (dst2-dst1)
 */
void update_field (FieldVariable* dst, int dst1, int dst2, 
                   FieldVariable* src, int src1, double dt)
{
    /* 
     * NOTE: the exp(-pow(sin(dt),)) factor is to make more operations per
     * memory access. This is to test the parallell performance. 
     */
    for (int i=dst1; i<dst2; ++i, ++src1)
        dst->value[i] += dt*exp(-pow(sin(dt),2))*
            (src->value[src1+1] - src->value[src1])/src->dx;
}

/* 
 * Leapfrog time update adapted to use return value and argument suitable
 * for pthreads.
 */
void* pthread_update_field (void* arg)
{
    FieldUpdateParam* p = (FieldUpdateParam*) arg;
    update_field (p->dst, p->dst1, p->dst2, p->src, p->src1, p->dt);
    return NULL;
}

/* 
 * Collect parameters
 */
FieldUpdateParam collect_param (FieldVariable* dst, int dst1, int dst2, 
                                FieldVariable* src, int src1, double dt)
{
    FieldUpdateParam param;
    param.dst = dst;
    param.dst1 = dst1;
    param.dst2 = dst2;
    param.src = src;
    param.src1 = src1;
    param.dt = dt;
    return param;
}

/*
 * Round up integer division
 */
int round_up_divide (int x, int y)
{
    return (x-1)/y + 1;
}

/* 
 * Main program 
 */
int main()
{
    /* Paramters */
    double length=1, cfl=1, T=0.2, c=1;

    /* Initialize */
    Field f;
    alloc_field(&f.p, NX);
    alloc_field(&f.u, NX+1);
    set_grid (&f.p, 0.5*length/NX, length-0.5*length/NX);
    set_grid (&f.u, 0, length);
    field_func (&f.p, gauss);
    field_func (&f.u, zero);

    /* Partition data */
    int cells_per_node = NX/NODES;
    assert(cells_per_node*NODES == NX);
    int start, end;
    pthread_t thr[NODES];
    FieldUpdateParam param[NODES];

    /* Timestep */
    f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
    f.Nt = T/f.dt;
    clock_t tic, toc;

    tic = clock();
    for (int n=0; n<f.Nt; ++n) {

        /* update the pressure (p) */
        for (int i=0; i<NODES; ++i) {
            start = i*cells_per_node;
            end = (i+1)*cells_per_node;
            param[i] = collect_param (&f.p, start, end, &f.u, start, f.dt);
            pthread_create (&thr[i], NULL, pthread_update_field, &param[i]);
            /*update_field (&f.p, start, end, &f.u, start, f.dt);*/
        }
        for (int i=0; i<NODES; ++i)
            pthread_join (thr[i], NULL);

        /* update the velocity (u) */
        for (int i=0; i<NODES; ++i) {
            start = 1 + i*cells_per_node;
            end = 1 + (i+1)*cells_per_node;
            if (i==NODES-1) /* last node gets one less */
                --end;
            param[i] = collect_param (&f.u, start, end, &f.p, start-1,f.dt);
            pthread_create (&thr[i], NULL, pthread_update_field, &param[i]);
            /*update_field (&f.u, start, end, &f.p, start-1, f.dt);*/
        }
        for (int i=0; i<NODES; ++i)
            pthread_join (thr[i], NULL);
        assert( (NODES-1 + 1)*cells_per_node == NX );
    }
    toc = clock();
    printf ("Elapsed: %f seconds\n", (double)(toc-tic)/CLOCKS_PER_SEC);

    /* write data to disk */
    write_to_disk(f.p, "output_p"); 
    write_to_disk(f.u, "output_u"); 

    free_field(f.p);
    free_field(f.u);

    return EXIT_SUCCESS;
}
