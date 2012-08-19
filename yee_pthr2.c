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
#include <getopt.h>

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

/* 
 * Generate grid from start to end with N number of points  
 */
void set_grid (FieldVariable* f, double start, double end)
{
    unsigned int i;

    /* NOTE: f->N is the number of grid points, not intervals */
    f->dx = (end-start)/(f->N - 1);
    for (i=0; i<f->N; ++i) 
        f->x[i] = start + i*f->dx;
}

/*
 * Vectorization of function.
 * Apply the function func on the array dst with argument arg.
 */
void vec_func (double* dst, double* arg, double (*func) (double), 
               unsigned long i1, unsigned long i2)
{
    unsigned long i;
    for (i=i1; i<i2; ++i)
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
    for (i=0; i<f.N; ++i)
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
     * NOTE: the pow(exp(-pow(sin(dt),2)),3.2) factor is to make more
     * operations per memory access. This is to test the parallell
     * performance. 
     */
    int i;
    for (i=dst1; i<dst2; ++i, ++src1)
        dst->value[i] += dt*
            /*pow(exp(-pow(sin(dt),2)),3.2)**/
            /*pow(exp(-pow(sin(dt),4)),1.2)**/
            /*pow(exp(-pow(sin(dt),2)),4.2)**/
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
    pthread_exit(NULL);
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
int main(int argc, char* argv[])
{
    double length, cfl, T, c;
    int NX, NODES, opt;
    int cells_per_node;
    int n, i, start, end;
    /*clock_t tic, toc;*/
    Field f;
    pthread_t* thr;
    pthread_attr_t attr;
    FieldUpdateParam* param_p;
    FieldUpdateParam* param_u;

    /* Parse parameters */
    NX=1024; 
    NODES=2;
    while ((opt = getopt(argc, argv, "t:n:")) != -1) {
        switch (opt) {
        case 't':
            NODES = atoi (optarg);
            break;
        case 'n':
            NX = atoi(optarg);
            break;
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-t threads] [-n intervals] name\n",
                    argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    printf("Running with: N=%d, threads=%d\n", NX, NODES);

    /* Numerical paramters */
    length = 1;
    cfl = 1; 
    T = 0.3; 
    c = 1;

    /* Initialize */
    alloc_field(&f.p, NX);
    alloc_field(&f.u, NX+1);
    set_grid (&f.p, 0.5*length/NX, length-0.5*length/NX);
    set_grid (&f.u, 0, length);
    field_func (&f.p, gauss);
    field_func (&f.u, zero);

    /* Depends on the numerical variables initialzed above */
    f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
    f.Nt = T/f.dt;

    /* Thread - partition data */
    pthread_attr_init (&attr);
    pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
    thr = malloc (sizeof(pthread_t)*NODES);
    param_p = malloc (sizeof(FieldUpdateParam)*NODES);
    param_u = malloc (sizeof(FieldUpdateParam)*NODES);
    if (!thr || !param_p || !param_u) {
        fprintf(stderr,"Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    cells_per_node = NX/NODES;
    assert(cells_per_node*NODES == NX);

    for (i=0; i<NODES; ++i) {
        start = i*cells_per_node;
        end = (i+1)*cells_per_node;
        param_p[i] = collect_param (&f.p, start, end, &f.u, start, f.dt);

        start = 1 + i*cells_per_node;
        end = 1 + (i+1)*cells_per_node;
        if (i==NODES-1) /* last node gets one less */
            --end;
        param_u[i] = collect_param (&f.u, start, end, &f.p, start-1,f.dt);
    }

    /* Timestep */
    /*tic = clock();*/
    for (n=0; n<f.Nt; ++n) {
        /* update the pressure (p) */
        for (i=0; i<NODES; ++i) 
            pthread_create (&thr[i], &attr, pthread_update_field, &param_p[i]);
        for (i=0; i<NODES; ++i)
            pthread_join (thr[i], NULL);

        /* update the velocity (u) */
        for (i=0; i<NODES; ++i)
            pthread_create (&thr[i], &attr, pthread_update_field, &param_u[i]);
        for (i=0; i<NODES; ++i)
            pthread_join (thr[i], NULL);
    }
    /*toc = clock();*/
    /*printf ("Elapsed: %f seconds\n", (double)(toc-tic)/CLOCKS_PER_SEC);*/

    /* write data to disk and free data */
    write_to_disk(f.p, "output_p"); 
    write_to_disk(f.u, "output_u"); 
    free_field(f.p);
    free_field(f.u);

    return EXIT_SUCCESS;
}
