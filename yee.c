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

#ifndef NX
#define NX 1024
#endif

#ifndef NODES
#define NODES 4
#endif

/* 
 * Data structures
 */
struct FieldVariable {
    double* value;
    double* x;
    double dx;
    unsigned long N;    /* number of points */
};

struct Field {
    struct FieldVariable p;
    struct FieldVariable u;
    double dt;
    double Nt;
};

typedef struct FieldVariable FieldVariable_t;
typedef struct Field Field_t;

void AllocField (struct FieldVariable* f, unsigned long N)
{
    f->N = N;
    f->value = (double*) malloc (sizeof(double)*f->N);
    f->x = (double*) malloc (sizeof(double)*f->N);
}

void FreeField (struct FieldVariable f)
{
    free(f.value);
    free(f.x);
}

/* 
 * Generate grid from start to end with N number of points  
 */
void SetGrid (struct FieldVariable* f, double start, double end)
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
void VecFunc (double* dst, double* arg, double (*func) (double), 
               unsigned long i1, unsigned long i2)
{
    for (unsigned long i=i1; i<i2; ++i)
        dst[i] = func (arg[i]);
}

/*
 * Vectorization of function applied to FieldVariable.
 */
void FieldFunc (FieldVariable_t* f, double (*func) (double))
{
    VecFunc (f->value, f->x, func, 0, f->N);
}

/* 
 * Some grid functions 
 */
double Gauss (double x) { return exp(-pow(x-0.5,2)/pow(0.05,2)); }
double Zero (double x) { return 0.0; }

int WriteToDisk (struct FieldVariable f, char* str)
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
void UpdateField (FieldVariable_t* dst, int dst1, int dst2, 
                  FieldVariable_t src,  int src1, double dt)
{
    for (int i=dst1; i<dst2; ++i, ++src1)
        dst->value[i] += dt*(src.value[src1+1] - src.value[src1])/src.dx;
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
    double L=1, cfl=1, T=0.2, c=1;

    /* Initialize */
    Field_t f;
    AllocField(&f.p, NX);
    AllocField(&f.u, NX+1);
    SetGrid (&f.p, 0.5*L/NX, L-0.5*L/NX);
    SetGrid (&f.u, 0, L);
    FieldFunc (&f.p, Gauss);
    FieldFunc (&f.u, Zero);

    /* Partition data */
    int cells_per_node = NX/NODES;
    int start, end;
    assert(cells_per_node*NODES == NX);

    /* Timestep */
    f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
    f.Nt = T/f.dt;

    for (int n=0; n<f.Nt; ++n) {

        /* Split the update to make multithreading easier */

        for (int i=0; i<NODES; i++) {
            start = i*cells_per_node;
            end = (i+1)*cells_per_node;
            UpdateField (&f.p, start, end, f.u, start, f.dt);
        }

        for (int i=0; i<NODES; i++) {
            start = 1 + i*cells_per_node;
            end = 1 + (i+1)*cells_per_node;
            if (i==NODES-1) /* last node gets one less */
                --end;
            UpdateField (&f.u, start, end, f.p, start-1, f.dt);
        }
        assert( (NODES-1 + 1)*cells_per_node == NX );
    }

    /* write data to disk */
    WriteToDisk(f.p, "output_p"); 
    WriteToDisk(f.u, "output_u"); 

    FreeField(f.p);
    FreeField(f.u);

    return EXIT_SUCCESS;
}
