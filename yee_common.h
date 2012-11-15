/*
 * Data structures and functions common to the various implementations.
 */
#ifndef YEE_COMMON_H
#define YEE_COMMON_H

/***************************************************************************
 * Data structures
 **************************************************************************/
typedef struct FieldVariable FieldVariable;
typedef struct Field Field;
typedef struct UpdateParam UpdateParam;
typedef struct FieldPartition FieldPartition;

struct FieldVariable {
    double* value;
    double* x;
    double dx;
    unsigned long size;    /* number of points */
};

struct Field {
    FieldVariable p;
    FieldVariable u;
    double dt;
    double Nt;
};

struct UpdateParam {
    FieldVariable* dst;
    int dst1;
    int dst2;
    FieldVariable* src;
    int src1;
    double dt;
};

struct FieldPartition {
    int start_p;
    int end_p;
    int start_u;    /* should be start_p+1 */
    int end_u;      /* should be start_u+1, except for last node */
};

/***************************************************************************
 * Functions
 **************************************************************************/

/* 
 * Allocate and deallocate the memory used by the grid (field) 
 */
void alloc_field (FieldVariable*, unsigned long size);
void free_field (FieldVariable);

/* 
 * Generate grid from start to end with N number of points  
 */
void set_grid (FieldVariable*, double start, double end);

/*
 * Vectorization of function.
 * Apply the function func on the array dst with argument arg.
 */
void vec_func (double* dst, double* arg, double (*func) (double),
               unsigned long i1, unsigned long i2);

/*
 * Vectorization of function applied to FieldVariable.
 */
void field_func (FieldVariable* f, double (*func) (double));

/*
 * Leapfrog time update
 * Update the field points dst for the indices dst1,...,dst2, using the
 * field points src using indices src1,...,src1 + (dst2-dst1)
 */
void update_field (FieldVariable* dst, int dst1, int dst2, 
                   FieldVariable* src, int src1, double dt);

void update_field2 (FieldVariable* dst, int dst1, int dst2, 
                    FieldVariable* src, int src1, double dt);

/*
 * Divide the grid for the different threads.
 */
FieldPartition partition_grid (int current_node, int nodes,
                               int cells_per_node);

/*
 * Divide the grid for the different threads.
 * This is the second version, suited for the MPI implementation. Here we
 * divide the grid according to cells.
 */
FieldPartition partition_grid_into_cells (int current_node, int nodes,
                                          int cells_per_node);

/* 
 * Collect parameters
 */
UpdateParam collect_param (FieldVariable* dst, int dst1, int dst2, 
                           FieldVariable* src, int src1, double dt);

/*
 * Parse commandline argument and set the number of grid points as well as
 * the number of threads.
 */
void parse_cmdline (unsigned long* nx, unsigned long* nodes, 
                    int argc, char* argv[]);

/*
 * Write field to dist
 */
int write_to_disk (FieldVariable f, char* str);

/* 
 * Some grid functions 
 */
double gauss (double);
double zero (double);

/*
 * Round up integer division
 */
int round_up_divide (int x, int y);

/*
 * Get time in a way that is compatible with threading
 */
double gettime (void);

#endif
