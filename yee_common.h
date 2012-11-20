/*
 * Data structures and functions common to the various implementations.
 */
#ifndef YEE_COMMON_H
#define YEE_COMMON_H

#define field_func apply_func

/***************************************************************************
 * Data structures
 **************************************************************************/
typedef struct FieldVariable FieldVariable;
struct FieldVariable {
    double* value;
    double* x;
    double dx;
    unsigned long size;    /* number of points */
};

typedef struct Field Field;
struct Field {
    FieldVariable p;
    FieldVariable u;
    double dt;
    double Nt;
};

typedef struct UpdateParam UpdateParam;
struct UpdateParam {
    FieldVariable* dst;
    int dst1;
    int dst2;
    FieldVariable* src;
    int src1;
    double dt;
};

typedef struct FieldPartition FieldPartition;
struct FieldPartition {
    int start_p;
    int end_p;
    int start_u;   
    int end_u;    
};

typedef struct Partition Partition;
struct Partition {
    /* The grid points that is selected to be updated */
    /*FieldPartition inner;*/
    /* The inner grid points plus the extra grid points needed to update the
     * inner grid points */
    /*FieldPartition outer;*/
    int p[2];   /* start and end indices of internal domain */
    int u[2];
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

/*
 * Leapfrog time update
 * Update the size number of field points starting at the index idst, using
 * the field points starting at isrc.
 */
void update_field2 (FieldVariable* dst, int idst, int size, 
                    FieldVariable* src, int isrc, double dt);

void update_field3 (FieldVariable* dst, int dst1, int dst2, 
                    FieldVariable* src, int src1, double dt);

/*
 * Divide the grid for the different threads.
 */
FieldPartition 
partition_grid (int current_node, int nodes, int cells_per_node);

/*
 * Divide the grid for the different threads.
 * This is the second version, suited for the MPI implementation. Here we
 * divide the grid according to cells.
 */
Partition partition_grid2 (int current_node, int nodes, int cells_per_node);

/* 
 * Expand the partition struct into indices 
 */
void expand_indices (Partition partition, 
                     int* begin_p, int* end_p, int* size_p, 
                     int* begin_u, int* end_u, int* size_u);

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
