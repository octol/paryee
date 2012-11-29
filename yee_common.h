/*
 * Data structures and functions common to the various implementations.
 */
#ifndef YEE_COMMON_H
#define YEE_COMMON_H

#define MASTER 0                /* id of master process */
#define MINWORKER 1
#define MAXWORKER 8
#define NONE 0                  /* no neighbour */
#define STR_SIZE 256

/***************************************************************************
 * Data structures
 **************************************************************************/
struct FieldVariable {
    double *value;
    double *x;
    double dx;
    unsigned long size;         /* number of points */
};

struct Field {
    struct FieldVariable p;
    struct FieldVariable u;
    double dt;
    double Nt;
};

struct UpdateParam {
    struct FieldVariable *dst;
    int dst1;
    int dst2;
    struct FieldVariable *src;
    int src1;
    double dt;
};

struct Partition {
    int p[2];                   /* start end indices of internal domain */
    int u[2];
};

/***************************************************************************
 * Functions
 **************************************************************************/

/* 
 * Allocate and deallocate the memory used by the grid (field) 
 */
void alloc_field(struct FieldVariable *, unsigned long size);
void free_field(struct FieldVariable);

/* 
 * Generate grid from start to end with N number of points  
 */
void set_grid(struct FieldVariable *, double start, double end);

/*
 * Vectorization of function.
 * Apply the function func on the array dst with argument arg.
 */
void vec_func(double *dst, double *arg, double (*func) (double),
              unsigned long i1, unsigned long i2);

/*
 * Vectorization of function applied to FieldVariable.
 */
void apply_func(struct FieldVariable *f, double (*func) (double));

/* 
 * Leapfrog time update.
 * Update the size number of field points starting at the index idst, using
 * the field points starting at isrc. Note: _s in the name indicates the
 * input (size).
 */
void update_field_s(struct FieldVariable *dst, int idst, int size,
                    struct FieldVariable *src, int isrc, double dt);

/* 
 * Leapfrog time update.
 * Update the field points starting at the index dst1 and ending at dst2,
 * using the field points starting at src1. Note: _i in the name indicates
 * the input (index).
 */
void update_field_i(struct FieldVariable *dst, int dst1, int dst2,
                    struct FieldVariable *src, int src1, double dt);

/*
 * Divide the grid for the different threads.
 * We divide the grid according to cells.
 * Note: only the inner nodes are returned.
 */
struct Partition partition_grid(int current_node, int cells_per_node);

/* 
 * Expand the partition struct into indices 
 */
void expand_indices(struct Partition partition,
                    int *begin_p, int *end_p, int *size_p,
                    int *begin_u, int *end_u, int *size_u);

/*
 * Checks that the grid is according to the specs.
 */
void verify_grid_integrity(struct Partition partition, int tid,
                           unsigned long nx, int numworkers, int left);

/*
 * Compute the local start/end indices and local array sizes from the
 * partition sizes as well as if the we are on the left boundary.
 */
void set_local_index(int size_p, int size_u, int left,
                     int *local_begin_p, int *local_end_p,
                     int *local_size_p, int *local_begin_u,
                     int *local_end_u, int *local_size_u);

/* 
 * Collect parameters
 */
struct UpdateParam collect_param(struct FieldVariable *dst, int dst1,
                                 int dst2, struct FieldVariable *src,
                                 int src1, double dt);

/*
 * Parse commandline argument and set the number of grid points, * the
 * number of threads as well as the output filenames.
 */
void parse_cmdline(unsigned long *nx, unsigned long *nodes,
                   char *outfile_p, char *outfile_u,
                   int argc, char *argv[]);

/*
 * Write field to dist
 */
int write_to_disk(struct FieldVariable f, char *str);

/* 
 * Some grid functions 
 */
double gauss(double);
double zero(double);

/*
 * Round up integer division
 */
int round_up_divide(int x, int y);

/*
 * Get time in a way that is compatible with threading
 */
double gettime(void);

#endif
