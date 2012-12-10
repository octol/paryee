/* 
 * Some very basic unit tests.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <CUnit/Basic.h>

#include "yee_common.h"

#define TOL 1e-14

int init_testsuite1(void)
{
    return 0;
}

int clean_testsuite1(void)
{
    return 0;
}

/* void alloc_field(struct field_variable *, long size);*/
void test_alloc_field(void)
{
    struct field_variable fv;

    alloc_field(&fv, 8);
    CU_ASSERT(fv.size == 8);
    free_field(fv);
}

/* void set_grid(struct field_variable *, double start, double end);*/
void test_set_grid(void)
{
    const double start = 5.0;
    const double end = 8.0;
    const long cells = 8;
    const long nodes = cells + 1;

    /* basic case */
    struct field_variable fv;
    double dx = (end - start) / (double) cells;

    alloc_field(&fv, nodes);
    set_grid(&fv, start, end);

    CU_ASSERT(fv.size == nodes);
    CU_ASSERT_DOUBLE_EQUAL(fv.dx, dx, TOL);
    for (int i = 0; i < nodes; ++i)
        CU_ASSERT_DOUBLE_EQUAL(fv.x[i], start + i * dx, TOL);
    CU_ASSERT_DOUBLE_EQUAL(fv.x[nodes - 1], end, TOL);

    free_field(fv);
}

/* void vec_func(double *dst, double *arg, double (*func) (double), long i1, long i2);*/
void test_vec_func(void)
{
    const int size = 8;
    double *dst = malloc(sizeof(double) * size);
    double *arg = malloc(sizeof(double) * size);
    for (int i = 0; i < size; ++i)
        arg[i] = i;

    vec_func(dst, arg, zero, 0, size);
    for (int i = 0; i < size; ++i)
        CU_ASSERT_DOUBLE_EQUAL(dst[i], 0, TOL);

    vec_func(dst, arg, identity, 0, size);
    for (int i = 0; i < size; ++i)
        CU_ASSERT_DOUBLE_EQUAL(dst[i], arg[i], TOL);

    free(dst);
    free(arg);
}

/* void apply_func(struct field_variable *f, double (*func) (double));*/
void test_apply_func(void)
{
    struct field_variable fv;
    long size = 8;
    double start = 0;
    double end = 1;
    alloc_field(&fv, size);
    set_grid(&fv, start, end);

    apply_func(&fv, zero);
    for (int i = 0; i < size; ++i)
        CU_ASSERT_DOUBLE_EQUAL(fv.value[i], 0, TOL);

    apply_func(&fv, identity);
    for (int i = 0; i < size; ++i)
        CU_ASSERT_DOUBLE_EQUAL(fv.value[i], fv.x[i], TOL);

    free_field(fv);
}

/* struct field* init_acoustic_field(struct field *, long cells,*/
/*                                   double start, double end); */
void test_init_acoustic_field(void)
{
    struct field f;
    long cells = 8;
    double start = 0;
    double end = 1;

    f = init_acoustic_field(cells, start, end);

    CU_ASSERT_DOUBLE_EQUAL(f.p.size, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.size, cells + 1, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dx, f.u.dx, TOL);

    for (int i = 0; i < cells; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(f.p.x[i], f.u.x[i] + f.u.dx / 2, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i], 0, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i], 0, TOL);
    }
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[cells + 1], 0, TOL);

    free_acoustic_field(f);
}

/* void update_field_s(struct field_variable *restrict dst, int idst, int size,*/
/*         struct field_variable *restrict src, int isrc, double dt);          */
void test_update_field_s(void)
{
    double start = 0;
    double end = 1;
    long N = 8;                 /* number of cells, not nodes */
    struct field f = init_acoustic_field(N, start, end);

    CU_ASSERT_DOUBLE_EQUAL(f.p.dx, f.u.dx, TOL);

    /* 
     * Update p from u 
     */
    for (int i = 0; i < f.u.size; ++i)
        f.u.value[i] = i;
    /* simple test were we only update one node point */
    update_field_s(&f.p, 0, 1, &f.u, 0, 1.0f);
    CU_ASSERT_DOUBLE_EQUAL(f.p.value[0], 1 / f.u.dx, TOL);
    for (int i = 1; i < N; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i], 0, TOL);
    /* now update all node points */
    apply_func(&f.p, zero);
    update_field_s(&f.p, 0, f.p.size, &f.u, 0, 1.0f);
    for (int i = 0; i < N; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i], 1 / f.u.dx, TOL);
    /* now update all node points (with dt != 1) */
    apply_func(&f.p, zero);
    update_field_s(&f.p, 0, f.p.size, &f.u, 0, 42.2);
    for (int i = 0; i < N; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i], 1 * 42.2 / f.u.dx, TOL);

    /* 
     * Update u from p
     */
    apply_func(&f.u, zero);
    for (int i = 0; i < f.p.size; ++i)
        f.p.value[i] = i;

    /* simple test were we only update one node point */
    /* note that we leave the outer boundary point alone */
    update_field_s(&f.u, 1, 1, &f.p, 0, 1.0f);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[N], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[1], 1 / f.p.dx, TOL);
    for (int i = 2; i < N - 1; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i], 0, TOL);
    /* now update all (inner) node points */
    apply_func(&f.u, zero);
    update_field_s(&f.u, 1, f.u.size - 2, &f.p, 0, 1.0f);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[N], 0, TOL);
    for (int i = 1; i < N - 1; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i], 1 / f.p.dx, TOL);
    /* now update all node points (with dt != 1) */
    apply_func(&f.u, zero);
    update_field_s(&f.u, 1, f.u.size - 2, &f.p, 0, 42.2);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[N], 0, TOL);
    for (int i = 1; i < N - 1; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i], 1 * 42.2 / f.p.dx, TOL);

    free_acoustic_field(f);
}

/* void update_field_i(struct field_variable *restrict dst, int dst1, int dst2,*/
/*         struct field_variable *restrict src, int src1, double dt);          */
void test_update_field_i(void)
{
    double start = 0;
    double end = 1;
    long N = 8;                 /* number of cells, not nodes */
    struct field f = init_acoustic_field(N, start, end);

    /* 
     * Update p from u 
     */
    for (int i = 0; i < f.u.size; ++i)
        f.u.value[i] = i;
    /* simple test were we only update one node point */
    update_field_i(&f.p, 0, 1, &f.u, 0, 1.0f);
    CU_ASSERT_DOUBLE_EQUAL(f.p.value[0], 1 / f.u.dx, TOL);
    for (int i = 1; i < N; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i], 0, TOL);
    /* now update all node points */
    apply_func(&f.p, zero);
    update_field_i(&f.p, 0, f.p.size, &f.u, 0, 1.0f);
    for (int i = 0; i < N; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i], 1 / f.u.dx, TOL);
    /* now update all node points (with dt != 1) */
    apply_func(&f.p, zero);
    update_field_i(&f.p, 0, f.p.size, &f.u, 0, 42.2);
    for (int i = 0; i < N; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i], 1 * 42.2 / f.u.dx, TOL);

    /* 
     * Update u from p
     */
    apply_func(&f.u, zero);
    for (int i = 0; i < f.p.size; ++i)
        f.p.value[i] = i;

    /* simple test were we only update one node point */
    /* note that we leave the outer boundary point alone */
    update_field_i(&f.u, 1, 2, &f.p, 0, 1.0f);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[N], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[1], 1 / f.p.dx, TOL);
    for (int i = 2; i < N - 1; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i], 0, TOL);
    /* now update all (inner) node points */
    apply_func(&f.u, zero);
    update_field_i(&f.u, 1, f.u.size - 1, &f.p, 0, 1.0f);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[N], 0, TOL);
    for (int i = 1; i < N; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i], 1 / f.p.dx, TOL);
    /* now update all node points (with dt != 1) */
    apply_func(&f.u, zero);
    update_field_i(&f.u, 1, f.u.size - 1, &f.p, 0, 42.2);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[N], 0, TOL);
    for (int i = 1; i < N; ++i)
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i], 1 * 42.2 / f.p.dx, TOL);

    free_acoustic_field(f);
}

/* struct partition partition_grid(int current_thread, int cells_per_thread);*/
void test_partition_grid(void)
{
    long cells = 8;
    struct partition part;

    part = partition_grid(0, cells);
    CU_ASSERT(part.p[0] == 0);
    CU_ASSERT(part.p[1] == 7);
    CU_ASSERT(part.u[0] == 1);
    CU_ASSERT(part.u[1] == 7);

    part = partition_grid(0, cells / 2);
    CU_ASSERT(part.p[0] == 0);
    CU_ASSERT(part.p[1] == 3);
    CU_ASSERT(part.u[0] == 1);
    CU_ASSERT(part.u[1] == 3);

    part = partition_grid(1, cells / 2);
    CU_ASSERT(part.p[0] == 4);
    CU_ASSERT(part.p[1] == 7);
    CU_ASSERT(part.u[0] == 4);
    CU_ASSERT(part.u[1] == 7);

    part = partition_grid(0, cells / 4);
    CU_ASSERT(part.p[0] == 0);
    CU_ASSERT(part.p[1] == 1);
    CU_ASSERT(part.u[0] == 1);
    CU_ASSERT(part.u[1] == 1);

    part = partition_grid(1, cells / 4);
    CU_ASSERT(part.p[0] == 2);
    CU_ASSERT(part.p[1] == 3);
    CU_ASSERT(part.u[0] == 2);
    CU_ASSERT(part.u[1] == 3);

    part = partition_grid(2, cells / 4);
    CU_ASSERT(part.p[0] == 4);
    CU_ASSERT(part.p[1] == 5);
    CU_ASSERT(part.u[0] == 4);
    CU_ASSERT(part.u[1] == 5);

    part = partition_grid(3, cells / 4);
    CU_ASSERT(part.p[0] == 6);
    CU_ASSERT(part.p[1] == 7);
    CU_ASSERT(part.u[0] == 6);
    CU_ASSERT(part.u[1] == 7);
}

/* void expand_indices(struct partition partition, long *begin_p, long *end_p,*/
/*         long *size_p, long *begin_u, long *end_u, long *size_u);           */
void test_expand_indices(void)
{
    long cells = 8;
    struct partition part;
    long bp, ep, sp, bu, eu, su;

    part = partition_grid(0, cells);
    expand_indices(part, &bp, &ep, &sp, &bu, &eu, &su);
    CU_ASSERT(bp == 0);
    CU_ASSERT(ep == 7);
    CU_ASSERT(sp == 8);
    CU_ASSERT(bu == 1);
    CU_ASSERT(eu == 7);
    CU_ASSERT(su == 7);

    part = partition_grid(0, cells / 2);
    expand_indices(part, &bp, &ep, &sp, &bu, &eu, &su);
    CU_ASSERT(bp == 0);
    CU_ASSERT(ep == 3);
    CU_ASSERT(sp == 4);
    CU_ASSERT(bu == 1);
    CU_ASSERT(eu == 3);
    CU_ASSERT(su == 3);

    part = partition_grid(0, cells / 4);
    expand_indices(part, &bp, &ep, &sp, &bu, &eu, &su);
    CU_ASSERT(bp == 0);
    CU_ASSERT(ep == 1);
    CU_ASSERT(sp == 2);
    CU_ASSERT(bu == 1);
    CU_ASSERT(eu == 1);
    CU_ASSERT(su == 1);
}

/* void verify_grid_integrity(struct partition partition, int tid, long nx, int*/
/*         numworkers, int left);                                              */
void test_verify_grid(void)
{
    struct partition part;
    part = partition_grid(0, 8);
    verify_grid_integrity(part, 0, 8, 1, NONE);
    part = partition_grid(0, 4);
    verify_grid_integrity(part, 0, 4, 1, NONE);
    part = partition_grid(1, 4);
    verify_grid_integrity(part, 1, 4, 2, 0);
}

/* void set_local_index(long size_p, long size_u, long left, long *local_begin_p,*/
/*         long *local_end_p, long *local_size_p, long *local_begin_u, long      */
/*         *local_end_u, long *local_size_u);                                    */
void test_set_local_index(void)
{
    long local_begin_p;
    long local_end_p;
    long local_size_p;
    long local_begin_u;
    long local_end_u;
    long local_size_u;

    long size_p = 4;
    long size_u = 3;
    long left = NONE;
    set_local_index(size_p, size_u, left, &local_begin_p, &local_end_p,
                    &local_size_p, &local_begin_u, &local_end_u,
                    &local_size_u);
    CU_ASSERT(local_begin_p == 1);
    CU_ASSERT(local_end_p == 4);
    CU_ASSERT(local_size_p == 5);
    CU_ASSERT(local_begin_u == 0);
    CU_ASSERT(local_end_u == 2);
    CU_ASSERT(local_size_u == 4);
}

/* void parse_cmdline(long *nx, long *threads, char *outfile_p, char *outfile_u,*/
/*         int argc, char *argv[]);                                             */
void test_parse_cmdline(void)
{
    long nx = 0;
    long threads = 0;
    char outfile_p[STR_SIZE];
    char outfile_u[STR_SIZE];
    int argc = 3;
    char *argv[] = { "/usr/bin/yee", "-n", "8" };
    parse_cmdline(&nx, &threads, outfile_p, outfile_u, argc, argv);
    CU_ASSERT(nx == 8);
}

/* double zero(double);*/
void test_zero(void)
{
    CU_ASSERT(zero(1) == 0.0);
    CU_ASSERT(zero(45.5) == 0.0);
    CU_ASSERT(zero(0) == 0.0);
    CU_ASSERT(zero(-5.0f) == 0.0);
}

/* double identity(double);*/
void test_identity(void)
{
    CU_ASSERT(identity(1) == 1.0);
    CU_ASSERT(identity(45.5) == 45.5);
    CU_ASSERT(identity(0) == 0.0);
    CU_ASSERT(identity(-5.0f) == -5.0f);
}

/* int round_up_divide(int x, int y);*/
void test_round_up_divide(void)
{
    CU_ASSERT(round_up_divide(1, 3) == 1);
    CU_ASSERT(round_up_divide(2, 3) == 1);
    CU_ASSERT(round_up_divide(3, 3) == 1);
    CU_ASSERT(round_up_divide(4, 3) == 2);
    CU_ASSERT(round_up_divide(0, 3) == 1);
    CU_ASSERT(round_up_divide(7, 3) == 3);
}

int main()
{
    CU_pSuite pSuite = NULL;

    /* Init CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add suite to the registry */
    pSuite = CU_add_suite("yee_common", init_testsuite1, clean_testsuite1);
    if (!pSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Add the tests to the suite */
    if (!CU_add_test(pSuite, "alloc_field", test_alloc_field)
        || !CU_add_test(pSuite, "set_grid", test_set_grid)
        || !CU_add_test(pSuite, "vec_func", test_vec_func)
        || !CU_add_test(pSuite, "apply_func", test_apply_func)
        || !CU_add_test(pSuite, "init_acoustic_field",
                        test_init_acoustic_field)
        || !CU_add_test(pSuite, "update_field_s", test_update_field_s)
        || !CU_add_test(pSuite, "update_field_i", test_update_field_i)
        || !CU_add_test(pSuite, "partition_grid", test_partition_grid)
        || !CU_add_test(pSuite, "expand_indices", test_expand_indices)
        /*|| !CU_add_test(pSuite, "verify_grid_integrity", test_verify_grid) */
        /*|| !CU_add_test(pSuite, "set_local_index", test_set_local_index) */
        || !CU_add_test(pSuite, "parse_cmdline", test_parse_cmdline)
        || !CU_add_test(pSuite, "zero", test_zero)
        || !CU_add_test(pSuite, "identity", test_identity)
        || !CU_add_test(pSuite, "round_up_divide", test_round_up_divide)) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Run all the tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}
