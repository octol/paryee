/* 
 * Some very basic unit tests.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <CUnit/Basic.h>

#include "yee_common.h"

#define TOL 1e-14

int init_suite1(void)
{
    return 0;
}

int clean_suite1(void)
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
    const double start = 0.0;
    const double end = 1.0;
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

    init_acoustic_field(&f, cells, start, end);

    CU_ASSERT_DOUBLE_EQUAL(f.p.dx, f.u.dx, TOL);
    for (int i = 0; i < cells; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(f.p.x[i], f.u.x[i] + f.u.dx / 2, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i], 0, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i], 0, TOL);
    }
    CU_ASSERT_DOUBLE_EQUAL(f.u.value[cells + 1], 0, TOL);
}

/* void update_field_s(struct field_variable *restrict dst, int idst, int size,*/
/*         struct field_variable *restrict src, int isrc, double dt);          */
void test_update_field_s(void)
{
    struct field_variable p;
    struct field_variable u;
    double start = 0;
    double end = 1;
    long N = 8;                 /* number of cells, not nodes */

    alloc_field(&p, N);
    alloc_field(&u, N + 1);
    set_grid(&u, start, end);
    set_grid(&p, start + u.dx / 2, end - u.dx / 2);
    apply_func(&p, zero);
    apply_func(&u, zero);

    CU_ASSERT_DOUBLE_EQUAL(p.dx, u.dx, TOL);

    /* first a simple test were we only update one node point */
    for (int i = 0; i < N; ++i)
        u.value[i] = i;
    update_field_s(&p, 0, 1, &u, 0, 1.0f);

    CU_ASSERT_DOUBLE_EQUAL(p.value[0], 1 / u.dx, TOL);
    for (int i = 1; i < N; ++i)
        CU_ASSERT_DOUBLE_EQUAL(p.value[i], 0, TOL);

    free_field(p);
    free_field(u);
}

/* void update_field_i(struct field_variable *restrict dst, int dst1, int dst2,*/
/*         struct field_variable *restrict src, int src1, double dt);          */
void test_update_field_i(void)
{
    CU_FAIL('Not implemented yet');
}

/* struct partition partition_grid(int current_thread, int cells_per_thread);*/
void test_partition_grid(void)
{
    CU_FAIL('Not implemented yet');
}

/* void expand_indices(struct partition partition, long *begin_p, long *end_p,*/
/*         long *size_p, long *begin_u, long *end_u, long *size_u);           */
void test_expand_indices(void)
{
    CU_FAIL('Not implemented yet');
}

/* void verify_grid_integrity(struct partition partition, int tid, long nx, int*/
/*         numworkers, int left);                                              */
void test_verify_grid(void)
{
    CU_FAIL('Not implemented yet');
}

/* void set_local_index(long size_p, long size_u, long left, long *local_begin_p,*/
/*         long *local_end_p, long *local_size_p, long *local_begin_u, long      */
/*         *local_end_u, long *local_size_u);                                    */
void test_set_local_index(void)
{
    CU_FAIL('Not implemented yet');
}

/* void parse_cmdline(long *nx, long *threads, char *outfile_p, char *outfile_u,*/
/*         int argc, char *argv[]);                                             */
void test_parse_cmdline(void)
{
    CU_FAIL('Not implemented yet');
}

/* int write_to_disk(struct field_variable f, char *str);*/
void test_write_to_disk(void)
{
    CU_FAIL('Not implemented yet');
}

/* double gauss(double);*/
void test_gauss(void)
{
    CU_FAIL('Not implemented yet');
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

/* double gettime(void);*/
void test_gettime(void)
{
    CU_FAIL('Not implemented yet');
}

int main()
{
    CU_pSuite pSuite = NULL;

    /* Init CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add suite to the registry */
    pSuite = CU_add_suite("yee_common", init_suite1, clean_suite1);
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
        || !CU_add_test(pSuite, "verify_grid_integrity", test_verify_grid)
        || !CU_add_test(pSuite, "set_local_index", test_set_local_index)
        || !CU_add_test(pSuite, "parse_cmdline", test_parse_cmdline)
        || !CU_add_test(pSuite, "write_to_disk", test_write_to_disk)
        || !CU_add_test(pSuite, "gauss", test_gauss)
        || !CU_add_test(pSuite, "zero", test_zero)
        || !CU_add_test(pSuite, "identity", test_identity)
        || !CU_add_test(pSuite, "round_up_divide", test_round_up_divide)
        || !CU_add_test(pSuite, "gettime", test_gettime)) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Run all the tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}
