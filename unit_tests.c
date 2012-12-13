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

void test_alloc_field(void)
{
    struct field_variable fv;

    alloc_field(&fv, 8, 9);
    CU_ASSERT(fv.size_x == 8);
    CU_ASSERT(fv.size_y == 9);
    free_field(fv);
}

void test_set_grid(void)
{
    struct field_variable fv;
    const double x[] = { 5.0, 8.0 };
    const double y[] = { 3.0, 6.0 };
    const unsigned long cells_x = 8;
    const unsigned long cells_y = 9;
    const unsigned long nodes_x = cells_x + 1;
    const unsigned long nodes_y = cells_y + 1;
    double dx = (x[1] - x[0]) / (double) cells_x;
    double dy = (y[1] - y[0]) / (double) cells_y;

    alloc_field(&fv, nodes_x, nodes_y);
    set_grid(&fv, x, y);

    CU_ASSERT_DOUBLE_EQUAL(fv.dx, dx, TOL);
    CU_ASSERT_DOUBLE_EQUAL(fv.dy, dy, TOL);
    for (unsigned long i = 0; i < nodes_x; ++i)
        CU_ASSERT_DOUBLE_EQUAL(fv.x[i], x[0] + i * dx, TOL);
    for (unsigned long i = 0; i < nodes_y; ++i)
        CU_ASSERT_DOUBLE_EQUAL(fv.y[i], y[0] + i * dy, TOL);
    CU_ASSERT_DOUBLE_EQUAL(fv.x[nodes_x - 1], x[1], TOL);
    CU_ASSERT_DOUBLE_EQUAL(fv.y[nodes_y - 1], y[1], TOL);

    free_field(fv);
}

void test_vec_func(void)
{
    const int size = 8;
    double *dst = malloc(sizeof(double) * size);
    double *arg = malloc(sizeof(double) * size);
    for (int i = 0; i < size; ++i)
        arg[i] = i;

    vec_func(dst, zero, arg, size);
    for (int i = 0; i < size; ++i)
        CU_ASSERT_DOUBLE_EQUAL(dst[i], 0, TOL);

    vec_func(dst, identity, arg, size);
    for (int i = 0; i < size; ++i)
        CU_ASSERT_DOUBLE_EQUAL(dst[i], arg[i], TOL);

    free(dst);
    free(arg);
}

void test_vec_func2d(void)
{
    const unsigned long size_x = 8;
    const unsigned long size_y = 9;
    double *arg_x = malloc(sizeof(double) * size_x);
    double *arg_y = malloc(sizeof(double) * size_y);

    for (unsigned long k = 0; k < size_x; ++k)
        arg_x[k] = k;
    for (unsigned long k = 0; k < size_y; ++k)
        arg_y[k] = k;

    double *dst = malloc(sizeof(double) * size_x * size_y);
    double *arg = malloc(sizeof(double) * size_x * size_y);

    vec_func2d(dst, zero2d, arg_x, size_x, arg_y, size_y);
    for (unsigned long i = 0; i < size_x; ++i)
        for (unsigned long j = 0; j < size_y; ++j)
            CU_ASSERT_DOUBLE_EQUAL(dst[i + j * size_x], 0, TOL);

    vec_func2d(dst, identity2d, arg_x, size_x, arg_y, size_y);
    for (unsigned long i = 0; i < size_x; ++i)
        for (unsigned long j = 0; j < size_y; ++j)
            CU_ASSERT_DOUBLE_EQUAL(dst[i + j * size_x], arg_x[i], TOL);

    free(dst);
    free(arg);
}

void test_apply_func(void)
{
    struct field_variable fv;
    unsigned long size = 8;
    double x[] = { 0, 3 };
    double y[] = { 3, 4 };
    alloc_field(&fv, size, size);
    set_grid(&fv, x, y);

    apply_func(&fv, zero2d);
    for (unsigned long i = 0; i < size; ++i)
        for (unsigned long j = 0; j < size; ++j)
            CU_ASSERT_DOUBLE_EQUAL(fv.value[i + j * size], 0, TOL);

    apply_func(&fv, identity2d);
    for (unsigned long i = 0; i < size; ++i)
        for (unsigned long j = 0; j < size; ++j)
            CU_ASSERT_DOUBLE_EQUAL(fv.value[i + j * size], fv.x[i], TOL);

    free_field(fv);
}

void test_init_acoustic_field_internal(double x[2], double y[2],
                                       unsigned long cells)
{
    struct field f = init_acoustic_field(cells, cells, x, y);

    CU_ASSERT_DOUBLE_EQUAL(f.p.size_x, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.size_y, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.size_x, cells + 1.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.size_y, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.size_x, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.size_y, cells + 1.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.x[0], f.p.dx / 2.0 + x[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.y[0], f.p.dy / 2.0 + y[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.x[0], x[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.y[0], f.u.dy / 2.0 + y[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.x[0], f.v.dx / 2.0 + x[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.y[0], y[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dx, f.u.dx, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dx, f.v.dx, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dy, f.u.dy, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dy, f.v.dy, TOL);

    for (unsigned long i = 0; i < cells; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(f.p.x[i], f.u.x[i] + f.u.dx / 2.0, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.y[i], f.u.y[i], TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.x[i], f.v.x[i], TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.y[i], f.v.y[i] + f.v.dy / 2.0, TOL);
    }

    for (unsigned long i = 0; i < cells; ++i)
        for (unsigned long j = 0; j < cells; ++j)
            CU_ASSERT_DOUBLE_EQUAL(f.p.value[i + j * cells], 0.0, TOL);

    for (unsigned long i = 0; i < cells + 1; ++i)
        for (unsigned long j = 0; j < cells; ++j)
            CU_ASSERT_DOUBLE_EQUAL(f.u.value[i + j * (cells + 1)], 0.0,
                                   TOL);

    for (unsigned long i = 0; i < cells; ++i)
        for (unsigned long j = 0; j < cells + 1; ++j)
            CU_ASSERT_DOUBLE_EQUAL(f.v.value[i + j * cells], 0.0, TOL);

    free_acoustic_field(f);
}

void test_init_acoustic_field(void)
{
    double x[2], y[2];
    unsigned long cells;

    cells = 4;
    x[0] = 0;
    x[1] = 1;
    y[0] = 0;
    y[1] = 1;
    test_init_acoustic_field_internal(x, y, cells);

    cells = 43;
    x[0] = 4.0;
    x[1] = 9.0;
    y[0] = 3.0;
    y[1] = 7.0;
    test_init_acoustic_field_internal(x, y, cells);
}

void test_assign_and_get(void)
{
    double x[] = { 4, 7 };
    double y[] = { 3, 5 };
    unsigned long N = 8;
    struct field f = init_acoustic_field(N, N, x, y);
    unsigned long i, j;
    double val;

    /* used by macro indexing */
    double *p = f.p.value;
    double *u = f.u.value;
    double *v = f.v.value;
    unsigned long nx = N;

    for (unsigned n = 0; n < 1000; ++n) {
        i = rand() % N;
        j = rand() % N;
        val = (double) rand();
        assign_to(f.p, i, j, val);
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i + j * N], val, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i + j * N], get_from(f.p, i, j),
                               TOL);
        CU_ASSERT_DOUBLE_EQUAL(p(i, j), val, TOL);
    }

    for (unsigned n = 0; n < 1000; ++n) {
        i = rand() % (N + 1);
        j = rand() % N;
        val = (double) rand();
        assign_to(f.u, i, j, val);
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i + j * (N + 1)], val, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i + j * (N + 1)],
                               get_from(f.u, i, j), TOL);
        CU_ASSERT_DOUBLE_EQUAL(u(i, j), val, TOL);
    }

    for (unsigned n = 0; n < 1000; ++n) {
        i = rand() % N;
        j = rand() % (N + 1);
        val = (double) rand();
        assign_to(f.v, i, j, val);
        CU_ASSERT_DOUBLE_EQUAL(f.v.value[i + j * N], val, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.v.value[i + j * N], get_from(f.v, i, j),
                               TOL);
        CU_ASSERT_DOUBLE_EQUAL(v(i, j), val, TOL);
    }

    free_acoustic_field(f);
}

void test_leapfrog(void)
{
    double x[] = { 0, 1 };
    double y[] = { 0, 1 };
    unsigned long N = 2;        /* number of cells, not nodes */
    struct field f = init_acoustic_field(N, N, x, y);

    /* used by the indexing macro */
    double *p = f.p.value;
    double *u = f.u.value;
    double *v = f.v.value;
    unsigned long nx = f.p.size_x;

    /* test 1 */
    for (unsigned long i = 0; i < N; ++i)
        for (unsigned long j = 0; j < N; ++j)
            p(i, j) = 1;
    leapfrog(&f);
    for (unsigned long i = 0; i < N + 1; ++i)
        for (unsigned long j = 0; j < N; ++j)
            CU_ASSERT_DOUBLE_EQUAL(u(i, j), 0, TOL);
    for (unsigned long i = 0; i < N; ++i)
        for (unsigned long j = 0; j < N + 1; ++j)
            CU_ASSERT_DOUBLE_EQUAL(v(i, j), 0, TOL);

    /* test 2 */
    apply_func(&f.p, zero2d);
    for (unsigned long i = 0; i < N + 1; ++i)
        for (unsigned long j = 0; j < N; ++j)
            u(i, j) = 1;
    for (unsigned long i = 0; i < N; ++i)
        for (unsigned long j = 0; j < N + 1; ++j)
            v(i, j) = 1;
    leapfrog(&f);
    for (unsigned long i = 0; i < N; ++i)
        for (unsigned long j = 0; j < N; ++j)
            CU_ASSERT_DOUBLE_EQUAL(p(i, j), 0, TOL);

    free_acoustic_field(f);
}

/*[> struct partition partition_grid(int current_thread, int cells_per_thread);<]         */
/*void test_partition_grid(void)                                                          */
/*{                                                                                       */
/*    long cells = 8;                                                                     */
/*    struct partition part;                                                              */

/*    part = partition_grid(0, cells);                                                    */
/*    CU_ASSERT(part.p[0] == 0);                                                          */
/*    CU_ASSERT(part.p[1] == 7);                                                          */
/*    CU_ASSERT(part.u[0] == 1);                                                          */
/*    CU_ASSERT(part.u[1] == 7);                                                          */

/*    part = partition_grid(0, cells / 2);                                                */
/*    CU_ASSERT(part.p[0] == 0);                                                          */
/*    CU_ASSERT(part.p[1] == 3);                                                          */
/*    CU_ASSERT(part.u[0] == 1);                                                          */
/*    CU_ASSERT(part.u[1] == 3);                                                          */

/*    part = partition_grid(1, cells / 2);                                                */
/*    CU_ASSERT(part.p[0] == 4);                                                          */
/*    CU_ASSERT(part.p[1] == 7);                                                          */
/*    CU_ASSERT(part.u[0] == 4);                                                          */
/*    CU_ASSERT(part.u[1] == 7);                                                          */

/*    part = partition_grid(0, cells / 4);                                                */
/*    CU_ASSERT(part.p[0] == 0);                                                          */
/*    CU_ASSERT(part.p[1] == 1);                                                          */
/*    CU_ASSERT(part.u[0] == 1);                                                          */
/*    CU_ASSERT(part.u[1] == 1);                                                          */

/*    part = partition_grid(1, cells / 4);                                                */
/*    CU_ASSERT(part.p[0] == 2);                                                          */
/*    CU_ASSERT(part.p[1] == 3);                                                          */
/*    CU_ASSERT(part.u[0] == 2);                                                          */
/*    CU_ASSERT(part.u[1] == 3);                                                          */

/*    part = partition_grid(2, cells / 4);                                                */
/*    CU_ASSERT(part.p[0] == 4);                                                          */
/*    CU_ASSERT(part.p[1] == 5);                                                          */
/*    CU_ASSERT(part.u[0] == 4);                                                          */
/*    CU_ASSERT(part.u[1] == 5);                                                          */

/*    part = partition_grid(3, cells / 4);                                                */
/*    CU_ASSERT(part.p[0] == 6);                                                          */
/*    CU_ASSERT(part.p[1] == 7);                                                          */
/*    CU_ASSERT(part.u[0] == 6);                                                          */
/*    CU_ASSERT(part.u[1] == 7);                                                          */
/*}                                                                                       */

/*[> void expand_indices(struct partition partition, long *begin_p, long *end_p,<]        */
/*[>         long *size_p, long *begin_u, long *end_u, long *size_u);           <]        */
/*void test_expand_indices(void)                                                          */
/*{                                                                                       */
/*    long cells = 8;                                                                     */
/*    struct partition part;                                                              */
/*    long bp, ep, sp, bu, eu, su;                                                        */

/*    part = partition_grid(0, cells);                                                    */
/*    expand_indices(part, &bp, &ep, &sp, &bu, &eu, &su);                                 */
/*    CU_ASSERT(bp == 0);                                                                 */
/*    CU_ASSERT(ep == 7);                                                                 */
/*    CU_ASSERT(sp == 8);                                                                 */
/*    CU_ASSERT(bu == 1);                                                                 */
/*    CU_ASSERT(eu == 7);                                                                 */
/*    CU_ASSERT(su == 7);                                                                 */

/*    part = partition_grid(0, cells / 2);                                                */
/*    expand_indices(part, &bp, &ep, &sp, &bu, &eu, &su);                                 */
/*    CU_ASSERT(bp == 0);                                                                 */
/*    CU_ASSERT(ep == 3);                                                                 */
/*    CU_ASSERT(sp == 4);                                                                 */
/*    CU_ASSERT(bu == 1);                                                                 */
/*    CU_ASSERT(eu == 3);                                                                 */
/*    CU_ASSERT(su == 3);                                                                 */

/*    part = partition_grid(0, cells / 4);                                                */
/*    expand_indices(part, &bp, &ep, &sp, &bu, &eu, &su);                                 */
/*    CU_ASSERT(bp == 0);                                                                 */
/*    CU_ASSERT(ep == 1);                                                                 */
/*    CU_ASSERT(sp == 2);                                                                 */
/*    CU_ASSERT(bu == 1);                                                                 */
/*    CU_ASSERT(eu == 1);                                                                 */
/*    CU_ASSERT(su == 1);                                                                 */
/*}                                                                                       */

/*[> void verify_grid_integrity(struct partition partition, int tid, long nx, int<]       */
/*[>         numworkers, int left);                                              <]       */
/*void test_verify_grid(void)                                                             */
/*{                                                                                       */
/*    struct partition part;                                                              */
/*    part = partition_grid(0, 8);                                                        */
/*    verify_grid_integrity(part, 0, 8, 1, NONE);                                         */
/*    part = partition_grid(0, 4);                                                        */
/*    verify_grid_integrity(part, 0, 4, 1, NONE);                                         */
/*    part = partition_grid(1, 4);                                                        */
/*    verify_grid_integrity(part, 1, 4, 2, 0);                                            */
/*}                                                                                       */

/*[> void set_local_index(long size_p, long size_u, long left, long *local_begin_p,<]     */
/*[>         long *local_end_p, long *local_size_p, long *local_begin_u, long      <]     */
/*[>         *local_end_u, long *local_size_u);                                    <]     */
/*void test_set_local_index(void)                                                         */
/*{                                                                                       */
/*    long local_begin_p;                                                                 */
/*    long local_end_p;                                                                   */
/*    long local_size_p;                                                                  */
/*    long local_begin_u;                                                                 */
/*    long local_end_u;                                                                   */
/*    long local_size_u;                                                                  */

/*    long size_p = 4;                                                                    */
/*    long size_u = 3;                                                                    */
/*    long left = NONE;                                                                   */
/*    set_local_index(size_p, size_u, left, &local_begin_p, &local_end_p,                 */
/*                    &local_size_p, &local_begin_u, &local_end_u,                        */
/*                    &local_size_u);                                                     */
/*    CU_ASSERT(local_begin_p == 1);                                                      */
/*    CU_ASSERT(local_end_p == 4);                                                        */
/*    CU_ASSERT(local_size_p == 5);                                                       */
/*    CU_ASSERT(local_begin_u == 0);                                                      */
/*    CU_ASSERT(local_end_u == 2);                                                        */
/*    CU_ASSERT(local_size_u == 4);                                                       */
/*}                                                                                       */

/* void parse_cmdline(long *nx, long *threads, char *outfile_p, char *outfile_u,*/
/*         int argc, char *argv[]);                                             */
void test_parse_cmdline(void)
{
    long nx = 0;
    long threads = 0;
    char outfile[STR_SIZE];
    int argc = 5;
    char *argv[] = { "/usr/bin/yee", "-n", "8", "-t", "4" };
    parse_cmdline(&nx, &threads, outfile, argc, argv);
    CU_ASSERT(nx == 8);
    CU_ASSERT(threads == 4);
}

void test_zero(void)
{
    CU_ASSERT(zero(1) == 0.0);
    CU_ASSERT(zero(45.5) == 0.0);
    CU_ASSERT(zero(0) == 0.0);
    CU_ASSERT(zero(-5.0f) == 0.0);
}

void test_zero2d(void)
{
    CU_ASSERT(zero2d(1, 2) == 0.0);
    CU_ASSERT(zero2d(45.5, 5.3) == 0.0);
    CU_ASSERT(zero2d(0, 0) == 0.0);
    CU_ASSERT(zero2d(-5.0f, 5.0f) == 0.0);
}

void test_identity(void)
{
    CU_ASSERT(identity(1) == 1.0);
    CU_ASSERT(identity(45.5) == 45.5);
    CU_ASSERT(identity(0) == 0.0);
    CU_ASSERT(identity(-5.0f) == -5.0f);
}

void test_identity2d(void)
{
    CU_ASSERT(identity2d(1, 5) == 1.0);
    CU_ASSERT(identity2d(45.5, 4.6) == 45.5);
    CU_ASSERT(identity2d(0, 0) == 0.0);
    CU_ASSERT(identity2d(-5.0f, 7.0f) == -5.0f);
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
        || !CU_add_test(pSuite, "vec_func2d", test_vec_func2d)
        || !CU_add_test(pSuite, "apply_func", test_apply_func)
        || !CU_add_test(pSuite, "init_acoustic_field",
                        test_init_acoustic_field)
        || !CU_add_test(pSuite, "assign_to, get_from", test_assign_and_get)
        || !CU_add_test(pSuite, "leapfrog", test_leapfrog)
        /*|| !CU_add_test(pSuite, "partition_grid", test_partition_grid) */
        /*|| !CU_add_test(pSuite, "expand_indices", test_expand_indices) */
        /*|| !CU_add_test(pSuite, "verify_grid_integrity", test_verify_grid) */
        /*|| !CU_add_test(pSuite, "set_local_index", test_set_local_index) */
        || !CU_add_test(pSuite, "parse_cmdline", test_parse_cmdline)
        || !CU_add_test(pSuite, "zero", test_zero)
        || !CU_add_test(pSuite, "zero2d", test_zero2d)
        || !CU_add_test(pSuite, "identity", test_identity)
        || !CU_add_test(pSuite, "identity2d", test_identity)
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
