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
    struct py_field_variable fv;

    py_alloc_field(&fv, 8, 9);
    CU_ASSERT(fv.size_x == 8);
    CU_ASSERT(fv.size_y == 9);
    py_free_field(fv);
}

void test_set_grid(void)
{
    struct py_field_variable fv;
    const double x[] = { 5.0, 8.0 };
    const double y[] = { 3.0, 6.0 };
    const unsigned long cells_x = 8;
    const unsigned long cells_y = 9;
    const unsigned long nodes_x = cells_x + 1;
    const unsigned long nodes_y = cells_y + 1;
    double dx = (x[1] - x[0]) / (double) cells_x;
    double dy = (y[1] - y[0]) / (double) cells_y;

    py_alloc_field(&fv, nodes_x, nodes_y);
    py_set_grid(&fv, x, y);

    CU_ASSERT_DOUBLE_EQUAL(fv.dx, dx, TOL);
    CU_ASSERT_DOUBLE_EQUAL(fv.dy, dy, TOL);
    for (unsigned long i = 0; i < nodes_x; ++i)
        CU_ASSERT_DOUBLE_EQUAL(fv.x[i], x[0] + i * dx, TOL);
    for (unsigned long i = 0; i < nodes_y; ++i)
        CU_ASSERT_DOUBLE_EQUAL(fv.y[i], y[0] + i * dy, TOL);
    CU_ASSERT_DOUBLE_EQUAL(fv.x[nodes_x - 1], x[1], TOL);
    CU_ASSERT_DOUBLE_EQUAL(fv.y[nodes_y - 1], y[1], TOL);

    py_free_field(fv);
}

void test_vec_func(void)
{
    const int size = 8;
    double *dst = malloc(sizeof(double) * size);
    double *arg = malloc(sizeof(double) * size);
    for (int i = 0; i < size; ++i)
        arg[i] = i;

    py_vec_func(dst, zero, arg, size);
    for (int i = 0; i < size; ++i)
        CU_ASSERT_DOUBLE_EQUAL(dst[i], 0, TOL);

    py_vec_func(dst, identity, arg, size);
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

    py_vec_func2d(dst, zero2d, arg_x, size_x, arg_y, size_y);
    for (unsigned long i = 0; i < size_x; ++i)
        for (unsigned long j = 0; j < size_y; ++j)
            CU_ASSERT_DOUBLE_EQUAL(dst[i + j * size_x], 0, TOL);

    py_vec_func2d(dst, identity2d, arg_x, size_x, arg_y, size_y);
    for (unsigned long i = 0; i < size_x; ++i)
        for (unsigned long j = 0; j < size_y; ++j)
            CU_ASSERT_DOUBLE_EQUAL(dst[i + j * size_x], arg_x[i], TOL);

    free(dst);
    free(arg);
}

void test_apply_func(void)
{
    struct py_field_variable fv;
    unsigned long size = 8;
    double x[] = { 0, 3 };
    double y[] = { 3, 4 };
    py_alloc_field(&fv, size, size);
    py_set_grid(&fv, x, y);

    py_apply_func(&fv, zero2d);
    for (unsigned long i = 0; i < size; ++i)
        for (unsigned long j = 0; j < size; ++j)
            CU_ASSERT_DOUBLE_EQUAL(fv.value[i + j * size], 0, TOL);

    py_apply_func(&fv, identity2d);
    for (unsigned long i = 0; i < size; ++i)
        for (unsigned long j = 0; j < size; ++j)
            CU_ASSERT_DOUBLE_EQUAL(fv.value[i + j * size], fv.x[i], TOL);

    py_free_field(fv);
}

void test_init_acoustic_field_internal(double x[2], double y[2],
                                       unsigned long cells)
{
    struct py_field f = py_init_acoustic_field(cells, cells, x, y);

    CU_ASSERT_DOUBLE_EQUAL(f.p.size_x, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.size_y, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.size_x, cells + 1.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.size_y, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.size_x, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.size_y, cells + 1.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dx, f.u.dx, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dx, f.v.dx, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dy, f.u.dy, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dy, f.v.dy, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.x[0], f.p.dx / 2.0 + x[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.y[0], f.p.dy / 2.0 + y[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.x[0], x[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.y[0], f.u.dy / 2.0 + y[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.x[0], f.v.dx / 2.0 + x[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.y[0], y[0], TOL);

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

    py_free_acoustic_field(f);
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

void test_init_local_acoustic_field_internal(double x[2], double y[2],
                                             unsigned long cells)
{
    struct py_field f = py_init_local_acoustic_field(cells, cells, x, y);

    CU_ASSERT_DOUBLE_EQUAL(f.p.size_x, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.size_y, cells + 1.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.size_x, cells + 1.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.size_y, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.size_x, cells, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.size_y, cells + 1.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dx, f.u.dx, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dx, f.v.dx, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dy, f.u.dy, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.dy, f.v.dy, TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.x[0], f.p.dx / 2.0 + x[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.p.y[0], -f.p.dy / 2.0 + y[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.x[0], x[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.u.y[0], f.u.dy / 2.0 + y[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.x[0], f.v.dx / 2.0 + x[0], TOL);
    CU_ASSERT_DOUBLE_EQUAL(f.v.y[0], y[0], TOL);


    CU_ASSERT_DOUBLE_EQUAL(f.p.y[0], f.u.y[0] - f.p.dy, TOL);
    for (unsigned long i = 0; i < cells; ++i) {
        CU_ASSERT_DOUBLE_EQUAL(f.p.x[i], f.u.x[i] + f.u.dx / 2.0, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.y[i + 1], f.u.y[i], TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.x[i], f.v.x[i], TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.y[i], f.v.y[i] - f.v.dy / 2.0, TOL);
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

    py_free_acoustic_field(f);
}

void test_init_local_acoustic_field(void)
{
    double x[2], y[2];
    unsigned long cells;

    cells = 2;
    x[0] = 0;
    x[1] = 1;
    y[0] = 0;
    y[1] = 1;
    test_init_local_acoustic_field_internal(x, y, cells);

    cells = 43;
    x[0] = 4.0;
    x[1] = 9.0;
    y[0] = 3.0;
    y[1] = 7.0;
    test_init_local_acoustic_field_internal(x, y, cells);
}


void test_assign_and_get(void)
{
    double x[] = { 4, 7 };
    double y[] = { 3, 5 };
    unsigned long N = 8;
    struct py_field f = py_init_acoustic_field(N, N, x, y);
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
        py_assign_to(f.p, i, j, val);
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i + j * N], val, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.p.value[i + j * N], py_get_from(f.p, i, j),
                               TOL);
        CU_ASSERT_DOUBLE_EQUAL(P(i, j), val, TOL);
    }

    for (unsigned n = 0; n < 1000; ++n) {
        i = rand() % (N + 1);
        j = rand() % N;
        val = (double) rand();
        py_assign_to(f.u, i, j, val);
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i + j * (N + 1)], val, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.u.value[i + j * (N + 1)],
                               py_get_from(f.u, i, j), TOL);
        CU_ASSERT_DOUBLE_EQUAL(U(i, j), val, TOL);
    }

    for (unsigned n = 0; n < 1000; ++n) {
        i = rand() % N;
        j = rand() % (N + 1);
        val = (double) rand();
        py_assign_to(f.v, i, j, val);
        CU_ASSERT_DOUBLE_EQUAL(f.v.value[i + j * N], val, TOL);
        CU_ASSERT_DOUBLE_EQUAL(f.v.value[i + j * N], py_get_from(f.v, i, j),
                               TOL);
        CU_ASSERT_DOUBLE_EQUAL(V(i, j), val, TOL);
    }

    py_free_acoustic_field(f);
}

void test_set_boundary(void)
{
    double x[] = { 0, 1 };
    double y[] = { 0, 1 };
    struct py_field f = py_init_acoustic_field(2, 4, x, y);

    /* used by the indexing macro */
    double *u = f.u.value;
    double *v = f.v.value;
    unsigned long nx = f.p.size_x;

    /*
     *   v   v
     * u p u p u
     *   v   v
     * u p u p u
     *   v   v
     * u p u p u
     *   v   v
     * u p u p u
     *   v   v
     */
    CU_ASSERT_DOUBLE_EQUAL(V(0, 0), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(1, 0), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(0, 1), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(1, 1), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(0, 2), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(1, 2), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(0, 3), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(1, 3), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(0, 4), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(1, 4), 0, TOL);

    CU_ASSERT_DOUBLE_EQUAL(U(0, 0), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(0, 1), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(0, 2), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(0, 3), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(1, 0), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(1, 1), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(1, 2), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(1, 3), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(2, 0), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(2, 1), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(2, 2), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(2, 3), 0, TOL);

    py_apply_func(&f.u, one2d);    /* initial data */
    py_apply_func(&f.v, one2d);    /* initial data */

    CU_ASSERT_DOUBLE_NOT_EQUAL(V(0, 0), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(1, 0), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(0, 1), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(1, 1), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(0, 2), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(1, 2), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(0, 3), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(1, 3), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(0, 4), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(1, 4), 0, TOL);

    CU_ASSERT_DOUBLE_NOT_EQUAL(U(0, 0), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(0, 1), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(0, 2), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(0, 3), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(1, 0), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(1, 1), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(1, 2), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(1, 3), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(2, 0), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(2, 1), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(2, 2), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(2, 3), 0, TOL);

    py_set_boundary(&f);

    CU_ASSERT_DOUBLE_EQUAL(V(0, 0), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(1, 0), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(0, 1), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(1, 1), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(0, 2), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(1, 2), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(0, 3), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(V(1, 3), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(0, 4), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(V(1, 4), 0, TOL);

    CU_ASSERT_DOUBLE_EQUAL(U(0, 0), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(0, 1), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(0, 2), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(0, 3), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(1, 0), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(1, 1), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(1, 2), 0, TOL);
    CU_ASSERT_DOUBLE_NOT_EQUAL(U(1, 3), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(2, 0), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(2, 1), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(2, 2), 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(U(2, 3), 0, TOL);
}

void test_leapfrog(void)
{
    double x[] = { 0, 1 };
    double y[] = { 0, 1 };
    unsigned long N = 2;        /* number of cells, not nodes */
    struct py_field f = py_init_acoustic_field(N, N, x, y);

    /* used by the indexing macro */
    double *p = f.p.value;
    double *u = f.u.value;
    double *v = f.v.value;
    unsigned long nx = f.p.size_x;

    /* test 1 */
    for (unsigned long i = 0; i < N; ++i)
        for (unsigned long j = 0; j < N; ++j)
            P(i, j) = 1;
    py_leapfrog(&f);
    for (unsigned long i = 0; i < N + 1; ++i)
        for (unsigned long j = 0; j < N; ++j)
            CU_ASSERT_DOUBLE_EQUAL(U(i, j), 0, TOL);
    for (unsigned long i = 0; i < N; ++i)
        for (unsigned long j = 0; j < N + 1; ++j)
            CU_ASSERT_DOUBLE_EQUAL(V(i, j), 0, TOL);

    /* test 2 */
    py_apply_func(&f.p, zero2d);
    for (unsigned long i = 0; i < N + 1; ++i)
        for (unsigned long j = 0; j < N; ++j)
            U(i, j) = 1;
    for (unsigned long i = 0; i < N; ++i)
        for (unsigned long j = 0; j < N + 1; ++j)
            V(i, j) = 1;
    py_leapfrog(&f);
    for (unsigned long i = 0; i < N; ++i)
        for (unsigned long j = 0; j < N; ++j)
            CU_ASSERT_DOUBLE_EQUAL(P(i, j), 0, TOL);

    py_free_acoustic_field(f);
}

void test_partition_grid(void)
{
    struct py_cell_partition *part;
    unsigned long cells_x;
    unsigned long threads;

    cells_x = 8;
    threads = 1;
    part = py_partition_grid(threads, cells_x);
    CU_ASSERT(part[0].begin == 0 && part[0].end == 7);
    free(part);

    cells_x = 8;
    threads = 2;
    part = py_partition_grid(threads, cells_x);
    CU_ASSERT(part[0].begin == 0 && part[0].end == 3);
    CU_ASSERT(part[1].begin == 4 && part[1].end == 7);
    free(part);

    cells_x = 8;
    threads = 4;
    part = py_partition_grid(threads, cells_x);
    CU_ASSERT(part[0].begin == 0 && part[0].end == 1);
    CU_ASSERT(part[1].begin == 2 && part[1].end == 3);
    CU_ASSERT(part[2].begin == 4 && part[2].end == 5);
    CU_ASSERT(part[3].begin == 6 && part[3].end == 7);
    free(part);

    cells_x = 3;
    threads = 2;
    part = py_partition_grid(threads, cells_x);
    CU_ASSERT(part[0].begin == 0 && part[0].end == 1);
    CU_ASSERT(part[1].begin == 2 && part[1].end == 2);
    free(part);

    cells_x = 9;
    threads = 2;
    part = py_partition_grid(threads, cells_x);
    CU_ASSERT(part[0].begin == 0 && part[0].end == 4);
    CU_ASSERT(part[1].begin == 5 && part[1].end == 8);
    free(part);

    cells_x = 7;
    threads = 3;
    part = py_partition_grid(threads, cells_x);
    CU_ASSERT(part[0].begin == 0 && part[0].end == 2);
    CU_ASSERT(part[1].begin == 3 && part[1].end == 5);
    CU_ASSERT(part[2].begin == 6 && part[2].end == 6);
    free(part);
}

void test_get_partition_coords(void)
{
    double x[] = { 0, 1 };
    double y[] = { 0, 1 };
    double y_part[2];
    long threads, n;
    struct py_field f;
    struct py_cell_partition *part;

    n = 4;
    threads = 1;
    f = py_init_acoustic_field(n, n, x, y);
    part = py_partition_grid(threads, n);
    py_get_partition_coords(part[0], &f, y_part);
    CU_ASSERT_DOUBLE_EQUAL(y_part[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(y_part[1], 1, TOL);
    py_free_acoustic_field(f);

    n = 4;
    threads = 2;
    f = py_init_acoustic_field(n, n, x, y);
    part = py_partition_grid(threads, n);
    py_get_partition_coords(part[0], &f, y_part);
    CU_ASSERT_DOUBLE_EQUAL(y_part[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(y_part[1], 0.5, TOL);
    py_free_acoustic_field(f);

    n = 4;
    threads = 2;
    f = py_init_acoustic_field(n, n, x, y);
    part = py_partition_grid(threads, n);
    py_get_partition_coords(part[1], &f, y_part);
    CU_ASSERT_DOUBLE_EQUAL(y_part[0], 0.5, TOL);
    CU_ASSERT_DOUBLE_EQUAL(y_part[1], 1, TOL);
    py_free_acoustic_field(f);

    n = 4;
    threads = 4;
    f = py_init_acoustic_field(n, n, x, y);
    part = py_partition_grid(threads, n);
    py_get_partition_coords(part[0], &f, y_part);
    CU_ASSERT_DOUBLE_EQUAL(y_part[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(y_part[1], 0.25, TOL);
    py_free_acoustic_field(f);

    f = py_init_acoustic_field(n, n, x, y);
    part = py_partition_grid(threads, n);
    py_get_partition_coords(part[1], &f, y_part);
    CU_ASSERT_DOUBLE_EQUAL(y_part[0], 0.25, TOL);
    CU_ASSERT_DOUBLE_EQUAL(y_part[1], 0.5, TOL);
    py_free_acoustic_field(f);

    f = py_init_acoustic_field(n, n, x, y);
    part = py_partition_grid(threads, n);
    py_get_partition_coords(part[2], &f, y_part);
    CU_ASSERT_DOUBLE_EQUAL(y_part[0], 0.5, TOL);
    CU_ASSERT_DOUBLE_EQUAL(y_part[1], 0.75, TOL);
    py_free_acoustic_field(f);

    f = py_init_acoustic_field(n, n, x, y);
    part = py_partition_grid(threads, n);
    py_get_partition_coords(part[3], &f, y_part);
    CU_ASSERT_DOUBLE_EQUAL(y_part[0], 0.75, TOL);
    CU_ASSERT_DOUBLE_EQUAL(y_part[1], 1, TOL);
    py_free_acoustic_field(f);
}

void test_parse_cmdline(void)
{
    long nx = 0;
    long threads = 0;
    char outfile[STR_SIZE];
    int argc = 6;
    char *argv[] = { "/usr/bin/yee", "-n", "8", "-t", "4", "-q" };
    int write = 1;
    parse_cmdline(&nx, &threads, outfile, &write, argc, argv);
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
        || !CU_add_test(pSuite, "init_local_acoustic_field",
                        test_init_local_acoustic_field)
        || !CU_add_test(pSuite, "assign_to, get_from", test_assign_and_get)
        || !CU_add_test(pSuite, "set_boundary", test_set_boundary)
        || !CU_add_test(pSuite, "leapfrog", test_leapfrog)
        || !CU_add_test(pSuite, "partition_grid", test_partition_grid)
        || !CU_add_test(pSuite, "get_partition_coords",
                        test_get_partition_coords)
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
