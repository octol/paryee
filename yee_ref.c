/*
 *    Copyright (C) 2012 Jon Haggblad
 *
 *    This file is part of ParYee.
 *
 *    ParYee is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ParYee is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ParYee.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Minimal reference implementation to compare overhead. This implementation
 * is completely independent, and does not depend on yee_common.c.
 */
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define STR_SIZE 256

#define p(x, y) (p[(x) + (y)*nx])
#define u(x, y) (u[(x) + (y)*(nx + 1)])
#define v(x, y) (v[(x) + (y)*nx])

void parse_cmdline(unsigned long *nx, char *outfile, int argc,
                   char *argv[])
{
    int opt;
    while ((opt = getopt(argc, argv, "n:o:")) != -1) {
        switch (opt) {
        case 'n':
            *nx = atoi(optarg);
            break;
        case 'o':
            strncpy(outfile, optarg, STR_SIZE);
            outfile[STR_SIZE - 1] = '\0';       /* force null termination */
            break;
        default:               /* '?' */
            fprintf(stderr, "Usage: %s ", argv[0]);
            fprintf(stderr, "[-n intervals] ");
            fprintf(stderr, "[-o outfile]\n");
            exit(EXIT_FAILURE);
        }
    }
}

int write_to_disk(double *f, double *f_x, double *f_y,
                  unsigned long size_x, unsigned long size_y, char *str)
{
    char fstr[256];
    FILE *fp;
    unsigned long i, j;
    strcpy(fstr, str);

    printf("Writing to: %s\n", fstr);

    fp = fopen(fstr, "w");
    if (fp == NULL) {
        perror("Error: can not write to disk");
        return EXIT_FAILURE;
    }
    for (i = 0; i < size_x; ++i)
        for (j = 0; j < size_y; ++j)
            fprintf(fp, "%e\t%e\t%e\n", f_x[i], f_y[i], f[i + j * size_x]);
    fclose(fp);
    return EXIT_SUCCESS;
}

double gettime(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1.0e-6 * t.tv_usec + t.tv_sec);
}

int main(int argc, char *argv[])
{
    /* Simulation parameters */
    double length = 1;
    double cfl = 1 / sqrt(2) * 0.99;
    double T = 0.3;
    double c = 1;
    unsigned long nx = 8;
    unsigned long ny = 8;

    unsigned long n, i, j;      /* time, x, y index */
    double dx, dy, dt, Nt;
    double tic, toc;
    double *p;
    double *u;
    double *v;
    double *p_x;
    double *p_y;
    double *u_x;
    double *u_y;
    double *v_x;
    double *v_y;
    char outfile[STR_SIZE] = "yee_ref.tsv";
    /*char outfile_u[STR_SIZE] = "yee_ref_u.tsv"; */
    /*char outfile_v[STR_SIZE] = "yee_ref_v.tsv"; */
    double radius_squared;      /* shorthand in computations */

    /* Parse parameters from commandline */
    parse_cmdline(&nx, outfile, argc, argv);
    printf("Running with: N=%lu\n", nx);

    p = malloc(sizeof(double) * nx * ny);
    u = malloc(sizeof(double) * (nx + 1) * ny);
    v = malloc(sizeof(double) * nx * (ny + 1));
    p_x = malloc(sizeof(double) * nx);
    p_y = malloc(sizeof(double) * ny);
    u_x = malloc(sizeof(double) * (nx + 1));
    u_y = malloc(sizeof(double) * ny);
    v_x = malloc(sizeof(double) * nx);
    v_y = malloc(sizeof(double) * (ny + 1));
    if (!p || !u || !v || !p_x || !p_y || !u_x || !u_y || !v_x || !v_y) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    dx = length / nx;
    dy = length / ny;
    dt = cfl * dx / c;
    Nt = T / dt;

    /* set grid (coordinates) */
    for (i = 0; i < nx; ++i) {
        p_x[i] = (i + 0.5) * dx;
        p_y[i] = (i + 0.5) * dy;
        u_x[i] = i * dx;
        u_y[i] = (i + 0.5) * dy;
        u_x[i] = (i + 0.5) * dx;
        u_y[i] = i * dy;
    }
    u_x[nx] = nx * dx;
    v_y[ny] = ny * dy;

    /* initial data */
    for (i = 0; i < nx; ++i) {
        for (j = 0; i < ny; ++j) {
            radius_squared = pow(p_x[i] - 0.5, 2) + pow(p_y[j] - 0.5, 2);
            p(i, j) = exp(-radius_squared / pow(0.05, 2));
        }
    }
    for (i = 0; i < nx + 1; ++i)
        for (j = 0; i < ny; ++j)
            u(i, j) = 0;
    for (i = 0; i < nx; ++i)
        for (j = 0; i < ny + 1; ++j)
            v(i, j) = 0;

    /* timestep */
    tic = gettime();
    for (n = 0; n < Nt; ++n) {

        /* update p */
        for (i = 0; i < nx; ++i) {
            for (j = 0; i < nx; ++j) {
                p(i, j) += dt * (1 / dx * (u(i + 1, j) - u(i, j))
                                 + 1 / dy * (v(i, j + 1) - v(i, j)));
            }
        }

        /* update u */
        for (i = 1; i < nx - 1; ++i) {
            for (j = 0; i < ny; ++j) {
                u(i, j) += dt / dx * (p(i, j) - p(i - 1, j));
            }
        }

        /* update v */
        for (i = 0; i < nx; ++i) {
            for (j = 1; i < ny - 1; ++j) {
                u(i, j) += dt / dx * (p(i, j) - p(i, j - 1));
            }
        }
    }
    toc = gettime();
    printf("Elapsed: %f seconds\n", toc - tic);

    /* write data to disk and free data */
    write_to_disk(p, p_x, p_y, nx, ny, outfile);
    /*write_to_disk(u, u_x, nx + 1, outfile_u); */
    free(p);
    free(u);
    free(v);
    free(p_x);
    free(u_x);
    free(v_x);
    free(p_y);
    free(u_y);
    free(v_y);

    return EXIT_SUCCESS;
}
