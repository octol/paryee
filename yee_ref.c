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

void parse_cmdline(unsigned long *nx, char *outfile_p, char *outfile_u,
                   int argc, char *argv[])
{
    int opt;
    while ((opt = getopt(argc, argv, "n:p:u:")) != -1) {
        switch (opt) {
        case 'n':
            *nx = atoi(optarg);
            break;
        case 'p':
            strncpy(outfile_p, optarg, STR_SIZE);
            outfile_p[STR_SIZE - 1] = '\0';     /* force null termination */
            break;
        case 'u':
            strncpy(outfile_u, optarg, STR_SIZE);
            outfile_u[STR_SIZE - 1] = '\0';     /* force null termination */
            break;
        default:               /* '?' */
            fprintf(stderr, "Usage: %s ", argv[0]);
            fprintf(stderr, "[-n intervals] ");
            fprintf(stderr, "[-p outfile_p] [-u outfile_u]\n");
            exit(EXIT_FAILURE);
        }
    }
}

int write_to_disk(double *f, double *f_x, unsigned long size, char *str)
{
    /* append file suffix */
    char fstr[256];
    FILE *fp;
    unsigned int i;
    strcpy(fstr, str);
    /*strcat(fstr,".tsv"); */

    printf("Writing to: %s\n", fstr);

    fp = fopen(fstr, "w");
    if (fp == NULL) {
        perror("Error: can not write to disk");
        return EXIT_FAILURE;
    }
    for (i = 0; i < size; ++i)
        fprintf(fp, "%e\t%e\n", f_x[i], f[i]);
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
    double length = 1, cfl = 1, T = 0.3, c = 1;
    unsigned long nx = 2048, n, i;
    double dx, dt, Nt;
    double tic, toc;
    double *p;
    double *u;
    double *p_x;
    double *u_x;
    char outfile_p[STR_SIZE] = "yee_ref_p.tsv";
    char outfile_u[STR_SIZE] = "yee_ref_u.tsv";

    /* Parse parameters from commandline */
    parse_cmdline(&nx, outfile_p, outfile_u, argc, argv);
    printf("Running with: N=%lu\n", nx);

    p = malloc(sizeof(double) * nx);
    u = malloc(sizeof(double) * (nx + 1));
    p_x = malloc(sizeof(double) * nx);
    u_x = malloc(sizeof(double) * (nx + 1));
    if (!p || !u || !p_x || !u_x) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    dx = length / nx;
    dt = cfl * dx / c;
    Nt = T / dt;

    /* set grid (coordinates) and initial data */
    for (i = 0; i < nx; ++i) {
        u_x[i] = i * dx;
        p_x[i] = (i + 0.5) * dx;
        u[i] = 0;
        p[i] = exp(-pow(p_x[i] - 0.5, 2) / pow(0.05, 2));
    }
    u_x[nx] = nx * dx;
    u[nx] = 0;

    /* timestep */
    tic = gettime();
    for (n = 0; n < Nt; ++n) {
        /* update p */
        for (i = 0; i < nx; ++i)
            p[i] += dt / dx * (u[i + 1] - u[i]);

        /* update u */
        for (i = 1; i < nx - 1; ++i)
            u[i] += dt / dx * (p[i] - p[i - 1]);
    }
    toc = gettime();
    printf("Elapsed: %f seconds\n", toc - tic);

    /* write data to disk and free data */
    write_to_disk(p, p_x, nx, outfile_p);
    write_to_disk(u, u_x, nx + 1, outfile_u);
    free(p);
    free(u);
    free(p_x);
    free(u_x);

    return EXIT_SUCCESS;
}
