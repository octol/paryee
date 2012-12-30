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
 * Basic MPI distributed memory parallelizaton
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <unistd.h> // sleep

#include "mpi.h"

#include "yee_common.h"

#define BEGIN 1                 /* message tag */
#define COLLECT 2               /* message tag */
#define UTAG 3                  /* communication tag */
#define PTAG 4                  /* communication tag */

void send_param_data(int taskid, long left, long right, struct cell_partition part)
{
    /* Send out neighbour information */
    MPI_Send(&left, 1, MPI_LONG, taskid, BEGIN, MPI_COMM_WORLD);
    MPI_Send(&right, 1, MPI_LONG, taskid, BEGIN, MPI_COMM_WORLD);

    /* Send out partition information */
    MPI_Send(&part.begin, 1, MPI_LONG, taskid, BEGIN, MPI_COMM_WORLD);
    MPI_Send(&part.end, 1, MPI_LONG, taskid, BEGIN, MPI_COMM_WORLD);

#ifdef DEBUG
    printf("Sent to task %i:  ", taskid);
    printf("left=%li  right=%li  ", left, right);
    printf("begin=%li  end=%li\n", part.begin, part.end);
#endif
}


void receive_param_data(long *left, long *right, long *begin, long *end, MPI_Status *status)
{
    /* Data on neighbours */
    MPI_Recv(left, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, status);
    MPI_Recv(right, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, status);

    /* Partition information */
    MPI_Recv(begin, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, status);
    MPI_Recv(end, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, status);
}

void send_grid_data(int taskid, struct field *f)
{
    MPI_Send(&f->p.dx, 1, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);
    /*MPI_Send(&f->u.dx, 1, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->v.dx, 1, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->p.dy, 1, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->u.dy, 1, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->v.dy, 1, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->p.x[0], f->p.size_x, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->u.x[0], f->u.size_x, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->v.x[0], f->v.size_x, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->p.y[begin_yp], size_yp, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->u.y[begin_yu], size_yu, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->v.y[begin_yv], size_yv, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/

#ifdef DEBUG
    printf("Sent to task %i:  ", taskid);
    printf("dx=%f\n", f->p.dx);
#endif
}

void receive_grid_data(struct field *f, MPI_Status *status)
{
    MPI_Recv(&f->p.dx, 1, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, status);
    /*MPI_Recv(&f.u.dx, 1, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);*/
}

/*void send_field_data()*/
/*{*/
    /* Send out field data.
     * Note that MASTER deals with the taskid=0 special case, i.e., we do
     * not have to worry about it here. */
    /*long begin = part.begin;    [> for convenience <]*/
    /*long size = part.end - part.begin;*/

    /*long nx = f->p.size_x;*/
    /*long ny = f->p.size_y;*/
    /*long begin_p = begin*nx;*/
    /*long begin_u = begin*(nx + 1);*/
    /*long begin_v = begin*nx;*/
    
    /*MPI_Send(&f->p.value[begin_p], size, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->u.value[begin_u], size, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/
    /*MPI_Send(&f->v.value[begin_v], size, MPI_DOUBLE, taskid, BEGIN, MPI_COMM_WORLD);*/

    /* Send out grid data. 
     * This needs to take into consideration that the clients need to
     * know the grid data for the nodes that they receive data for. */
    /*long begin_yp = begin - 1;*/
    /*long begin_yu = begin - 1;*/
    /*long begin_yv = begin;*/
    /*long size_yp = size + 1;*/
    /*long size_yu = size + 1;*/
    /*long size_yv = size + 1;*/

/*}*/

void receive_field_data()
{
    /* Field data */
    /*MPI_Recv(f->p.value[begin, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, status);*/
    /*MPI_Recv(end, 1, MPI_LONG, MASTER, BEGIN, MPI_COMM_WORLD, status);*/
}

void collect_data(int taskid, struct cell_partition part, struct field *f, MPI_Status *status)
{
    long begin = part.begin;    /* for convenience */
    long size = part.end - part.begin;

    MPI_Recv(&f->p.value[begin], size, MPI_DOUBLE, taskid, COLLECT, MPI_COMM_WORLD, status);
    MPI_Recv(&f->u.value[begin], size, MPI_DOUBLE, taskid, COLLECT, MPI_COMM_WORLD, status);
    MPI_Recv(&f->v.value[begin], size, MPI_DOUBLE, taskid, COLLECT, MPI_COMM_WORLD, status);

#ifdef DEBUG
    printf("Results collected from task %d\n", taskid);
#endif
}

/* 
 * Main program 
 */
int main(int argc, char *argv[])
{
    /* Parameters */
    double x[] = { 0, 1 };
    double y[] = { 0, 1 };
    double cfl = 0.99 / sqrt(2);        /* CFL condition: c*dt/dx = cfl <= 1/sqrt(2) */
    double T = 0.3;
    double c = 1;
    long nx = 32;
    long ny = 32;
    char outfile[STR_SIZE] = "yee_mpi.tsv";
    long left = 0;
    long right = 0;

    /* Parse parameters from commandline */
    parse_cmdline(&nx, NULL, outfile, argc, argv);
    ny = nx;                    /* square domain */
    printf("Domain: %li x %li\n", nx, ny);

    /* Setup MPI and initialize variables */
    MPI_Status status;
    int rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    int numtasks;
    int taskid;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    /* Field variables - data structures for the simulation */
    /* Contains field data (pressure and velocity) */
    struct field f;
    /* Array of domain partitions */
    struct cell_partition *part;

    if (taskid == MASTER) {
        /********************* Master code *********************/
        /* MASTER deals with the setting up, sending out, receiving the result,
         * writing to disk. 
         * It also deals with the bottom partition. */
        
        /* Initialize */
        printf("MPI processes: %d\n", numtasks);
        f = init_acoustic_field(nx, ny, x, y);
        apply_func(&f.p, gauss2d);      /* initial data */

        /* Compute partition, neighbours and then send out data to the
         * workers */
        /* NOTE: partitioning across ny! (not nx) */
        part = partition_grid(numtasks, ny);

        for (long taskid = 1; taskid < numtasks; ++taskid) {
            /* compute neighbours */
            left = taskid - 1;
            right = (taskid == numtasks - 1) ? NONE : taskid + 1;

            /* send out data to the other workers / computational nodes */
            send_param_data(taskid, left, right, part[taskid]);
            send_grid_data(taskid, &f);
            /*send_field_data();*/
        }

        /* timing */
        /*double tic, toc;*/
        /*tic = gettime();*/

        /*[> Compute on the leftmost partition.  <]*/

        /*[> Wait for returning data from the workers <]*/
        /*for (long taskid = 1; taskid < numtasks; ++taskid) {*/
            /*collect_data(taskid, part[taskid], &f, &status);*/
        /*}*/

        /*toc = gettime();*/
        /*printf("Elapsed: %f seconds\n", toc - tic);*/

        sleep(2);

        write_to_disk(f.p, outfile);
        free(part);
        free_acoustic_field(f);
        MPI_Finalize();

    } else {
        /********************* Worker code *********************/
        /*part = malloc(sizeof(struct cell_partition));*/
        /*if (!part) {*/
            /*fprintf(stderr, "Memory allocation failed\n");*/
            /*exit(EXIT_FAILURE);*/
        /*}*/

        /* Receive grid and node parameters from master */

        /* Partition data */
        long begin, end;
        receive_param_data(&left, &right, &begin, &end, &status);
        long size = end - begin;

#ifdef DEBUG
        printf("Task %i received:  ", taskid);
        printf("left=%li  right=%li  ", left, right);
        printf("begin=%li  end=%li\n", begin, end);
#endif

        receive_grid_data(&f, &status);

        /* Given the partition data we compute the corresponding local array
         * indices in the array that is expanded to also contain the
         * neighbouring points needed to update. */
        /*long begin_p = 1*nx;*/
        /*long begin_u = 0*(nx+1);*/
        /*long begin_v = 1*nx;*/
        /*long size_p = (size + 1)*nx;*/
        /*long size_u = (size + 1)*nx;*/
        /*long size_v = (size + 1)*nx;*/

/*#ifdef DEBUG*/
        /*printf("taskid=%i:", taskid);*/
        /*printf("begin_p=%li  begin_u=%li  begin_v=%li  ", begin_p, begin_u, begin_v);*/
        /*printf("size_p=%li  size_u=%li  size_v=%li  ", size_p, size_u, size_v);*/
/*#endif*/

        /* Now that we have the grid info we can setup the data structures
         * containing the field. */
        /*alloc_field(&f.p, lsp);*/
        /*alloc_field(&f.u, lsu);*/
        /*apply_func(&f.p, zero); [> set everything to zero <]*/
        /*apply_func(&f.u, zero); [> set everything to zero <]*/

        /*[> Field data <]*/
        /*[> NOTE: we only recieve sp/su number of values, not lsp/lsu. <]*/
        /*MPI_Recv(&f.p.value[lbp], sp, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);*/
        /*MPI_Recv(&f.u.value[lbu], su, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);*/
        /*[> Grid data <]*/
        /*MPI_Recv(&f.p.x[lbp], sp, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);*/
        /*MPI_Recv(&f.u.x[lbu], su, MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);*/
/*#ifdef DEBUG*/
        /*printf("Task=%d received:  left=%ld  right=%ld", taskid, left, right);*/
        /*printf("  bp=%ld  ep=%ld", bp, ep);*/
        /*printf("  bu=%ld  eu=%ld\n", bu, eu);*/
/*#endif*/

        /*[> Depends on the numerical variables initialized above <]*/
        /*[> CFL condition is: c*dt/dx = cfl <= 1 <]*/
        /*f.dt = cfl * f.p.dx / c;*/
        /*f.Nt = T / f.dt;*/

        /*long i;*/
        /*for (long n = 0; n < f.Nt; ++n) {*/
            /*[> Communicate <]*/
            /*[> send u to the left <]*/
            /*[> receive u from the right <]*/
            /*if (left != NONE) {*/
                /*MPI_Send(&f.u.value[lbu], 1, MPI_DOUBLE, left, UTAG, MPI_COMM_WORLD);*/
            /*}*/
            /*if (right != NONE) {*/
                /*MPI_Recv(&f.u.value[leu + 1], 1, MPI_DOUBLE, right, UTAG, MPI_COMM_WORLD, &status);*/
            /*}*/

            /*[> update the pressure (p) <]*/
            /*[>update_field_s(&f.p, lbp, sp, &f.u, 0, f.dt); <]*/
            /*for (i = lbp; i <= lbp + sp - 1; ++i)*/
                /*f.p.value[i] += f.dt / f.p.dx*/
                    /** (f.u.value[i + 1 - lbp] - f.u.value[i - lbp]);*/

            /*[> Communicate <]*/
            /*[> receive p from the left <]*/
            /*[> send p to the right <]*/
            /*if (left != NONE) {*/
                /*MPI_Recv(&f.p.value[lbp - 1], 1, MPI_DOUBLE, left, PTAG, MPI_COMM_WORLD, &status);*/
            /*}*/
            /*if (right != NONE) {*/
                /*MPI_Send(&f.p.value[lep], 1, MPI_DOUBLE, right, PTAG, MPI_COMM_WORLD);*/
            /*}*/

            /*[> update the velocity (u) <]*/
            /*[>update_field_s(&f.u, lbu, su, &f.p, 0, f.dt); <]*/
            /*for (i = lbu; i <= lbu + su - 1; ++i)*/
                /*f.u.value[i] += f.dt / f.p.dx*/
                    /** (f.p.value[i - lbu + 1] - f.p.value[i - lbu]);*/
        /*}*/

        /*[> Send back data to master <]*/
        /*MPI_Send(&f.p.value[lbp], sp, MPI_DOUBLE, MASTER, COLLECT, MPI_COMM_WORLD);*/
        /*MPI_Send(&f.u.value[lbu], su, MPI_DOUBLE, MASTER, COLLECT, MPI_COMM_WORLD);*/

        /* Remember to deallocate */
        /*free(part);*/
        /*free_acoustic_field(f);*/

        MPI_Finalize();
    }

    return 0;
}
