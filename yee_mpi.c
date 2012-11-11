/* 
 * Basic implementation of 1D wave equation on system form
 *
 *      p_t = a u_x
 *      u_t = b p_x
 *
 * This formulation is sometimes used in acoustics, where then p is the
 * pressure and u is the velocity field in the direction of the x-axis.
 * 
 * We use a staggared grid with u defined on the endpoints, and thus on full
 * integer indices. N denotes the number of interval lengths, i.e., we have
 * N+1 number of nodes for u
 *
 * The outer (x=0, x=L) boundary is set to u=0. This means we start the
 * enumeration in the second u point when we partition the domain. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mpi.h"

#include "yee_common.h"

#define MASTER 0        /* id of master process */
#define MINWORKER 1
#define MAXWORKER 4
#define NONE 0          /* no neighbour */
#define BEGIN 1         /* message tag */

int main(int argc, char* argv[])
{
    /* MPI */
    int rc, numtasks, taskid, numworkers, left, right;
    MPI_Status status;

    /* Paramters */
    double length=1, cfl=1, T=0.3, c=1;
    unsigned long nx=2048;
    unsigned int n, i, cells_per_worker;
    double tic, toc;
    Field f;
    FieldPartition* part;

    /* Parse parameters from commandline */
    parse_cmdline (&nx, NULL, argc, argv);

    /* Setup MPI and initialize variables */
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    numworkers = numtasks - 1;

    if (taskid == MASTER) {
        /********************* Master code *********************/ 
        if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)) {
            printf ("ERROR: number of tasks must be between %d and %d.\n",
                    MINWORKER+1,MAXWORKER+1);
            printf ("Quitting...\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(1);
        }

        printf ("Starting yee_mpi with %i workers\n",numworkers);
        printf ("Running with: N=%ld\n", nx);

        /* Initialize */
        alloc_field(&f.p, nx);
        alloc_field(&f.u, nx+1);
        set_grid (&f.p, 0.5*length/nx, length-0.5*length/nx);
        set_grid (&f.u, 0, length);
        field_func (&f.p, gauss); /* initial data */
        field_func (&f.u, zero);  /* initial data */

        /* Compute partition, neighbours and then send out data to the
         * workers */
        cells_per_worker = nx/numworkers;
        assert (cells_per_worker*numworkers == nx);
        part = malloc (sizeof(FieldPartition)*numworkers);
        if (!part) {
            fprintf(stderr,"Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (i=1; i<=numworkers; ++i) {
            /* partition_grid() required the nodes to be enumerated starting
             * from 0: 0,...,total_nodes */ 
            part[i-1] = partition_grid (i-1,numworkers, cells_per_worker);
            
            /* compute neighbours */
            if (i==1)
                left = NONE;
            else
                left = i-1;
            if (i==numworkers)
                right = NONE;
            else
                right = i+1;

            /* Send out data */
            MPI_Send (&left, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&right, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&part[i-1].start_p,1,MPI_INT,i,BEGIN, MPI_COMM_WORLD);
            MPI_Send (&part[i-1].end_p, 1, MPI_INT,i,BEGIN, MPI_COMM_WORLD);
            MPI_Send (&part[i-1].start_u,1,MPI_INT,i,BEGIN, MPI_COMM_WORLD);
            MPI_Send (&part[i-1].end_u,1,MPI_INT, i,BEGIN, MPI_COMM_WORLD);
            MPI_Send (&f.p.value[part[i-1].start_p], 
                      part[i-1].end_p-part[i-1].start_p+1,
                      MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&f.u.value[part[i-1].start_u], 
                      part[i-1].end_u-part[i-1].start_u+1,
                      MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);
            printf ("Sent to task %d: left=%d, right=%d\n",i,left,right);
        }

        /* Wait for returning data from the workers */

    } else {
        /********************* Worker code *********************/ 
        i = taskid;

        /* Initialize */
        alloc_field(&f.p, nx);
        alloc_field(&f.u, nx+1);
        field_func (&f.p, zero);  /* set everything to zero */
        field_func (&f.u, zero);  /* set everything to zero */

        part = malloc (sizeof(FieldPartition));
        if (!part) {
            fprintf(stderr,"Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }

        /* Receive data from master */
        MPI_Recv (&left, 1, MPI_INT, MASTER,BEGIN,MPI_COMM_WORLD,&status);
        MPI_Recv (&right, 1, MPI_INT,MASTER,BEGIN,MPI_COMM_WORLD,&status);
        MPI_Recv (&part->start_p,1,
                  MPI_INT,MASTER,BEGIN, MPI_COMM_WORLD,&status);
        MPI_Recv (&part->end_p, 1, 
                  MPI_INT,MASTER,BEGIN, MPI_COMM_WORLD,&status);
        MPI_Recv (&part->start_u,1,
                  MPI_INT,MASTER,BEGIN, MPI_COMM_WORLD,&status);
        MPI_Recv (&part->end_u,1, 
                  MPI_INT, MASTER,BEGIN, MPI_COMM_WORLD,&status);
        MPI_Recv (&f.p.value[part->start_p], 
                  part->end_p-part->start_p+1,
                  MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv (&f.u.value[part->start_u], 
                  part->end_u-part->start_u+1,
                  MPI_DOUBLE, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        printf ("Task %d received: left=%d, right=%d\n",taskid,left,right);

        /* Depends on the numerical variables initialized above */
        f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
        f.Nt = T/f.dt;

        for (n=0; n<f.Nt; ++n) {
            /* Take one timestep */

            /* Communicate */
        }

        /* Send back data to master */
    }

    /* Remember to deallocate */

    MPI_Finalize ();
    return 0;
}
