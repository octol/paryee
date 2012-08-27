#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mpi.h"

#include "yee_common.h"

#define MASTER 0
#define MAXWORKER 4
#define MINWORKER 1

int main(int argc, char* argv[])
{
    /* MPI */
    int rc, numtasks, taskid, numworkers;

    /* Paramters */
    double length=1, cfl=1, T=0.3, c=1;
    unsigned long nx=2048;
    unsigned int n, i, cells_per_workers;
    double tic, toc;
    Field f;
    FieldPartition* part;

    /* Parse parameters from commandline */
    parse_cmdline (&nx, NULL, argc, argv);
    printf("Running with: N=%ld\n", nx);

    /* Initialize */
    /*alloc_field(&f.p, nx);*/
    /*alloc_field(&f.u, nx+1);*/
    /*set_grid (&f.p, 0.5*length/nx, length-0.5*length/nx);*/
    /*set_grid (&f.u, 0, length);*/
    /*field_func (&f.p, gauss); [> initial data <]*/
    /*field_func (&f.u, zero);  [> initial data <]*/

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
            printf("ERROR: number of tasks must be between %d and %d.\n",
                    MINWORKER+1,MAXWORKER+1);
            printf("Quitting...\n");
            MPI_Abort(MPI_COMM_WORLD, rc);
            exit(1);
        }


    } else {
        /********************* Worker code *********************/ 
    }

    printf ("numtasks: %i\n", numtasks);
    /*printf ("taskid: %i\n", taskid);*/

    MPI_Finalize ();
    return 0;
}
