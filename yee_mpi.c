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
#define MAXWORKER 8
#define NONE 0          /* no neighbour */
#define BEGIN 1         /* message tag */
#define COLLECT 2       /* message tag */
#define LTAG 3
#define RTAG 4

int main(int argc, char* argv[])
{
    /* MPI */
    int rc, numtasks, taskid, numworkers, left, right;
    MPI_Status status;

    /* Parameters */
    double length=1, cfl=1, T=0.3, c=1;
    unsigned long nx=2048;
    unsigned int n, i, cells_per_worker;
    Field f;
    FieldPartition* part;

    /* Parse parameters from commandline */
    parse_cmdline (&nx, NULL, argc, argv);

    /* Setup MPI and initialize variables */
    rc = MPI_Init (&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program. Terminating.\n");
        MPI_Abort (MPI_COMM_WORLD, rc);
    }
    MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank (MPI_COMM_WORLD, &taskid);
    numworkers = numtasks - 1;

    if (taskid == MASTER) {
        /********************* Master code *********************/ 
        if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)) {
            printf ("ERROR: number of tasks must be between %d and %d.\n",
                    MINWORKER+1, MAXWORKER+1);
            printf ("Quitting...\n");
            MPI_Abort (MPI_COMM_WORLD, rc);
            exit (EXIT_FAILURE);
        }

        printf ("Starting yee_mpi with %i workers\n",numworkers);
        printf ("Running with: N=%ld\n", nx);

        /* Initialize */
        alloc_field (&f.p, nx);
        alloc_field (&f.u, nx+1);
        set_grid (&f.p, 0.5*length/nx, length-0.5*length/nx);
        set_grid (&f.u, 0, length);
        field_func (&f.p, gauss); /* initial data */
        field_func (&f.u, zero);  /* initial data */

        /* Compute partition, neighbours and then send out data to the
         * workers */
        cells_per_worker = nx/numworkers;
        assert (cells_per_worker*numworkers == nx);
        /*printf ("nx=%i  numworkers=%i  cells_per_worker=%i\n",*/
        /*        nx,numworkers,cells_per_worker);              */
        part = malloc (sizeof(FieldPartition)*numworkers);
        if (!part) {
            fprintf (stderr,"Memory allocation failed\n");
            exit (EXIT_FAILURE);
        }
        /* remember that i=0 is the master */
        for (i=1; i<=numworkers; ++i) {
            int bp,bu,ep,eu,sp,su;  /* array indices and sizes */

            /* partition_grid() requires the nodes to be enumerated starting
             * from 0: 0,...,total_nodes */ 
            part[i-1] = partition_grid_into_cells (i-1, numworkers, 
                                                   cells_per_worker);
            /* Make field data more compact, so we can more easily use as 
             * array indices. For convenience only. */
            bp = part[i-1].start_p; /* begin index */
            bu = part[i-1].start_u; 
            ep = part[i-1].end_p;   /* end index */
            eu = part[i-1].end_u; 
            sp = ep - bp + 1;       /* size */
            su = eu - bu + 1;  

            /* compute neighbours */
            left = (i==1) ? NONE : i-1;
            right = (i==numworkers) ? NONE : i+1;

            /* Send out neighbour information */
            MPI_Send (&left, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&right,1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);

            /* Send out partition information */
            MPI_Send (&bp, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&ep, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&bu, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&eu, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);

            /* Send out field data */
            MPI_Send (&f.p.value[bp],sp,MPI_DOUBLE,i,BEGIN,MPI_COMM_WORLD);
            MPI_Send (&f.u.value[bu],su,MPI_DOUBLE,i,BEGIN,MPI_COMM_WORLD);

            /* Send out grid data */
            MPI_Send (&f.p.x[bp], sp, MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&f.u.x[bu], su, MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&f.p.dx, 1, MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);
            MPI_Send (&f.u.dx, 1, MPI_DOUBLE, i, BEGIN, MPI_COMM_WORLD);

            printf ("Sent to task %d: left=%d, right=%d\n",i,left,right);
        }

        /* Wait for returning data from the workers */
        for (i=1; i<=numworkers; ++i) {
            /* Make field data more compact, so we can more easily use as 
             * array indices. For convenience only. */
            int bp,bu,ep,eu,sp,su;  /* array indices and sizes */
            bp = part[i-1].start_p; /* begin index */
            bu = part[i-1].start_u; 
            ep = part[i-1].end_p;   /* end index */
            eu = part[i-1].end_u; 
            sp = ep - bp + 1;       /* size */
            su = eu - bu + 1;  

            MPI_Recv (&f.p.value[bp],sp, MPI_DOUBLE, i, COLLECT, 
                      MPI_COMM_WORLD, &status);
            MPI_Recv (&f.u.value[bu],su, MPI_DOUBLE, i, COLLECT, 
                      MPI_COMM_WORLD, &status);
            printf ("Results collected from task %d\n",i);
        }

        write_to_disk(f.p, "output_p"); 
        write_to_disk(f.u, "output_u"); 

        free_field(f.p);
        free_field(f.u);

        MPI_Finalize ();

    } else {
        /********************* Worker code *********************/ 
        int bp,bu,ep,eu,sp,su;  /* array indices and sizes */
        i = taskid;

        part = malloc (sizeof(FieldPartition));
        if (!part) {
            fprintf(stderr,"Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }

        /* Receive grid and node parameters from master */
        
        /* Data on neighbours */
        MPI_Recv (&left, 1, MPI_INT, MASTER,BEGIN,MPI_COMM_WORLD,&status);
        MPI_Recv (&right,1, MPI_INT, MASTER,BEGIN,MPI_COMM_WORLD,&status);
        /* Partition data */
        MPI_Recv (&bp, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv (&ep, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv (&bu, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
        MPI_Recv (&eu, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);

        /* From indices we compute the sizes */
        sp = ep - bp + 1;   
        su = eu - bu + 1;  

        /* The expanded domain to store the additional data that we receive
         * from the neighbours */
        int sp_expanded = sp + 1;   /* p receives data from the left */
        int su_expanded = su + 1;   /* u receives data from the right */
        if (left == NONE) {
            /* p[-1] is not needed since we dont update u[0] */
            --sp_expanded;  
            /* add the additional bc point since p depends on it */
            ++su_expanded;  
        }

        /* Compute local indices. These specify the local indices of the
         * points we update in the expanded local array. */
        int lbp = 1;    /* since p is padded by one point to the left */
        int lbu = 0;    /* since u is not padded on the left */
        if (left == NONE) {
            lbp = 0;
            lbu = 1;
        }
        int lep = ep - bp + lbp;
        int leu = eu - bu + lbu;

        /* Now that we have the grid info we can setup the data structures
         * containing the field. */
        alloc_field (&f.p, sp_expanded);
        alloc_field (&f.u, su_expanded);
        field_func (&f.p, zero);  /* set everything to zero */
        field_func (&f.u, zero);  /* set everything to zero */

        /* Field data */
        MPI_Recv (&f.p.value[lbp], sp_expanded, MPI_DOUBLE, MASTER, BEGIN, 
                  MPI_COMM_WORLD, &status);
        MPI_Recv (&f.u.value[lbu], su_expanded, MPI_DOUBLE, MASTER, BEGIN, 
                  MPI_COMM_WORLD, &status);
        /* Grid data */
        MPI_Recv (&f.p.x[lbp], sp_expanded, MPI_DOUBLE, MASTER, BEGIN, 
                  MPI_COMM_WORLD, &status);
        MPI_Recv (&f.u.x[lbu], su_expanded, MPI_DOUBLE, MASTER, BEGIN,
                  MPI_COMM_WORLD, &status);
        MPI_Recv (&f.p.dx,1,MPI_DOUBLE,MASTER,BEGIN,MPI_COMM_WORLD,&status);
        MPI_Recv (&f.u.dx,1,MPI_DOUBLE,MASTER,BEGIN,MPI_COMM_WORLD,&status);
        printf ("Task=%d received:  left=%d  right=%d",taskid,left,right);
        printf ("  bp=%d  ep=%d",bp,ep);
        printf ("  bu=%d  eu=%d\n",bu,eu);

        /* Depends on the numerical variables initialized above */
        f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
        f.Nt = T/f.dt;

        for (n=0; n<f.Nt; ++n) {
            /*printf ("taking a timestep\n");*/
            /* Take one timestep */
            /*update_field (&f.p, 0, nx, &f.u, 0, f.dt);*/
            /*update_field (&f.u, 1, nx-1, &f.p, 0, f.dt);*/
            /*update_field (&f.p, bp, ep, &f.u, 0, f.dt);*/
            /*update_field (&f.u, bu, eu, &f.p, 0, f.dt);*/

            /* Communicate */
            /* send u to the left */
            /* receive u from the right */
TODO: get this right
            if (left != NONE) {
                MPI_Send (&f.u.value[lbu], sp, MPI_DOUBLE, left, LTAG, 
                          MPI_COMM_WORLD);
            }
            if (right != NONE) {
                MPI_Recv (&f.u.value[leu+1], sp, MPI_DOUBLE, left, LTAG, 
                          MPI_COMM_WORLD);
            }
            update_field (&f.p, lbp, lep+1, &f.u, 0, f.dt);

            /* Communicate */
            update_field (&f.u, lbu, leu+1, &f.p, 0, f.dt);
        }

        /* Send back data to master */
        MPI_Send (&f.p.value[bp], sp, MPI_DOUBLE, MASTER, COLLECT, 
                  MPI_COMM_WORLD);
        MPI_Send (&f.u.value[bu], su, MPI_DOUBLE, MASTER, COLLECT, 
                  MPI_COMM_WORLD);

        /* Remember to deallocate */
        free_field(f.p);
        free_field(f.u);

        MPI_Finalize ();
    }

    return 0;
}
