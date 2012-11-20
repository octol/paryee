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

#define MASTER 0            /* id of master process */
#define MINWORKER 1 
#define MAXWORKER 8
#define NONE 0              /* no neighbour */
#define BEGIN 1             /* message tag */
#define COLLECT 2           /* message tag */
#define UTAG 3              /* communication tag */
#define PTAG 4              /* communication tag */
#define MCW MPI_COMM_WORLD  /* abbreviation */

int main(int argc, char* argv[])
{
    /* Simulation parameters */
    double length = 1;
    double cfl = 1;
    double T = 0.3;
    double c = 1;

    /* Numerical paramters */
    unsigned long nx = 2048;

    int cells_per_worker = 0;
    unsigned int left = 0; 
    unsigned int right = 0;

    /* Parse parameters from commandline */
    parse_cmdline (&nx, NULL, argc, argv);

    /* Setup MPI and initialize variables */
    MPI_Status status;
    int rc = MPI_Init (&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program. Terminating.\n");
        MPI_Abort (MPI_COMM_WORLD, rc);
    }
    int numtasks;
    int taskid;
    MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank (MPI_COMM_WORLD, &taskid);
    int numworkers = numtasks - 1;

    /* Field variables - data structures for the simulation */
    Field f;            /* Contains field data (pressure and velocity) */
    Partition* part;    /* Array of domain partitions */

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
        apply_func (&f.p, gauss); /* initial data */
        apply_func (&f.u, zero);  /* initial data */

        /* Compute partition, neighbours and then send out data to the
         * workers */
        cells_per_worker = nx/numworkers;
        if (cells_per_worker*numworkers != nx) {
            fprintf (stderr,"Only 2^n workers supported\n");
            exit (EXIT_FAILURE);
        }
        part = malloc (sizeof(Partition)*numworkers);
        if (!part) {
            fprintf (stderr,"Memory allocation failed\n");
            exit (EXIT_FAILURE);
        }

        /* remember that tid=0 is the master */
        for (int tid=1; tid<=numworkers; ++tid) {
            /* partition_grid() requires the nodes to be enumerated starting
             * from 0: 0,...,total_nodes */ 
            part[tid-1]=partition_grid2(tid-1,numworkers,cells_per_worker);
            /* compute neighbours */
            left = (tid==1) ? NONE : tid-1;
            right = (tid==numworkers) ? NONE : tid+1;

            /* Send out neighbour information */
            MPI_Send (&left, 1, MPI_INT, tid, BEGIN, MCW);
            MPI_Send (&right,1, MPI_INT, tid, BEGIN, MCW);

            /* Make field data more compact, so we can more easily use as 
             * array indices. For convenience only. */
            int bp, ep, sp, bu, eu, su;
            expand_indices (part[tid-1], &bp, &ep, &sp, &bu, &eu, &su);

            /* Send out partition information */
            MPI_Send (&bp, 1, MPI_INT, tid, BEGIN, MCW);
            MPI_Send (&ep, 1, MPI_INT, tid, BEGIN, MCW);
            MPI_Send (&bu, 1, MPI_INT, tid, BEGIN, MCW);
            MPI_Send (&eu, 1, MPI_INT, tid, BEGIN, MCW);

            /* Send out field data */
            MPI_Send(&f.p.value[bp],sp,MPI_DOUBLE,tid,BEGIN,MCW);
            MPI_Send(&f.u.value[bu],su,MPI_DOUBLE,tid,BEGIN,MCW);

            /* Send out grid data */
            MPI_Send (&f.p.x[bp], sp, MPI_DOUBLE, tid,BEGIN,MCW);
            MPI_Send (&f.u.x[bu], su, MPI_DOUBLE, tid,BEGIN,MCW);
            MPI_Send (&f.p.dx, 1, MPI_DOUBLE, tid, BEGIN, MCW);
            MPI_Send (&f.u.dx, 1, MPI_DOUBLE, tid, BEGIN, MCW);

#ifdef DEBUG
            printf ("Sent to task %d: left=%d, right=%d\n",i,left,right);
#endif
        }

        /* Wait for returning data from the workers */
        for (int tid=1; tid<=numworkers; ++tid) {
            /* Make field data more compact, so we can more easily use as 
             * array indices. For convenience only. */
            int bp, ep, sp, bu, eu, su;
            expand_indices (part[tid-1], &bp, &ep, &sp, &bu, &eu, &su);

            MPI_Recv (&f.p.value[bp],sp,MPI_DOUBLE,tid,COLLECT,MCW,&status);
            MPI_Recv (&f.u.value[bu],su,MPI_DOUBLE,tid,COLLECT,MCW,&status);
#ifdef DEBUG
            printf ("Results collected from task %d\n",i);
#endif
        }

        write_to_disk(f.p, "output_p"); 
        write_to_disk(f.u, "output_u"); 

        free (part);
        free_field(f.p);
        free_field(f.u);

        MPI_Finalize ();

    } else {
        /********************* Worker code *********************/ 
        int bp,bu,ep,eu,sp,su;  /* array indices and sizes */
        int lbp,lbu,lep,leu,lsp,lsu;  /* local array indices and sizes */
        /*int i = taskid;*/

        part = malloc (sizeof(Partition));
        if (!part) {
            fprintf(stderr,"Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }

        /* Receive grid and node parameters from master */
        
        /* Data on neighbours */
        MPI_Recv (&left, 1, MPI_INT, MASTER,BEGIN,MCW,&status);
        MPI_Recv (&right,1, MPI_INT, MASTER,BEGIN,MCW,&status);
        /* Partition data */
        MPI_Recv (&bp, 1, MPI_INT, MASTER, BEGIN, MCW, &status);
        MPI_Recv (&ep, 1, MPI_INT, MASTER, BEGIN, MCW, &status);
        MPI_Recv (&bu, 1, MPI_INT, MASTER, BEGIN, MCW, &status);
        MPI_Recv (&eu, 1, MPI_INT, MASTER, BEGIN, MCW, &status);

        /* From indices we compute the sizes */
        sp = ep - bp + 1;   
        su = eu - bu + 1;  

        /* Given the partition data we compute the corresponding local array
         * indices in the array that is expanded to also contain the
         * neighbouring points needed to update. */
        lbp = 1;    /* since p is padded by one point to the left */
        lbu = 0;    /* since u is not padded on the left */
        if (left == NONE) {
            lbp = 0;
            lbu = 1;
        }
        lep = lbp + sp - 1;
        leu = lbu + su - 1;
        /* sizes */
        lsp = lep + 1;
        lsu = leu + 2; 

#ifdef DEBUG
        printf ("taskid=%i", taskid);  
        printf ("bp=%i  ep=%i  bu=%i  eu=%i  ", bp, ep, bu, eu);
        printf ("sp=%i  su=%i  ",sp, su);
        printf ("lbp=%i  lep=%i  lbu=%i  leu=%i  ", lbp, lep, lbu, leu);
        printf ("lsp=%i  lsu=%i\n",lsp, lsu);
#endif

        /* Now that we have the grid info we can setup the data structures
         * containing the field. */
        alloc_field (&f.p, lsp);
        alloc_field (&f.u, lsu);
        field_func (&f.p, zero);  /* set everything to zero */
        field_func (&f.u, zero);  /* set everything to zero */

        /* Field data */
        /* NOTE: we only recieve sp/su number of values, not lsp/lsu. */
        MPI_Recv (&f.p.value[lbp], sp, MPI_DOUBLE,MASTER,BEGIN,MCW,&status);
        MPI_Recv (&f.u.value[lbu], su, MPI_DOUBLE,MASTER,BEGIN,MCW,&status);
        /* Grid data */
        MPI_Recv (&f.p.x[lbp], sp, MPI_DOUBLE, MASTER, BEGIN, MCW, &status);
        MPI_Recv (&f.u.x[lbu], su, MPI_DOUBLE, MASTER, BEGIN, MCW, &status);
        MPI_Recv (&f.p.dx,1,MPI_DOUBLE,MASTER,BEGIN,MCW,&status);
        MPI_Recv (&f.u.dx,1,MPI_DOUBLE,MASTER,BEGIN,MCW,&status);
#ifdef DEBUG
        printf ("Task=%d received:  left=%d  right=%d",taskid,left,right);
        printf ("  bp=%d  ep=%d",bp,ep);
        printf ("  bu=%d  eu=%d\n",bu,eu);
#endif

        /* Depends on the numerical variables initialized above */
        f.dt = cfl*f.p.dx/c; /* CFL condition is: c*dt/dx = cfl <= 1 */
        f.Nt = T/f.dt;

        for (unsigned long n=0; n<f.Nt; ++n) {
            /* Communicate */
            /* send u to the left */
            /* receive u from the right */
            if (left != NONE) {
                MPI_Send (&f.u.value[lbu], 1, MPI_DOUBLE, left, UTAG, MCW);
            }
            if (right != NONE) {
                MPI_Recv (&f.u.value[leu+1],1,MPI_DOUBLE,right,UTAG,MCW,&status);
            }
            /*update_field (&f.p, lbp, lep+1, &f.u, 0, f.dt);*/
            update_field2 (&f.p, lbp, sp, &f.u, 0, f.dt);
            /*update_field3 (&f.p, lbp, lep, &f.u, 0, f.dt);*/
            /*printf ("lbp=%i  lep=%i  lbp+sp=%i\n",lbp,lep,lbp+sp);*/

            /* Communicate */
            /* receive p from the left */
            /* send p to the right */
            if (left != NONE) {
                MPI_Recv (&f.p.value[lbp-1],1,MPI_DOUBLE,left,PTAG,MCW,&status);
            }
            if (right != NONE) {
                MPI_Send (&f.p.value[lep], 1, MPI_DOUBLE, right, PTAG, MCW);
            }
            /*update_field (&f.u, lbu, leu+1, &f.p, 0, f.dt);*/
            update_field2 (&f.u, lbu, su, &f.p, 0, f.dt);
            /*update_field3 (&f.u, lbu, leu, &f.p, 0, f.dt);*/
            /*printf ("lbu=%i  leu=%i  lbu+su=%i\n",lbu,leu,lbu+su);*/
        }

        /* Send back data to master */
        MPI_Send (&f.p.value[lbp], sp, MPI_DOUBLE, MASTER, COLLECT, MCW);
        MPI_Send (&f.u.value[lbu], su, MPI_DOUBLE, MASTER, COLLECT, MCW);

        /* Remember to deallocate */
        free (part);
        free_field (f.p);
        free_field (f.u);

        MPI_Finalize ();
    }

    return 0;
}
