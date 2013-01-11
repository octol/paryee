README - ParYee
===============

Basic parallel implementation of the finite-difference time-domain (FDTD)
method applied to the 2D wave equation on system form. 

        p_t = a u_x,
        u_t = b p_x,
        v_t = b p_y.

This formulation is sometimes used in acoustics, where p is the pressure and
(u,v) are the x- and y-components of the velocity field.

A number of different implementations are included, of varying complexity.

- yee_mpi  - Distributed memory parallelization using MPI.
- yee_mpi2 - Distributed memory parallelization using MPI (Non-blocking).
- yee_omp  - Shared memory parallelization using OpenMP.
- yee_pthr - Shared memory parallelization using POSIX threads.
- yee      - Single thread reference implementation to compare against.

Visualizations
--------------

![Initial field](https://raw.github.com/octol/paryee/master/figures/yee0.png)
![Field afterwards](https://raw.github.com/octol/paryee/master/figures/yee.png)

Performance tests
-----------------

### Intel Core i7 920, 4 cores (8 threads) @ 2.80 GHz
![Performance tests](https://raw.github.com/octol/paryee/master/tests_swiftsure/tests_perf.png)
![Scaling tests](https://raw.github.com/octol/paryee/master/tests_swiftsure/tests_scaling.png)

### 2 x SGI Origin, MIPS R14000, 4+2 CPUs @ 600 MHz
![Performance tests](https://raw.github.com/octol/paryee/master/tests_asuka/tests_perf.png)
![Scaling tests](https://raw.github.com/octol/paryee/master/tests_asuka/tests_scaling.png)

### Intel Xeon E5310, 4+4 cores  @ 1.60 GHz
![Performance tests](https://raw.github.com/octol/paryee/master/tests_europa/tests_perf.png)
![Scaling tests](https://raw.github.com/octol/paryee/master/tests_europa/tests_scaling.png)


Grid and simulation information
-------------------------------

The grid is discretized in a staggered fashion with p placed in the middle of
the cells, and u and v placed on the middle of the sides.

The number of cells in each coordinate direction is denoted by nx and ny. We
note that we have nx+1 number of u points along the x-axis, and ny+1 number of
v points along the y-axis. 

The field points are grouped such that each cell consists of a single p, u
and v value, where the u and v value is located just before the p points along
its respective axis. 

Thus for the cell indexed by (i,j), we have the coordinates

        p: x = (i + 1/2)*dx, y = (j + 1/2)*dy,
        u: x = i*dx,         y = (j + 1/2)*dy,
        v: x = (i + 1/2)*dx, y = j*dy.

The last u points on the x-axis and the last v points on
y-axis are treated as auxiliary points. This does not present a practical
problem as these are boundary cells anyway, which are not updated.

The outer boundary is set to u=0 at x=0 and x=L, and v=0 at y=0 and y=L. Thus
we only update the inner points. 

The time axis is also staggered, with (u,v) defined at t=0 and p at t=-dt/2 as
initial conditions. 

### Example: nx=4, ny=4

          v   v   v   v     
        u p u p u p u p u   
          v   v   v   v
        u p u p u p u p u
          v   v   v   v 
        u p u p u p u p u
          v   v   v   v 
        u p u p u p u p u
          v   v   v   v

Partition of grid (shared memory)
---------------------------------

The grid is split along the x-axis so that the different threads stay close in
memory to each other.

### Example: nx=4, ny=4, 2 partitions

          v   v      v   v     
        u p u p    u p u p u   
          v   v      v   v
        u p u p    u p u p u
          v   v      v   v
        u p u p    u p u p u
          v   v      v   v
        u p u p    u p u p u
          v   v      v   v

### Example: nx=4, ny=4, 4 partitions

          v      v      v      v     
        u p    u p    u p    u p u   
          v      v      v      v
        u p    u p    u p    u p u
          v      v      v      v
        u p    u p    u p    u p u
          v      v      v      v
        u p    u p    u p    u p u
          v      v      v      v


Partition of grid (distributed memory)
--------------------------------------

Note that the points surrounded by a parenthesis in the examples below are
ghost points that needs to be sent between the computational nodes.

The grid is split along the y-axis to make it at bit more convenient to send
data as single segments.

### Example: nx=4, ny=4, 2 partitions

          v   v   v   v     
        u p u p u p u p u   
          v   v   v   v
        u p u p u p u p u
          v   v   v   v
         (p) (p) (p) (p)
        
         (v) (v) (v) (v)
        u p u p u p u p u
          v   v   v   v
        u p u p u p u p u
          v   v   v   v

### Example: nx=4, ny=4, 4 partitions

          v   v   v   v     
        u p u p u p u p u  
          v   v   v   v
         (p) (p) (p) (p)

         (v) (v) (v) (v)
        u p u p u p u p u
          v   v   v   v
         (p) (p) (p) (p)

         (v) (v) (v) (v)
        u p u p u p u p u
          v   v   v   v
         (p) (p) (p) (p)

         (v) (v) (v) (v)
        u p u p u p u p u
          v   v   v   v

Author
------

Jon Haggblad <jon@haeggblad.com>


