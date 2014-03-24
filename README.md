README (paryee)
===============

Toy implementation of the finite-difference time-domain (FDTD) method to
illustrate different parallel implementations. The FDTD is sometimes referred
to as the Yee scheme, and is used to simulate wave propagation.

The system of equations under consideration is:

        p_t = a u_x,
        u_t = b p_x,
        v_t = b p_y.

This formulation is sometimes used in acoustics, where p is the pressure and
(u,v) are the x- and y-components of the velocity field. 

This is a very simple example of how high performance data-parallelism can done
in C for typical numerical codes.  A number of different implementations (of
varying complexity) are included.

- yee_mpi  - Distributed memory parallelization using MPI.
- yee_mpi2 - Distributed memory parallelization using MPI (Non-blocking).
- yee_omp  - Shared memory parallelization using OpenMP.
- yee_pthr - Shared memory parallelization using POSIX threads.
- yee      - Single thread reference implementation to compare against.

Installation
------------

Requires GNU Make to build. CUnit is required to build the unittests. numdiff is required for the integration tests.

Visualizations
--------------

![Initial field](https://raw.github.com/octol/paryee/master/figures/yee0.png)
![Field afterwards](https://raw.github.com/octol/paryee/master/figures/yee.png)

Performance tests
-----------------

### Intel Core i7 930, 4 cores (8 threads) @ 2.80 GHz

    Linux 3.5.0-22-generic #34-Ubuntu SMP Tue Jan 8 21:47:00 UTC 2013 x86_64 x86_64 x86_64 GNU/Linux.
    gcc 4.7.2

![Performance tests](https://raw.github.com/octol/paryee/master/performance-testing/saved/swiftsure/tests_perf_8.png)
![Scaling tests](https://raw.github.com/octol/paryee/master/performance-testing/saved/swiftsure/tests_scaling.png)

### 2 x SGI Origin, MIPS R14000, 4+2 CPUs @ 600 MHz

    IRIX64 6.5 IP35
    MIPSpro Compilers: Version 7.4.4m

![Performance tests](https://raw.github.com/octol/paryee/master/performance-testing/saved/asuka/tests_perf_6.png)
![Scaling tests](https://raw.github.com/octol/paryee/master/performance-testing/saved/asuka/tests_scaling.png)

### Intel Xeon E5310, 4+4 cores  @ 1.60 GHz

    SunOS 5.10 i86pc i386 i86pc
    Sun Studio 12 Update 1

![Performance tests](https://raw.github.com/octol/paryee/master/performance-testing/saved/europa/tests_perf_8.png)
![Scaling tests](https://raw.github.com/octol/paryee/master/performance-testing/saved/europa/tests_scaling.png)

    SunOS 5.10 i86pc i386 i86pc
    gcc 4.3.3

![Performance tests](https://raw.github.com/octol/paryee/master/performance-testing/saved/europa-gcc/tests_perf_8.png)
![Scaling tests](https://raw.github.com/octol/paryee/master/performance-testing/saved/europa-gcc/tests_scaling.png)

### Comparison
![Comparison](https://raw.github.com/octol/paryee/master/figures/tests_perf_comparison.png)


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


