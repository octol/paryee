set terminal png truecolor size 400,300 font "Helvetica" 10
set out "test_scaling.png"
set key top right
set xlabel "Threads"
set ylabel "Time [s]"
set title "1000 x 1000 cells"
plot \
    "test_scaling.tsv" using 1:2 title "OpenMP" with linespoints lw 2, \
    "test_scaling.tsv" using 1:3 title "Pthread" with linespoints lw 2, \
    "test_scaling.tsv" using 1:4 title "MPI" with linespoints lw 2, \
    "test_scaling.tsv" using 1:5 title "MPI (Non-blocking)" with linespoints lw 2

set logscale xy
set out "test_scaling_log.png"
replot
