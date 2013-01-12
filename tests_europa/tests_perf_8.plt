set terminal png truecolor size 400,300 font "Helvetica" 10
set out "tests_perf_8.png"
set key top left
set xlabel "Cells per axis"
set ylabel "Time [s]"
set title "8 threads"
plot \
    "tests_perf_8.tsv" using 1:2 title "OpenMP" with linespoints lw 2, \
    "tests_perf_8.tsv" using 1:3 title "Pthread" with linespoints lw 2, \
    "tests_perf_8.tsv" using 1:4 title "MPI" with linespoints lw 2, \
    "tests_perf_8.tsv" using 1:5 title "MPI (Non-blocking)" with linespoints lw 2

set logscale xy
set out "tests_perf2_8.png"
replot
