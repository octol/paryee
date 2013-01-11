set terminal png truecolor size 400,300 font "Helvetica" 10
set out "tests_perf.png"
set key top left
set xlabel "Cells"
set ylabel "Time [s]"
set title "4 threads"
plot \
    "tests_perf.tsv" using 1:3 title "OpenMP" with lines lw 2, \
    "tests_perf.tsv" using 1:4 title "Pthread" with lines lw 2, \
    "tests_perf.tsv" using 1:5 title "MPI" with lines lw 2, \
    "tests_perf.tsv" using 1:6 title "MPI (Non-blocking)" with lines lw 2

set logscale xy
set out "tests_perf2.png"
replot
