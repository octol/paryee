set terminal png truecolor size 700,350 font "Helvetica" 10
set out "test_perf_4.png"
set key top right out
set xlabel "Cells per axis"
set ylabel "Time [s]"
set title "4 threads"
set xrange [500:1000]
plot \
    "test_perf_4.tsv" using 1:2 title "OpenMP (naive)" with lines lw 1.5, \
    "test_perf_4.tsv" using 1:3 title "OpenMP (blocked)" with lines lw 1.5, \
    "test_perf_4.tsv" using 1:4 title "OpenMP (stride: 1)" with lines lw 1.5, \
    "test_perf_4.tsv" using 1:5 title "Pthread (blocked)" with lines lw 1.5, \
    "test_perf_4.tsv" using 1:6 title "Pthread (stride: 1)" with lines lw 1.5, \
    "test_perf_4.tsv" using 1:7 title "MPI (blocking)" with lines lw 1.5, \
    "test_perf_4.tsv" using 1:8 title "MPI (non-blocking)" with lines lw 1.5

set logscale xy
set grid
#set xrange [100:1000]
set out "test_perf_log_4.png"
replot
