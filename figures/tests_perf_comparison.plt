#set terminal png truecolor size 400,300 font "Helvetica" 10
set terminal png truecolor size 800,450 font "Helvetica" 10
set out "tests_perf_comparison.png"
set key top left
set xlabel "Cells per axis"
set ylabel "Time [s]"
set title "Comparison, MPI (Non-blocking)"
set xrange [500:1000]
set yrange [0:6]
plot \
    "tests_perf_comparison2.tsv" using 1:2 title "i7 930 (4 cores)" with lines lw 1, \
    "tests_perf_comparison2.tsv" using 1:3 title "MIPS R14000 (10 cpus)" with lines lw 1, \
    "tests_perf_comparison2.tsv" using 1:4 title "Xeon E5310 (sun) (2x4 cores)" with lines lw 1, \
    "tests_perf_comparison2.tsv" using 1:5 title "Xeon E5310 (gcc) (2x4 cores)" with lines lw 1

set logscale xy
set xrange [100:1000]
set yrange [0.001:10]
set out "tests_perf_comparison2.png"
replot
