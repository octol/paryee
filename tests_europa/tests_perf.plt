set terminal png truecolor font "Helvetica"
set out "tests_perf.png"
set key top left
#set logscale xy
set xlabel "Cells"
set ylabel "Time [s]"
set title "4 threads"
plot \
    "tests_perf.tsv" using 1:3 title "OpenMP" with lines, \
    "tests_perf.tsv" using 1:4 title "Pthread" with lines, \
    "tests_perf.tsv" using 1:5 title "MPI" with lines
