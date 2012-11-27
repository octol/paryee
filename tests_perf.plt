set terminal png truecolor font "Helvetica"
set out "tests_perf.png"
set key bottom right
set logscale xy
set xlabel "Cells"
set ylabel "Time [s]"
plot \
    "tests_perf.tsv" using 1:3 title "Ref" with linespoints, \
    "tests_perf.tsv" using 1:2 title "Single thread" with linespoints, \
    "tests_perf.tsv" using 1:4 title "OpenMP" with linespoints, \
    "tests_perf.tsv" using 1:6 title "Pthread" with linespoints, \
    "tests_perf.tsv" using 1:7 title "MPI" with linespoints

    #"tests_perf.tsv" using 1:5 title "Pthread (naive)" with linespoints, \
