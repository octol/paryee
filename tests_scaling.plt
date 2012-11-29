set terminal png truecolor font "Helvetica"
set out "tests_scaling.png"
set key bottom right
#set logscale xy
set xlabel "Nodes"
set ylabel "Time [s]"
plot \
    "tests_scaling.tsv" using 1:3 title "Ref" with linespoints, \
    "tests_scaling.tsv" using 1:2 title "Single thread" with linespoints, \
    "tests_scaling.tsv" using 1:4 title "OpenMP" with linespoints, \
    "tests_scaling.tsv" using 1:5 title "Pthread" with linespoints, \
    "tests_scaling.tsv" using 1:6 title "MPI" with linespoints
