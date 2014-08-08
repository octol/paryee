set terminal png truecolor size 400,300 font "Helvetica" 10
set out "tests_scaling.png"
set key top right
set xlabel "Threads"
set ylabel "Time [s]"
set title "1000 x 1000 cells"
set yrange [*:35]
plot \
    "tests_scaling.tsv" using 1:3 title "Pthread" with linespoints lw 2, \
    "tests_scaling.tsv" using 1:4 title "MPI" with linespoints lw 2, \
    "tests_scaling.tsv" using 1:5 title "MPI (Non-blocking)" with linespoints lw 2

set logscale xy
set out "tests_scaling2.png"
replot
