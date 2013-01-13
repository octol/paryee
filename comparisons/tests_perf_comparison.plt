set terminal png truecolor size 400,300 font "Helvetica" 10
set out "tests_perf_comparison.png"
set key top left
set xlabel "Cells per axis"
set ylabel "Time [s]"
set title "Comparison"
plot \
    "tests_perf_comparison.tsv" using 1:2 title "i7 920" with lines lw 1, \
    "tests_perf_comparison.tsv" using 1:3 title "MIPS R14000" with lines lw 1, \
    "tests_perf_comparison.tsv" using 1:4 title "Xeon E5310" with lines lw 1

set logscale xy
set out "tests_perf_comparison2.png"
replot
