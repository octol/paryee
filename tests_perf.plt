set terminal png truecolor
set out "tests_perf.png"
plot "tests_perf.tsv" with lines
