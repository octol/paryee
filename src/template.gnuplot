set terminal png truecolor size 350,300 font Sans 9
set out "file.png"
set contour
set hidden3d
set cbrange [-0.3:0.8]
set zrange [-0.5:1]
splot "file.tsv" with pm3d  t "file"
