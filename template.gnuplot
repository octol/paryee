set terminal png truecolor font Helvetica 9
set out "file.png"
set contour
set hidden3d
set cbrange [-0.5:1]
set zrange [-0.5:1]
splot "file.tsv" with pm3d  t "file"
