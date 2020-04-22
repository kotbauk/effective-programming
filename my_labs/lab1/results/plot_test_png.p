nx=100
ny=100
filename="results_matrix.bin"
set terminal png size 700,600
set output filename.".png"
set xrange[-1:nx]
set yrange[-1:ny]
set palette gray
set title filename
plot filename binary array=(ny,nx) format="%lf" with image
