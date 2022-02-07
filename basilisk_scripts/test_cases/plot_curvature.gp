set grid
set key top right
set xlabel 'y'
set ylabel 'kappa' 
# f(x)=(x > 0. ? 100.*x : 1e1000)
# set yrange [:100]
plot '../raw_data/curvature_flat.txt' u 2:3 w lines smooth unique t 'Circle curvature', \
    '../raw_data/curvature_curved.txt' u 2:3 w lines smooth unique t 'Actual curvature', \
    '../raw_data/curvature_adjusted.txt' u 2:3 w lines smooth unique t 'Adjusted curvature'

set term png
set output "curvature_mag_0.25_height_adjustment.png"
replot