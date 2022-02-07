set grid
set key top right
set xlabel 'y'
set ylabel 'heights' 
# f(x)=(x > 0. ? 100.*x : 1e1000)
# set yrange [:100]
plot '../raw_data/heights_y_mag_0.txt' u 2:3 w lines smooth unique t 'Flat x height', \
    '../raw_data/heights_y_mag_0.25.txt' u 2:3 w lines smooth unique t 'Curved x height', \
# plot '../raw_data/heights_y_mag_0.txt' u 2:3 w lines smooth unique t 'Flat y height', \
    # '../raw_data/heights_y_mag_0.25.txt' u 2:3 w lines smooth unique t 'Curved y height'

set term png
set output "heights_plot.png"
replot