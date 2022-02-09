set logscale
set grid
set key top right
set xlabel 'Diameter (grid points)'
set ylabel 'Relative curvature error / percentage' 
f(x)=(x > 0. ? 100.*x : 1e1000)
set yrange [:100]
plot 2./(x*x) t '2/x^{2}', '../raw_data/log' u 1:4 w lp t 'Max', '' u 1:3 w lp t 'RMS', \
  'popinet.csv' u ($1*2):2 w lp t 'Popinet (2009)', \
  '../raw_data/log' u 1:(f($5)) w lp t 'HF', '' u 1:(f($6)) w lp t 'HF fit', \
  '' u 1:(f($7)) w lp t 'Average', '' u 1:(f($8)) w lp t 'Centroids'
set term png
set output "curvature_errorplot.png"
replot