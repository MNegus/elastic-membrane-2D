set xrange [0:1.5];
set yrange [-2:15];
do for [ii=0:4000] {
        print ii
       	plot '../membrane_outputs/p_'.ii.'.txt' with lines smooth unique 
        pause 0.01
}
pause -1       