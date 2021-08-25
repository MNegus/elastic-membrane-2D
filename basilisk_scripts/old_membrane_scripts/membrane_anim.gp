#set yrange [-0.05:0.15];
do for [ii=0:4000] {
        print ii
       	plot '../membrane_outputs/w_'.ii.'.txt' with lines smooth unique 
        pause 0.01
}
pause -1       
