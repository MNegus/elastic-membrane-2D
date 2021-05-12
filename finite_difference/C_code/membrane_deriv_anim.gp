# set yrange [-2:2];
do for [ii=0:1000] {
        print ii
        plot 'outputs/w_deriv_'.ii.'.txt' using 1:2 with lines smooth unique
        pause 0.1
}
pause -1
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                  
