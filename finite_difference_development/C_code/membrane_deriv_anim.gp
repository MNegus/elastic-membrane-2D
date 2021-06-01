# set yrange [-2:2];
do for [ii=0:2500] {
        print ii
        plot 'mitchell_outputs/w_deriv_'.ii.'.txt' using 1:2 with lines smooth unique
        pause 0.01
}
pause -1
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                  
