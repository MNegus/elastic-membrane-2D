# set yrange [-2:2];
do for [ii=0:1000] {
        print ii
        plot 'outputs/w_'.ii.'.txt' using 1:2 with lines smooth unique
        pause 0.001
}
pause -1
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                  
