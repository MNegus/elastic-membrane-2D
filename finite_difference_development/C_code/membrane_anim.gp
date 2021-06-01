set yrange [-0.05:0.15];
DIRNAME=ARG1
do for [ii=0:250] {
        print ii
       	plot DIRNAME.'_outputs/w_'.ii.'.txt' with lines smooth unique 
        pause 0.01
}
pause -1
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                  
