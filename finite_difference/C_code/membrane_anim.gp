#set yrange [-2:2];
DIRNAME=ARG1
do for [ii=0:1000] {
        print ii
       	plot DIRNAME.'_outputs/w_'.ii.'.txt' with lines smooth unique 
        pause 0.001
}
pause -1
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                  
