#set yrange [-2:2];
DIRNAME=ARG1
do for [ii=0:2500] {
        print ii
       	plot DIRNAME.'_outputs/w_'.ii.'.txt' with lines smooth unique 
        pause 0.0001
}
pause -1
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                  
