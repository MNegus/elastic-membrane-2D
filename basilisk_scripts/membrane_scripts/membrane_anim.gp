#set yrange [-0.05:0.15];
do for [ii=0:250] {
        print ii
       	plot 'test_outputs/w_'.ii.'.txt' with lines smooth unique 
        pause 0.1
}
pause -1
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                                                                                                                                 
~                  
