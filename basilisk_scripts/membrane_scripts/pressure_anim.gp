#set yrange [-2:2];
DIRNAME=ARG1
do for [ii=0:250] {
        print ii
        plot 'test_outputs/p_'.ii.'.txt' with lines smooth unique
        pause 0.1
}
pause -1
~                                                                                
~               
