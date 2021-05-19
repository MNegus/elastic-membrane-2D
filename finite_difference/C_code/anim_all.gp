#set yrange[-1.8:1.8];
do for [ii=0:1000] {
        print ii
       	# plot 'explicit_outputs/w_'.ii.'.txt' with lines smooth unique, \
        #      'implicit_outputs/w_'.ii.'.txt' with lines smooth unique, \
	#      'exact_outputs/w_'.ii.'.txt' with lines smooth unique
        plot 'mitchell_outputs/w_'.ii.'.txt' with lines smooth unique, \
	     'exact_outputs/w_'.ii.'.txt' with lines smooth unique
        pause 0.0001
}
pause -1
