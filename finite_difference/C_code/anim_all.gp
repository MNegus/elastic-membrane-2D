do for [ii=0:10000] {
        print ii
       	# plot 'explicit_outputs/w_'.ii.'.txt' with lines smooth unique, \
        #      'implicit_outputs/w_'.ii.'.txt' with lines smooth unique, \
	#      'exact_outputs/w_'.ii.'.txt' with lines smooth unique
        plot 'implicit_outputs/w_'.ii.'.txt' with lines smooth unique, \
	     'exact_outputs/w_'.ii.'.txt' with lines smooth unique
        pause 0.0001
}
pause -1
