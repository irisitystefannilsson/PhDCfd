#xcog-file-type-1
reset-xcog
make-mapping square-grid
	cartesian-mapping
		xmin -1.0
                xmax 3.0
		ymin -1.0
                ymax 0.0
		set-r-lines 421
		set-s-lines 121
	exit
curve-label square-grid
	low-r -1
	high-r 1
	low-s 1
	high-s 1
exit
boundary-condition square-grid
	low-r 1
	high-r 21
	low-s 1
	high-s 1 
exit
make-mapping square-grid_2
	cartesian-mapping
		xmin -2.0
                xmax 1.0
		ymin -0.5
                ymax 0.0
		set-r-lines 331
		set-s-lines 91
	exit
curve-label square-grid_2
	low-r 1
#        high-r 1
	low-s -1
	high-s 1
exit
boundary-condition square-grid_2
	low-r 1
	high-r 0
	low-s 1
	high-s 1 
exit
overlap-parameters
		ghostpoints 1
		interpolation-width 3
		normal-width 3
	exit
compute-overlap
save-overlapping-grid
#	hdf-format lid_221x421.hdf
	ascii-format step_1xs.acg
	yes
