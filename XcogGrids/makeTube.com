#xcog-file-type-1
reset-xcog
make-mapping square-grid
	cartesian-mapping
		xmin -1.0
                xmax 12.0
		ymin -1.0
		set-r-lines 321
		set-s-lines 81
	exit
curve-label square-grid
	low-r 1
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
overlap-parameters
		ghostpoints 1
		interpolation-width 3
		normal-width 3
	exit
compute-overlap
save-overlapping-grid
#	hdf-format sharp_cGrid_2.hdf
	ascii-format tube_81x321.acg
	yes
