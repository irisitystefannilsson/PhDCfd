#xcog-file-type-1
reset-xcog
overlap-parameters
	normal-width 3
	interpolation-width 3
	ghostpoints 1
	exit
make-mapping square-grid
	cartesian-mapping
		xmin -1.0
		ymin -1.0
		set-r-lines 41
		set-s-lines 41
	exit
curve-label square-grid
	low-r 1
	high-r 1
	low-s 1
	high-s 1
exit
boundary-condition square-grid
	low-r 1
	high-r 1
	low-s 1
	high-s 1
exit
compute-overlap
save-overlapping-grid
#	hdf-format square.hdf
	ascii-format square3.acg
	yes
quit


