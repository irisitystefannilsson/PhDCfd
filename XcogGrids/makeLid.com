#xcog-file-type-1
reset-xcog
make-mapping square-grid
	cartesian-mapping
		xmin -1.0
                xmax 1.0
		ymin -1.0
                ymax 3.0
		set-r-lines 221
		set-s-lines 421
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
overlap-parameters
		ghostpoints 1
		interpolation-width 3
		normal-width 3
	exit
compute-overlap
save-overlapping-grid
	hdf-format lid_221x421.hdf
#	ascii-format lid_221x421.acg
	yes
