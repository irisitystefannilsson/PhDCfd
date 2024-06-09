#xcog-file-type-1
reset-xcog
overlap-parameters
	interpolation-width 3
	exit
make-curve cylinder
	circular-arc
		radius .4
		end-angle 360
	exit
make-mapping square-grid
	cartesian-mapping
		xmin -1.0
		ymin -1.0
		set-r-lines 101
		set-s-lines 101
	exit
make-mapping cylinder-grid
	normal-curve-mapping cylinder
		reverse-curve-parametrization
		constant-width .5
		set-r-lines 201
		set-s-lines 51
	exit
curve-label cylinder-grid
	low-s 1
	high-s 0
	low-r 3
	high-r 3
exit
curve-label square-grid
	low-r 2
	high-r 2
	low-s 2
	high-s 2
exit
boundary-condition square-grid
	low-r 1
	high-r 1
	low-s 1
	high-s 1
exit
boundary-condition cylinder-grid
	low-s 1
exit
compute-overlap
save-overlapping-grid
#	hdf-format cGrid.hdf
	ascii-format cGrid.acg
	yes


