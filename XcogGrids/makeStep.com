#xcog-file-type-1
reset-xcog
make-mapping square-grid
	cartesian-mapping
		xmin -1.0
                xmax 3.0
		ymin -1.0
                ymax 0.0
		set-r-lines 421
		set-s-lines 141
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
		set-r-lines 301
		set-s-lines 71
	exit
curve-label square-grid_2
	low-r 1
        high-r 0
	low-s -1
	high-s 1
exit
boundary-condition square-grid_2
	low-r 1
	high-r 0
	low-s 1
	high-s 1 
exit
make-curve lower-corner
	smooth-polygon
		enter-new-corners 3
			-1.1
			-0.5
			-1.005
			-0.501
			-1
			-0.6
		corner-sharpness 80
	exit
#pause
#
# The domain close to the lower fillet is gridded by taking normals out from the
# lower fillet curve.
#
make-mapping lower-corner-grid
	normal-curve-mapping lower-corner
		constant-width .05
		r-stretching
			layer-stretching
			change-width 3
			0.5
		exit
	exit
grid-lines lower-corner-grid
	set-r-lines 60
	set-s-lines 16
exit
curve-label lower-corner-grid
	low-s 1
exit
boundary-condition lower-corner-grid
	low-s 1
exit
overlap-parameters
		ghostpoints 1
		interpolation-width 3
	        normal-width 3
	        #interpolation-type implicit
	exit
compute-overlap
save-overlapping-grid
#	hdf-format lid_221x421.hdf
	ascii-format step_1xs.acg
	yes
