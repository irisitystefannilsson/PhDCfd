#xcog-file-type-1
reset-xcog
make-mapping simple-grid
	cartesian-mapping
		xmin -1.75
		xmax 10.
		ymin -1.75
		ymax 1.75
	exit
transform-mapping
	simple-grid
	rotation-angle
	0
	exit
boundary-condition simple-grid
	low-r 50
	high-r 50
	low-s 50
	high-s 50
exit
make-curve circle1
	circular-arc
		start-angle 0.
		end-angle 360.
		origin
0.
0.
		radius .55
		show
	exit
make-curve circle2
	circular-arc
		start-angle 0.
		end-angle 360.
		radius .25
		origin
0.
0.
	exit
grid-lines simple-grid
		set-s-lines 41
		set-r-lines 166
exit
make-mapping circle
	linear-interpolation-mapping
		circle2
		circle1
		set-s-lines 66
		set-r-lines 36
		s-stretching
		hyperbolic-tangent-stretching
		start-grid-size
		2.5e-2
		end-grid-size
		2.5e-1
		show
	exit
exit
boundary-condition circle
	low-r 3
	high-r 3
	low-s 50
	high-s 0
exit
curve-label circle
	low-s 2
	high-s 0
	low-r 0
	high-r 0
exit
curve-label simple-grid
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
#       hdf-format square.hdf



