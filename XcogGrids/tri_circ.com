#xcog-file-type-1
reset-xcog
make-mapping simple-grid
	cartesian-mapping
		xmin -2.75
		xmax 15.
		ymin -2.75
		ymax 2.75
	exit
transform-mapping
	simple-grid
	rotation-angle
	0
	exit
grid-lines simple-grid
		set-s-lines 261
		set-r-lines 661
exit
boundary-condition simple-grid
	low-r 1
	high-r 21
	low-s 5
	high-s 5
exit
make-curve circle1
	circular-arc
		start-angle 0.
		end-angle 360.
		origin
2.
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
2.
0.
	exit
make-mapping circle
	linear-interpolation-mapping
		circle2
		circle1
		set-s-lines 131
		set-r-lines 241
		s-stretching
		hyperbolic-tangent-stretching
		start-grid-size
		2.5e-2
		end-grid-size
		8.5e-1
		show
	exit
exit
boundary-condition circle
	low-r 3
	high-r 3
	low-s 1
	high-s 0
exit
curve-label circle
	low-s 2
	high-s 0
	low-r 0
	high-r 0
exit
copy-mapping circle 
	circle_2
exit
transform-mapping circle_2
	horizontal-translation -1
	vertical-translation -1
	
exit
copy-mapping circle 
	circle_3
exit
transform-mapping circle_3
	horizontal-translation -1
	vertical-translation 1
	
exit
curve-label circle_2
	low-s 3
	high-s 0
	low-r 0
	high-r 0
exit
curve-label circle_3
	low-s 4
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
#save-overlapping-grid
#       hdf-format square.hdf
#        ascii-format circle_tanh.acg
#        yes



