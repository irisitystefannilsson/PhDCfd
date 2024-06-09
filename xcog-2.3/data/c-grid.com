#xcog-file-type-1
reset-xcog
#
# We begin by making a cubic spline through node points on the wing profile.
#
make-curve wing
	cubic-spline
		read-node-file naca0012.spl
	exit
#
# We rotate the wing 10 degrees clockwise
#
pause
transform-curve wing
	rotation-angle -10
exit
#
# We make a grid outside of the wing by using the hyperbolic grid generator
#
pause
make-mapping wing-grid
	hyperbolic-mapping wing
		r-stretching
			layer-stretching
			change-width 18
			.1
			change-strength 18
			.2
		exit
		set-r-lines 90
		thickness 1
		velocity-threshold .15
		set-s-lines 21
		s-stretching
			exponential-stretching
			starting-grid-size 5.e-3
		exit
#
# Turn on the C-grid topology and specify the end-location of the branch cut
#
pause
		c-grid
		end-location
			1.000000e+00
			-.2
#
# Also specify the number of grid point along the branch cut.
#
pause
		branch-points 20
		compute-mapping
	exit
#
# We enclose the wing in a rectangular background grid
#
pause
make-mapping background
	cartesian-mapping
		xmin -1
		xmax 2.5
		ymin -1.5
		ymax 1
		set-r-lines 100
		set-s-lines 70
	exit
#
# In order to make an overlapping grid, we need to assign curve labels
#
pause
curve-label background
	low-r 1
	high-r 1
	low-s 1
	high-s 1
exit
curve-label wing-grid
	low-s 2
exit
show-curves-mappings
list-of-components
	show
	move-first wing-grid
exit
overlap-parameters
		show-physical-boundaries yes
	exit
#
# We are now ready to make the overlapping grid
#
pause
compute-overlap
#
# This is the hole-cutting boundary. Note that only the physical part of the
# boundary in the C-grid is included.
#
pause
	proceed
#
# Cleanup
#
pause
reset-xcog
#
# This completes the c-grid demo. Thank's for your attention!
#
pause


