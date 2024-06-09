#xcog-file-type-1
reset-xcog
#
# We begin by making a cubic spline through node points on the wing profile.
#
make-curve wing
	cubic-spline
		read-node-file wing.spl
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
		reverse-curve-parametrization
		thickness 1
		velocity-threshold .15
		set-s-lines 21
		s-stretching
			exponential-stretching
			starting-grid-size 5.e-3
		exit
		compute-mapping
		mapping-plot-mode
			all-grid-lines
		exit
	exit
#
# Let's turn on the plotting of all curves and all grid lines
#
pause
composite-plot-mode
	curves
	all-grid-lines
exit
# 
# Next we make a spline through an imaginary feature in the flowfield 
#
pause
make-curve feature
	cubic-spline
		enter-new-nodes 3
			0.6
			-0.05
			.7
			.15
			0.8
			0.4
		change-node 1
			.56
			-5.000000e-02
	exit
#
# We make a grid around the feature by the hyperbolic mapping again. But this time,
# we grow the grid in both directions from the curve.
#
pause
make-mapping refine
	hyperbolic-mapping feature
		thickness .2
		both-sided-grid
		compute-mapping
		mapping-plot-mode
			all-grid-lines
		exit
		set-s-lines 11
		compute-mapping
		curvature-coefficient .1
		r-stretching
			exponential-stretching
			starting-grid-size 5.e-3
		exit
		thickness .15
#
# The curve can be edited from within the hyperbolic mapping
#
pause
		change-curve
		change-node 2
			.65
			1.500000e-01
		exit
#
# Next, we compute the mapping
#
pause
		compute-mapping
#
# To make the side close to the wing surface follow the wing exactly, we 
# project the grid onto the wing.
#
pause
		project-side wing
			project-low-r
	exit
#
# In order to make an overlapping grid, we need to assign curve labels
#
pause
curve-label wing-grid
	low-s 1
	high-s 2
exit
curve-label refine
	low-r 1
exit
list-of-components
	show
exit
#
# We are now ready to make the overlapping grid
#
pause
compute-overlap
#
# Cleanup
#
pause
reset-xcog
#
# This completes the hyperbolic-wing demo. Thank's for your attention!
#
pause
