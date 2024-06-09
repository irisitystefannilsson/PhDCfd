#xcog-file-type-1
reset-xcog
#
# The first curve describes the paper surface. We use a smoothed polygon curve 
# in order to achieve a non-uniform distribution of grid points along the 
# straight line.
#
make-curve paper
	smooth-polygon
		enter-new-corners 4
			-12
			0
			-1.5
			0
			0
			0
			4.3
			0
	exit
pause
#
# The second curve describes the free surface. We also describe this curve by a 
# smoothed polygon, but it could equally well have been a cubic spline.
#
#
make-curve free-surface
	smooth-polygon
		enter-new-corners 7
			-6.98192
			6.974
			-4.923
			4
			-4.8
			3.5
			-4.8
			3
			-5.3
			2.25
			-6
			2.0
			-12
			2
pause
#
# Make the corners smoother than the default value.
#
		corner-sharpness 3
pause
	exit
#
# The third curve describes the front face of the paper knife. We use a smoothed
# polygon to properly resolve the domain around the convex corner and also attract 
# grid points to this area.
#
make-curve paper-knife
	smooth-polygon
		enter-new-corners 4
			4.3
			0.15
			0
			0.15
			-1
			1.5
			-5.5
			8
pause
#
# Make the corners sharper than the default value.
#
		corner-sharpness 7
pause
	exit
#
# The domain close to the paper surface is discretized by taking normals of 
# variable width out from the paper surface. Because this surface is treated as
# a no-slip surface in the solver, we will use an exponential stretching function
# to attract grid points close to the surface.
#
make-mapping paper-grid
	normal-curve-mapping paper
		variable-width
			1.2
			1.2
			1.2
			.1
			.1
			.1
		width-change-sharpness 50
		set-r-lines 110
		set-s-lines 20
		r-stretching
			layer-stretching
			change-strength 2
			0
			change-strength 3
			.4
			change-width 3
			0.1
		exit
		s-stretching
			exponential-stretching
			starting-grid-size 0.01
		exit
	exit
pause
#
# The domain close to the free surface is discretized by taking normals of 
# constant width out from the free surface. Because this surface is treated as
# a slip surface in the solver, no boundary layer can be expected here.
#
make-mapping free-surface-grid
	normal-curve-mapping free-surface
		constant-width 1.2
		set-r-lines 60
		set-s-lines 8
	exit
pause
#
# The domain close to the paper knife is discretized by taking normals of 
# variable width out from the knife. Because this surface is treated as
# a no-slip surface in the solver, we will use an exponential stretching 
# function to attract grid points close to the boundary.
#
make-mapping paper-knife-grid
	normal-curve-mapping paper-knife
		variable-width
			.1
			.1
			.1
			1.2
			1.2
			1.2
		width-change-sharpness 100
		set-r-lines 100
		set-s-lines 20
pause
#
# It is important to set the exponential stretching after the number of gridlines 
# has been set, because the exponential stretching parameter is the relative increment 
# in grid step, so the coefficients in the stretching function will depend on the number 
# of grid points in the normal direction.
#
		r-stretching
			layer-stretching
			change-strength 2
			.5
			change-width 2
			0.2
			change-strength 3
			0
		exit
		s-stretching
			exponential-stretching
			starting-grid-size 0.01
		exit
pause
		mapping-plot-mode
			grid-tickmarks
		exit
	exit
#
# We fill the hole in the middle of the domain with a Cartesian grid.
#
make-mapping base-grid
	cartesian-mapping
		xmin -6.0
		xmax -1.5
		ymin 0.5
		ymax 3.8
		set-r-lines 30
		set-s-lines 25
	exit
pause
#
# Label the curves that are aligned with the boundary of the computational domain.
# This is a simply connected domain, so we assign the same value of the curve-label
# to all those sides.
#
curve-label paper-grid
	low-r 1
	high-r 1
	low-s 1
exit
curve-label free-surface-grid
	low-r 1
	high-r 1
	low-s 1
exit
curve-label paper-knife-grid
	low-r 1
	high-r 1
	low-s 1
exit
pause
#
# To make the grid useful for a PDE solver, we must assign boundary condition labels
# to all sides that are aligned with the boundary of the computational domain.
#
boundary-condition paper-grid
	low-r 2
	high-r 6
	low-s 1
exit
boundary-condition free-surface-grid
	low-r 4
	high-r 2
	low-s 3
exit
boundary-condition paper-knife-grid
	low-r 6
	high-r 4
	low-s 5
exit
pause
#
# Change the priority of the component grids.
#
list-of-components
#
# Move the base-grid last.
#
	move-last base-grid
pause
exit
#
# Compute the overlapping grid
#
compute-overlap
pause
#
# Save the state of xcog, i.e. all curves, mappings and the value of all overlap
# parameters. After the program has been restarted or after a `reset-xcog', these
# quantities can be recovered by the command `read-state'. To save disc space, the
# overlapping grid is not saved and has to be recomputed.
#
save-state paper-knife.xcog
pause
#
# Cleanup
#
reset-xcog
#
# End of the paper coating example. Thank's for your attention!
#
pause
