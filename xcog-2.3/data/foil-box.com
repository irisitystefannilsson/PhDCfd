#xcog-file-type-1
reset-xcog
#
# We start by defining the mapping around a naca0012 airfoil
#
make-mapping foil
	theodorsen-garrick-mapping naca0012.spl
pause
#
# Change the radius of the trailing edge
#
		trailing-edge-radius 1.e-4
pause
#
# We introduce an exponential stretching function to better resolve
# the expected boundary layer close to the wing
#
		s-stretching
			exponential-stretching
		exit
pause
#
# We introduce a layer stretching to attract grid points to the sharp trailing edge.
#
		r-stretching
			layer-stretching
		exit
pause
#
# Change the number of grid points
#
		set-s-lines 15
		set-r-lines 75
	exit
pause
#
# We proceed by making a Cartesian background grid
#
make-mapping box
	cartesian-mapping
		xmin -0.4
		xmax 1.4
		ymin -0.7
		ymax 0.5
		set-r-lines 40
		set-s-lines 25
	exit
pause
#
# Identify the boundary of the computational domain. The domain is in this case
# double connected, so we assign different curve-labels to the two parts of the 
# boundary.
#
curve-label box
	low-r 1
	high-r 1
	high-s 1
	low-s 1
exit
curve-label foil
	low-s 2
exit
pause
#
# Assign boundary condition values to the sides that are aligned with the 
# boundary of the computational domain. This information is not used by the overlap
# algorithm, but indicates the type of boundary-condition to apply in a PDE solver.
#
boundary-condition box
	low-r 2
	high-r 3
	high-s 4
	low-s 4
exit
boundary-condition foil
	low-s 1
exit
pause
# 
# Give the `foil' grid the highest priority in the overlapping grid
#
list-of-components
	move-first foil
exit
pause
#
# Compute the overlapping grid
#
compute-overlap
pause
#
# We change the angle of attack of the wing by rotating that mapping
#
transform-mapping foil
	rotation-angle -10
exit
pause
#
# Compute a new overlapping grid.
#
compute-overlap
pause
#
# Save the state of xcog.
#
save-state foil-box.xcog
pause
#
# Cleanup
#
reset-xcog
pause
#
# End of the foil-box demo. Thank's for your attention!
#
pause

