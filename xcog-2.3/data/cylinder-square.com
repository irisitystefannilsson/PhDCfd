#xcog-file-type-1
reset-xcog
#
# The first curve is a circular arc with radius 0.4.
#
make-curve cylinder
	circular-arc
		radius .4
		end-angle 360
	exit
pause
#
# The background grid is a Cartesian grid.
#
make-mapping square-grid
	cartesian-mapping
		xmin -1.0
		ymin -1.0
		set-r-lines 25
		set-s-lines 25
	exit
pause
#
# The grid around the circle is constructed by taking normals out from the circle.
#
make-mapping cylinder-grid
	normal-curve-mapping cylinder
		reverse-curve-parametrization
		constant-width .5
	exit
pause
#
# This is a doubly connected domain and it is necessary to identify both boundary 
# curves to the algorithm. We do this by assigning curve-labels to the sides that
# are aligned with the boundary of the domain.
#
curve-label cylinder-grid
	low-s 1
exit
curve-label square-grid
	low-r 2
	high-r 2
	low-s 2
	high-s 2
exit
pause
#
# To make this grid useful for a PDE solver, we also need to assign boundary
# condition values for the sides that are aligned with the boundary.
#
boundary-condition square-grid
	low-r 1
	high-r 1
	low-s 1
	high-s 1
exit
boundary-condition cylinder-grid
	low-s 1
exit
pause
#
# We will now make a composite grid from the overlapping components.
#
compute-overlap
pause
#
# Save the composite grid on an ascii file.
#
save-overlapping-grid
	ascii-format cylinder-square.acg
	yes
pause
#
# Save the state of xcog, i.e. all curves, mappings and the value of all overlap
# parameters. After the program has been restarted or after a `reset-xcog', these
# quantities can be recovered by the command `read-state'. To save disc space, the
# overlapping grid is not saved and has to be recomputed.
#
save-state cylinder-square.xcog
pause
#
# Cleanup
#
reset-xcog
pause
#
# End of the cylinder in a square example. Thank's for your attention! 
#
pause