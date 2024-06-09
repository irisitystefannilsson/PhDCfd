#xcog-file-type-1
reset-xcog
#
# This demo illustrates how to incorporate externally defined grids into Xcog.
# We begin by reading in a grid over the coastal line of the mexican gulf. The
# grid is in PLOT3D ascii format. This will take a little time, because the grid
# contains many grid points.
#
make-mapping mexico
	discrete-point-mapping
	plot3d-ascii mexico.plot3d
	exit
pause
#
# We fill the atlantic with a Cartesian grid.
#
make-mapping atlantic
	cartesian-mapping
		xmax 320
		xmin 30
		ymax 267.5
		ymin -28
	exit
pause
#
# Set the number of grid lines.
#
grid-lines mexico
	set-s-lines 7
exit
grid-lines atlantic
	set-r-lines 70
	set-s-lines 75
exit
pause
#
# Make external gaps in the mixed physical / exterior grids in `atlantic'
#
mark-mapping-boundary atlantic
	low-s
		0.0
		0.59
	high-s
		0.0
		0.877
exit
pause
#
# Set the curve label to identify the boundary of the computational domain.
# Since the domain is simply connected, we use the same value for the curve-label
# for all sides of all components that are aligned with the boundary
#
curve-label mexico
	low-r 1
	high-r 1
	low-s 1
exit
curve-label atlantic
	low-s 1
	high-s 1
	high-r 1
exit
pause
#
# Set the boundary condition
#
boundary-condition mexico
	low-s 1
	low-r 1
	high-r 1
exit
boundary-condition atlantic
	low-s 1
	high-s 1
	high-r 1
exit
pause
#
# We give the coastal grid the highest priority in the overlapping grid.
#
pause
list-of-components
	move-first mexico
exit
#
# We proceed by making an overlapping grid for this simplified (Fidel wouldn't 
# like this) version of the Mexican Gulf.
#
compute-overlap
pause
#
# Save the state of xcog, i.e. all curves, mappings and the value of all overlap
# parameters. After the program has been restarted or after a `reset-xcog', these
# quantities can be recovered by the command `read-state'. To save disc space, the
# overlapping grid is not saved and has to be recomputed.
#
save-state mex-gulf.xcog
pause
#
# Cleanup
#
reset-xcog
pause
#
# End of the external grid example. Thank's for your attention!
#
pause

