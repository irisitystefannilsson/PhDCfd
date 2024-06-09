#xcog-file-type-1
reset
#
# The main compartment is discretized with a Cartesian grid.
#
make-mapping main-compartment
	cartesian-mapping
	exit
grid-lines main-compartment
	set-r-lines 50
	set-s-lines 50
exit
pause
#
# The outflow pipe is also discretized with a Cartesian grid.
#
make-mapping outflow-pipe
	cartesian-mapping
		xmin .8
		xmax 1.5
		ymin .375
		ymax .625
	exit
grid-lines outflow-pipe
	set-r-lines 40
	set-s-lines 15
exit
pause
#
# The corners are smoothed out with fillets. We start with the lower corner.
#
make-curve lower-corner
	smooth-polygon
		enter-new-corners 3
			1
			0.175
			1
			0.375
			1.2
			0.375
		corner-sharpness 20
	exit
pause
#
# The domain close to the lower fillet is gridded by taking normals out from the
# lower fillet curve.
#
make-mapping lower-corner-grid
	normal-curve-mapping lower-corner
		constant-width .1
		r-stretching
			layer-stretching
			change-width 2
			0.25
		exit
	exit
grid-lines lower-corner-grid
	set-r-lines 30
	set-s-lines 8
exit
pause
#
# We make the upper corner fillet and grid by copying the lower corner grid and
# rotating and translating that grid.
#
copy-mapping lower-corner-grid
upper-corner-grid
	rotation-angle 90
	horizontal-translation 1.375
	vertical-translation -.375
exit
pause
#
# In applications like the present one, where the component grids stick
# into each other, the boundaries must be treated slightly differently compared
# to the normal case. It is therefore necessary to give Xcog some extra 
# information so that it can handle this case. 
#
# All boundaries that in part are physical boundaries and in part interpolation
# boundaries must be given a NEGATIVE CURVE LABEL, where the absolute value equals 
# the curve label of that part of the physical boundary. In the present case, the 
# domain is simply connected, and we assign a 1 to all physical sides, and a -1 to all
# mixed sides
#
curve-label main-compartment
	low-r 1
	high-s 1
	low-s 1
	high-r -1
exit
curve-label outflow-pipe
	high-r 1
	high-s -1
	low-s -1
exit
curve-label upper-corner-grid
	low-s 1
exit
curve-label lower-corner-grid
	low-s 1
exit
pause
#
# Lets set the boundary condition to make grid useful for the PDE solver.
#
boundary-condition main-compartment
	low-r 2
	high-r 1
	low-s 1
	high-s 1
exit
boundary-condition outflow-pipe
	high-r 2
	low-s 1
	high-s 1
exit
boundary-condition upper-corner-grid
	low-s 1
exit
boundary-condition lower-corner-grid
	low-s 1
exit
pause
#
# We will now make a composite grid from the overlapping components.
#
# We will use the default order of the component grids, i.e. the mappings are inserted
# first in the hierarcy when they are created. Therefore, the most recently made grid
# will be first.
#
compute-overlap
pause
#
# Save all curves and mappings
#
save-state internal-flow.xcog
pause
#
# It is also possible to stick one grid into another without smoothing out the
# corner with a fillet. Because the intersection of the two boundary curves not
# necessarily coincide with grid points on both boundary curves, the corner
# will not have a very well defined position in the discrete approximation.
# However, this may be the only practical solution in some complicated cases and
# is therefore allowed.
#
# To demonstrate the idea, we simply remove the lower fillet grid from the list 
# of components in the composite grid.
#
list-of-components
	remove-from-list lower-corner-grid
exit
pause
#
# Lets recompute the overlapping grid.
#
compute-overlap
pause
#
# Cleanup
#
reset
pause
#
# End of the internal flow example. Thank's for your attention!
#
pause
