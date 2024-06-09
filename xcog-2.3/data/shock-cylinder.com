#xcog-file-type-1
reset-xcog
#
# The first curve describes the bow shock. We use a cubic spline to interpolate
# the given location of the spline.
#
make-curve bow-shock
	cubic-spline
		enter-new-nodes 3
			-.896542
			0
			3.10967
			4.476148
			8.4
			7.82295
pause
#
# Symmetry implies that the bow shock should be vertical at y=0. To achieve this we
# clamp the x-derivative of the parametric spline to zero at the start point.
#
		set-dx/dt(start) 0
pause
	exit
#
# The second curve describes the surface of the cylinder.
#
make-curve half-cylinder
	circular-arc
		end-angle 180
		radius .5
	exit
pause
#
# The third curve is a help curve for making a grid close to the bow shock. 
# Again, we use a cubic spline to interpolate the given location of the spline.
#
make-curve inner-curve
	cubic-spline
		enter-new-nodes 3
			-.6
			0.000000e+00
			1.725
			2.6
			8.4
			5.5
pause
#
# Symmetry implies that also the help curve should be vertical at y=0. 
# To achieve this we clamp the x-derivative of the parametric spline to zero at 
# the start point.
#
		set-dx/dt(start) 0
pause
	exit
#
# We use linear interpolation to discretize the domain close to the bow shock.
#
make-mapping shock-grid
	linear-interpolation-mapping
		inner-curve
		bow-shock
		set-r-lines 80
		set-s-lines 15
pause
#
# Redistribute the grid lines along the inner-curve by a hyperbolic tangent 
# stretching function.
#
		stretch-s=0-curve
			hyperbolic-tangent-stretching
			start-grid-size .03
			end-grid-size .2
			show-parameters
		exit
pause
#
# Also redistribute the grid lines along the bow-shock curve by a hyperbolic tangent 
# stretching function.
#
		stretch-s=1-curve
			hyperbolic-tangent-stretching
			start-grid-size .03
			end-grid-size .25
			show-parameters
		exit
		mapping-plot-mode
			all-grid-lines
		exit
	exit
pause
#
# The grid around the circle is constructed by taking normals out from the circle.
#
make-mapping cylinder-grid
	normal-curve-mapping half-cylinder
		reverse-curve-parametrization
		constant-width .3
		set-r-lines 60
		set-s-lines 25
pause
#
# Redistribute the grid lines along the inner-curve by an exponential
# stretching function.
#
		s-stretching
			exponential-stretching
			starting-grid-size 0.007
		exit
	exit
pause
#
# The solution behind the cylinder can be expected to be rather smooth, so a
# coarse Cartesian grid will resolve the solution adequately.
#
make-mapping base-grid
	cartesian-mapping
		xmin 0.9
		xmax 8.4
		ymin 0.0
		ymax 6.3
		set-r-lines 60
		set-s-lines 50
	exit
pause
#
# We fill the remaining hole with a Cartesian grid. To take advantage of the 
# fine resolution in the grids close to the shock and the cylinder, we will use
# a comparable resolution in this grid.
#
make-mapping middle-grid
	cartesian-mapping
		xmin -0.7
		xmax 1.5
		ymin 0.0
		ymax 2.9
		set-r-lines 45
		set-s-lines 55
	exit
pause
#
# 
# When a PDE is solved on the grid, the accuracy depends on the smoothness of the
# grid. A number of measures of the grid smoothness are evaluated by the 
# grid-quality command.
grid-quality shock-grid
pause
# 
# As can be seen by these numbers, this grid is rather smooth. The main problem is
# the large deviation from orthogonality between the grid lines.
#
# To more easily see the marked grid points, one can zoom in those regions of the
# Mapping window. To do this you place the cursor above and to the left of the 
# cylinder in the Mapping window. You then press the left mouse button and while
# keeping the mouse button down, drag the cursor to the right of and below the 
# cylinder, where you release the mouse button. The area that will be zoomed in is 
# enclose by "rubber bands". To reset the view to its original size, you press the 
# middle mouse button while the cursor is in the Mapping window. You can also 
# zoom out by using the left instead of the right mouse button and proceed as when 
# zooming in.
#
pause
#
# Label the sides that represent the boundary of the computational domain. 
# In this case, the domain is simply connected, so we use the same value of 
# the curve-label for all those sides.
#
curve-label shock-grid
	high-s 1
	low-r 1
exit
curve-label cylinder-grid
	low-r 1
	high-r 1
	low-s 1
exit
curve-label base-grid
	low-s 1
	high-r 1
exit
curve-label shock-grid
	high-r 1
exit
curve-label middle-grid
	low-s 1
exit
pause
#
# Make an external gap in the middle grid
#
mark-mapping-boundary middle-grid
	low-s
		0.08
		0.55
exit
pause
#
# Set the boundary condition for all boundaries on all component grids.
#
boundary-condition shock-grid
	high-s 1
	low-r 1
	high-r 1
exit
boundary-condition cylinder-grid
	low-r 1
	low-s 1
	high-r 1
exit
boundary-condition middle-grid
	low-s 1
exit
boundary-condition base-grid
	low-s 1
	high-r 1
exit
pause
#
# Modify the list of component grids in the overlapping grid to give the shock-grid
# the highest priority and the cylinder-grid the second highest priority.
#
list-of-components
	move-first cylinder-grid
	move-first shock-grid
	show
exit
pause
#
# We now construct the overlapping grid.
#
compute-overlap
pause
#
# Save the state of xcog, i.e. all curves, mappings and the value of all overlap
# parameters. After the program has been restarted or after a `reset-xcog', these
# quantities can be recovered by the command `read-state'. To save disc space, the
# overlapping grid is not saved and has to be recomputed.
#
save-state shock-cylinder.xcog
pause
#
# Save the composite grid on a HDF-file.
#
save-overlapping-grid 
	hdf-format shock-cylinder.hdf
pause
#
# Don't save the Jacobian of the transformation at all grid points.
#
	no
pause
#
# To demonstrate how Xcog behaves when there are not enough grid points in the overlap
# region, we will coarsen the grid around the cylinder and the middle-grid.
#
grid-lines middle-grid
	set-r-lines 25
	set-s-lines 25
exit
grid-lines cylinder-grid
	set-r-lines 30
	set-s-lines 5
exit
pause
#
# We will now attempt at making a composite grid using the same order of the component
# grids as in the previous composite grid. This should fail.
#
compute-overlap
pause
#
# The inconsistent grid points are in the bow-shock grid where it overlaps the 
# cylinder-grid and can more easily be seen by zooming in that region. 
#
# The overlap algorithm failed because there were not enough grid points in the overlap
# region between the cylinder grid and the bow-shock grid, which made it impossible
# to ensure explicit interpolation.
#
pause
#
# The number of grid points in the overlap region is still enough for implicit 
# interpolation and we will now change the interpolation type
#
overlap-parameters
	interpolation-type 
		implicit
exit
pause
#
# Making a coarse overlapping grid with implicit interpolation.
#
compute-overlap
pause
#
# Cleanup
#
reset-xcog
pause
#
# End of the shock and cylinder example. Thank's for your attention!
#
pause
