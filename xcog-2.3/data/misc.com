#xcog-file-type-1
reset-xcog
#
# Read in a saved state. The file `cylinder-square.xcog' is created by the
# cylinder-in-square tutorial, and we assume that this already has been done. If not,
# you can stop this tutorial now by giving the `break' command and try again after the
# file has been created.
#
pause
read-state cylinder-square.xcog
#
# Compute the overlap, so that we have something interesting to look at.
#
pause
compute-overlap
#
# Save the contents of the Composite grid window on a postscript file. Observe that
# it is possible to zoom in or out to modify the plot.
#
pause
composite-plot-mode
	postscript-copy intro.ps
#
# Save color information. On some black and white printers, the result gets better
# if the color information is NOT saved.
#
pause
	yes
#
# Clear the Composite grid window.
#
pause
	clear-graphics
#
# Show the defined curves with tickmarks that indicate the parametrization.
#
pause
	curves
	curve-tickmarks
#
# Turn off the coordinate axis
#
pause
	coordinate-axis
#
# Turn on the title
#
pause
	title My_xcog_data
#
# Clear the Composite grid window.
#
pause
	clear-graphics
#
# Show the mappings with directional arrows that shows the direction of
# the parametrization.
#
pause
	grid-directional-arrows
	all-grid-lines
#
# Clear the Composite grid window.
#
pause
	clear-graphics
#
# Show only the interpolation points
#
pause
	interpolation-points
#
# Turn on the grid points in the overlapping grid
#
pause
	overlapping-grid
#
# Turn on the coordinate axis
#
pause
	coordinate-axis
#
# End of the tour of plot-mode.
#
pause	
exit
#
# To look at a specific mapping and to save it on a postscript file, we need to first
# specify the mapping. 
#
pause
change-mapping cylinder-grid
#
# To save the mapping on a postscript file we enter plot-mode.
#
pause
		mapping-plot-mode
#
# Turn on all grid lines to have something more substantial to look at.
#
pause
			all-grid-lines
#
# Save the contents of the Mapping window on a postscript file. Observe that
# it is possible to zoom in or out to modify the plot.
#
pause
			postscript-copy mapping.ps
#
# Save color information. On some black and white printers, the result gets better
# if the color information is NOT saved.
#
pause
			yes
#
# Turn on the title which displays the name of the mapping
#
pause
			title
		exit
	exit
#
# To look at a specific curve and to save it on a postscript file, we need to first
# specify the curve. 
#
pause
change-curve cylinder
#
# To save the mapping on a postscript file we enter plot-mode.
#
pause
		curve-plot-mode
#
# Save the contents of the Curve window on a postscript file. Observe that
# it is possible to zoom in or out to modify the plot.
#
pause
			postscript-copy curve.ps
#
# Save color information. On some black and white printers, the result gets better
# if the color information is NOT saved.
#
pause
			yes
		exit
	exit
#
# End of demo.
#
# Cleanup
#
reset-xcog
