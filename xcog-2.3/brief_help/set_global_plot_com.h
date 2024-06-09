enum{CURVES, CURVE_DIR, CURVE_TICK, CURVE_CONTROL_PTS, GRID_BOUNDARY, GRID_ARROWS,
     GRID_TICKS, ALL_GRID_LINES, BOUNDARY_COND, CURVE_LABEL, OVERLAPPING_GRID, INTERP_PTS,
     BAD_PTS, PHYS_BNDRY, NO_PHYS_BNDRY, FLAG, TITLE, C_AXIS, POSTSCRIPT, SHOW, CLEAR, 
     HELP, EXIT};
const int LAST_COM = EXIT, LEVEL=1, NO_INDENT=0;
char *command[23], *brief_help[23];
int argument[23], *save_on_copy=NULL;

command[CURVES]  ="curves";
argument[CURVES] = 0;
brief_help[CURVES] = "Toggle the plotting of the curve.";

command[CURVE_DIR]  ="curve-direction";
argument[CURVE_DIR] = 0;
brief_help[CURVE_DIR] = "Toggle the plotting of the arrow that shows the "
"direction of the parametrization of the curve.";

command[CURVE_TICK]  ="curve-tickmarks";
argument[CURVE_TICK] = 0;
brief_help[CURVE_TICK] = "Toggle the plotting of the tickmarks on the curve that "
"indicate the density of grid lines along the curve when the curve is used for "
"constructing a grid.";

command[CURVE_CONTROL_PTS]  ="curve-control-points";
argument[CURVE_CONTROL_PTS] = 0;
brief_help[CURVE_CONTROL_PTS] = "Toggle the plotting of the curve's control points. "
"The appearence of the control points depend on the curve type.";

command[GRID_BOUNDARY]  ="grid-boundary";
argument[GRID_BOUNDARY] = 0;
brief_help[GRID_BOUNDARY] = "Toggle the plotting of the outline of the grid.";

command[GRID_ARROWS]  ="grid-directional-arrows";
argument[GRID_ARROWS] = 0;
brief_help[GRID_ARROWS] = "Toggle the plotting of the arrows and labels that "
"indicate the direction of the r and s parameters of the mapping.";

command[GRID_TICKS]  ="grid-tickmarks";
argument[GRID_TICKS] = 0;
brief_help[GRID_TICKS] = "Toggle the plotting of the tickmarks that show "
"where the grid lines start from the boundary.";

command[ALL_GRID_LINES]  ="all-grid-lines";
argument[ALL_GRID_LINES] = 0;
brief_help[ALL_GRID_LINES] = "Toggle the plotting of the grid lines.";

command[BOUNDARY_COND]  ="boundary-condition";
argument[BOUNDARY_COND] = 0;
brief_help[BOUNDARY_COND] = "Toggle the labeling of the sides with the value of the "
"boundary condition.";

command[CURVE_LABEL]  ="curve-label";
argument[CURVE_LABEL] = 0;
brief_help[CURVE_LABEL] = "Toggle the labeling of the sides with the value of the curve "
"label, that is used to assist the overlap algorithm in setting up the "
"interpolation formula close to boundaies.";

command[OVERLAPPING_GRID] ="overlapping-grid";
argument[OVERLAPPING_GRID] = 0;
brief_help[OVERLAPPING_GRID] = "Toggle the plotting of the composite grid, i.e. all "
"grid lines used in the overlapping grid.";

command[INTERP_PTS] ="interpolation-points";
argument[INTERP_PTS] = 0;
brief_help[INTERP_PTS] = "Toggle the marking of all interpolation points in the composite "
"grid.";

command[BAD_PTS] ="bad-points";
argument[BAD_PTS] = 0;
brief_help[BAD_PTS] = "Toggle the marking of all bad points in the overlapping grid, "
"i.e. the points that don't satisfy the requirements for discretization or "
"interpolation points  in the following way:\n"
"A BOX indicates a bad discretization point.\n"
"A CROSS is a dead interpolation location.\n"
"A FILLED CIRCLE represents a non-explicit interpolation location.\n"
"A TRIANGLE marks a bad interpolation point because of a dead interpolation "
"location.\n"
"A DOWNWARD TRIANGLE indicates a bad interpolation point because of a non explicit "
"interpolation location.";

command[PHYS_BNDRY] ="physical-boundary";
argument[PHYS_BNDRY] = TRUE;
brief_help[PHYS_BNDRY] = "Show the physical boundary corresponding to one curve-label.";

command[NO_PHYS_BNDRY] ="no-physical-boundary";
argument[NO_PHYS_BNDRY] = FALSE;
brief_help[NO_PHYS_BNDRY] = "Turn off the plotting of physical boundaries.";

command[FLAG] ="flag-array";
argument[FLAG] = 0;
brief_help[FLAG] = "Toggle the labeling of all grid points used in the composite "
"grid by the value of the flag array.";

command[TITLE] ="title";
argument[TITLE] = 0;
brief_help[TITLE] = "Toggle the labeling of the `Composite Grid' window with a title. If the title is turned on, also enter the new new title string.";

command[C_AXIS] ="coordinate-axis";
argument[C_AXIS] = 0;
brief_help[C_AXIS] = "Toggle the plotting of coordinate axes.";

command[POSTSCRIPT] ="postscript-copy";
argument[POSTSCRIPT] = 1;
brief_help[POSTSCRIPT] = "Save the present content of the `Composite Grid' window in a "
"named file in postscript format.";

command[SHOW] ="show-plot-mode";
argument[SHOW] = 0;
brief_help[SHOW] = "Show the enabled and disabled plot options.";

command[CLEAR] ="clear-graphics";
argument[CLEAR] = 0;
brief_help[CLEAR] = "Turn off the plotting of all curves, mappings and the "
"composite grid.";

command[HELP] ="help";
argument[HELP] = 0;
brief_help[HELP] = NULL;

command[EXIT] ="exit";
argument[EXIT] = 0;
brief_help[EXIT] = "Stop modifying the plot options and exit from this command level.";
