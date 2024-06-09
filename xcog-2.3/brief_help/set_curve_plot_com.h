#define ncom 9
char *command[ncom+1], *brief_help[ncom+1];
int argument[ncom+1], *save_on_copy=NULL;

command[0] ="curve";                     argument[0]=0;
command[1] ="directional-arrow";         argument[1]=0;
command[2] ="parametrization-tickmarks"; argument[2]=0;

command[3] ="control-points";            argument[3]=0;
command[4] ="coordinate-axis";           argument[4]=0;
command[5] ="title";                     argument[5]=0;

command[6] ="postscript-copy";           argument[6]=1;
command[7] ="show";                      argument[7]=0;
command[8] ="help";                      argument[8]=0;

command[9] ="exit";                      argument[9]=0;

brief_help[0] = "Toggle the plotting of the curve.";

brief_help[1] = "Toggle the plotting of the arrow that shows the direction of the "
"parametrization of the curve.";

brief_help[2] = "Toggle the plotting of the tickmarks on the curve that indicate "
"the density of grid lines along the curve when the curve is used for "
"constructing a grid.";

brief_help[3] = "Toggle the plotting of the control points along the curve. For "
"a straight line or a circular arc, the 2 control points are the endpoints. For "
"a cubic spline, the control points are the node points of the interpolating "
"spline and for a smoothed polygon, the control points are the corners of the "
"polygon.";

brief_help[4] = "Toggle the plotting of coordinate axes.";

brief_help[5] = "Toggle the plotting of a title that shows the name of the curve.";

brief_help[6] = "Save the present content of the `Curve' window in a named file "
"in postscript format.";

brief_help[7] = "Show the enabled and disabled plot options.";

brief_help[8] = NULL;

brief_help[9] = "Stop modifying the plot options and exit from this command level.";
