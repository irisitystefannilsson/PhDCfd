#define ncom 15
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] = "change-curve";           argument[0] = 0;
command[1] = "reverse-curve-parametrization"; argument[1] = 0;
command[2] = "constant-width";         argument[2] = 1;

command[3] = "variable-width";         argument[3] = 0;
command[4] = "change-variable-width";  argument[4] = 1;
command[5] = "width-change-sharpness"; argument[5] = 1;

command[6] = "s-stretching";           argument[6] = 0;  
/*file[6] = "choose_stretching_com.h"; */

command[7] = "no-s-stretching";        argument[7] = 0;  

command[8] = "r-stretching";           argument[8] = 0;
/*file[8] = "choose_stretching_com.h"; */

command[9] = "no-r-stretching";        argument[9] = 0;
command[10] = "set-r-lines";           argument[10] = 1;
command[11] = "set-s-lines";           argument[11] = 1;

command[12] = "mapping-plot-mode";             argument[12] = 0;
/*file[12] = "set_grid_plot_com.h"; */

command[13] = "show";                  argument[13] = 0;
command[14] = "help";                  argument[14] = 0;

command[15] = "exit";                  argument[15] = 0;

brief_help[0] = "Change the shape of the boundary curve.";

brief_help[1] = "Flip the curve parametrization to make the normals grow out in the "
"opposite direction.";

brief_help[2] = "Let the normals grow out a constant distance from the curve and "
"specify this distance."; 

brief_help[3] = "Let the length of the normals vary between each control "
"point of the specified curve and specify the length at each control point. "
"If the curve is a straight line, the 2 control points "
"are the starting and ending points. If the curve is a cubic spline, the control "
"points are the node points. If the curve is a smoothed polygon, the control points "
"are the corners of the polygon and if the curve is a circular arc, the 2 control "
"points are the starting and ending angles.";

brief_help[4] = "Change the width of the normals in a specified interval I, i.e. "
"between control points I and I+1.";

brief_help[5] = "Change the parameter that describes how sharply the length of the "
"normals should change when the `variable-width' option is choosen. The value should "
"normally be between 1 and 10, where a larger value gives a sharper transition.";

brief_help[6] = "Enable stretching of the grid points in the s-direction, i.e., "
"normal to the curve.";

brief_help[7] = "Disable stretching in the s-direction.";

brief_help[8] = "Enable stretching of the grid points in the r-direction, i.e., "
"tangentially along the curve.";

brief_help[9] = "Disable the stretching in the r-direction.";

brief_help[10] = "Change the number of grid lines in the r-direction of the grid.";

brief_help[11] = "Change the number of grid lines in the s-direction of the grid.";

brief_help[12] = "Set the graphical appearence in the `Mapping' window or save the "
"graphical contents on a postscript file.";

brief_help[13] = "Display the value of the parameters describing the Cartesian "
"mapping.";

brief_help[14] = NULL;

brief_help[15] = "Stop modifying the normal curve mapping and exit from this command "
"level.";

