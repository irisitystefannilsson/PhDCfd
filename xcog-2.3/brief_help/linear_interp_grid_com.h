#define ncom 15
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] = "change-s=0-curve";   argument[0] = 0;
command[1] = "reverse-s=0-curve";  argument[1] = 0;
command[2] = "stretch-s=0-curve";  argument[2] = 0;
/*file[2] = "choose_stretching_com.h"; */

command[3] = "no-s=0-stretching";  argument[3] = 0;

command[4] = "change-s=1-curve";   argument[4] = 0;
command[5] = "reverse-s=1-curve";  argument[5] = 0;

command[6] = "stretch-s=1-curve";  argument[6] = 0;
/*file[6] = "choose_stretching_com.h"; */

command[7] = "no-s=1-stretching";  argument[7] = 0;

command[8] = "s-stretching";       argument[8] = 0;
/*file[8] = "choose_stretching_com.h"; */

command[9] = "no-s-stretching";    argument[9] = 0;
command[10] = "set-r-lines";       argument[10] = 1;
command[11] = "set-s-lines";       argument[11] = 1;

command[12] = "mapping-plot-mode";         argument[12] = 0;
/*file[12] = "set_grid_plot_com.h"; */

command[13] = "show";              argument[13] = 0;
command[14] = "help";              argument[14] = 0;

command[15] = "exit";              argument[15] = 0;

brief_help[0] = "Change the shape of the s=0 curve.";

brief_help[1] = "Flip the curve parametrization of the s=0 curve.";

brief_help[2] = "Enable stretching of the grid points in the r-direction along "
"the s=0 curve.";

brief_help[3] = "Disable stretching of the s=0-curve.";

brief_help[4] = "Change the shape of the s=1 curve.";

brief_help[5] = "Flip the curve parametrization of the s=1 curve.";

brief_help[6] = "Enable stretching of the grid points in the r-direction along "
"the s=1 curve.";

brief_help[7] = "Disable stretching of the s=1-curve.";

brief_help[8] = "Enable stretching of the grid points in the s-direction between "
"the boundary curves.";

brief_help[9] = "Disable stretching in the s-direction.";

brief_help[10] = "Change the number of grid lines in the r-direction of the grid.";

brief_help[11] = "Change the number of grid lines in the s-direction of the grid.";

brief_help[12] = "Set the graphical appearence in the `Mapping' window or save the "
"graphical contents on a postscript file.";

brief_help[13] = "Display the value of the parameters describing the linear "
"interpolation mapping.";

brief_help[14] = NULL;

brief_help[15] = "Stop modifying the linear interpolation mapping and exit from "
"this command level.";
