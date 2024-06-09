#define ncom 13
char *command[ncom+1], *brief_help[ncom+1];
int argument[ncom+1], *save_on_copy=NULL;

for (i=0; i<=ncom; i++) argument[i]=0;

command[0] = "xmin";         argument[0] = 1;
command[1] = "xmax";         argument[1] = 1;
command[2] = "ymin";         argument[2] = 1;

command[3] = "ymax";         argument[3] = 1;
command[4] = "r-stretching";     
/*file[4] = "choose_stretching_com.h"; */
command[5] = "no-r-stretching";

command[6] = "s-stretching";              
/*file[6] = "choose_stretching_com.h"; */
command[7] = "no-s-stretching";
command[8] = "set-r-lines";  argument[8] = 1;

command[9] = "set-s-lines";  argument[9] = 1;
command[10] = "mapping-plot-mode";
/*file[10] = "set_grid_plot_com.h"; */

command[11] ="show";    
command[12] ="help";
command[13] ="exit";     

brief_help[0] = "Change the horizontal starting coordinate of the Cartesian grid.";

brief_help[1] = "Change the horizontal ending coordinate of the Cartesian grid.";

brief_help[2] = "Change the vertical starting coordinate of the Cartesian grid.";

brief_help[3] = "Change the vertical ending coordinate of the Cartesian grid.";

brief_help[4] = "Enable stretching of the grid points in the r-direction.";

brief_help[5] = "Disable the stretching in the r-direction.";

brief_help[6] = "Enable stretching of the grid points in the s-direction.";

brief_help[7] = "Disable stretching in the s-direction.";

brief_help[8] = "Change the number of grid lines in the r-direction of the grid.";

brief_help[9] = "Change the number of grid lines in the s-direction of the grid.";

brief_help[10] = "Set the graphical appearence in the `Mapping' window or save the "
"graphical contents on a postscript file.";

brief_help[11] = "Display the value of the parameters describing the Cartesian "
"mapping.";

brief_help[12] = NULL;

brief_help[13] = "Stop modifying the Cartesian mapping and exit from this command "
"level.";

