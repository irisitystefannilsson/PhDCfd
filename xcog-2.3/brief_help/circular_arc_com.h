#define ncom 7
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] = "start-angle"; argument[0] = 1;
command[1] = "end-angle";   argument[1] = 1;
command[2] = "radius";      argument[2] = 1;

command[3] = "origin";      argument[3] = 0;
command[4] = "curve-plot-mode";   argument[4] = 0;
/*file[4] = "set_curve_plot_com.h"; */

command[5] = "show";        argument[5] = 0;

command[6] = "help";        argument[6] = 0;
command[7] = "exit";        argument[7] = 0;

brief_help[0] = "Change the starting angle of the circular arc. The angle is "
"measured in degrees counter clockwise from the horizontal x-axis.";

brief_help[1] = "Change the ending angle of the circular arc. The angle is "
"measured in degrees counter clockwise from the horizontal x-axis.";

brief_help[2] = "Change the radius of the circular arc.";

brief_help[3] = "Change the Cartesian coordinates of the origin of the circular arc.";

brief_help[4] = "Set the graphical appearence in the `Curve' window or save the "
"graphical contents on a postscript file.";

brief_help[5] = "Display the value of the parameters describing the circular arc.";

brief_help[6] = NULL;

brief_help[7] = "Stop modifying the circular arc and exit from this command "
"level.";
