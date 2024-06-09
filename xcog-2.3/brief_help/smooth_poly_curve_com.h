#define ncom 6
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] = "enter-new-corners";     argument[0] = 1;
command[1] = "change-corner";         argument[1] = 1;
command[2] = "corner-sharpness";      argument[2] = 1;

command[3] = "curve-plot-mode";             argument[3] = 0;
/*file[3] = "set_curve_plot_com.h"; */
command[4] = "show";                  argument[4] = 0;
command[5] = "help";                 argument[5] = 0;

command[6] = "exit";                  argument[6] = 0;

brief_help[0] = "Enter a new set of corners for the smoothed polygon. "
"The old set of corners will be discarded.";

brief_help[1] = "Choose one existing corner and change its Cartesian coordinates.";

brief_help[2] = "Change the parameter that describes how close the smoothed polygon "
"is to the straight lines between the corners. A larger value means that the "
"smoothed polygon will be closer to the straight line polygon. "
"Normally, the parameter should have a value between 1 and 10.";

brief_help[3] = "Set the graphical appearence in the `Curve' window or save the "
"graphical contents on a postscript file.";

brief_help[4] = "Display the value of the parameters describing the smoothed polygon.";

brief_help[5] = NULL;

brief_help[6] = "Stop modifying the smoothed polygon and exit from this command "
"level.";
