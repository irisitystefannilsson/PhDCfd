#define ncom 5
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] ="start-point"; argument[0]=1;
command[1] ="end-point";   argument[1]=1;
command[2] ="curve-plot-mode";   argument[2] = 0;
/*file[2] = "set_curve_plot_com.h"; */

command[3] ="show";        argument[3] = 0;
command[4] ="help";        argument[4] = 0;
command[5] ="exit";        argument[5] = 0;

brief_help[0] = "Change the Cartesian coordinate of the starting point of the "
"straight line.";

brief_help[1] = "Change the Cartesian coordinate of the ending point of the "
"straight line.";

brief_help[2] = "Set the graphical appearence in the `Curve' window or save the "
"graphical contents on a postscript file.";

brief_help[3] = "Display the value of the parameters describing the straight line.";

brief_help[4] = NULL;

brief_help[5] = "Stop modifying the straight line and exit from this command "
"level.";
