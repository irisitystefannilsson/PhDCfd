#define ncom 10
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] = "enter-new-nodes";  argument[0] = 1;
command[1] = "read-node-file";   argument[1] = 1;
command[2] = "change-node";     argument[2] = 1;

command[3] = "set-dx/dt(start)"; argument[3] = 1;
command[4] = "set-dy/dt(start)"; argument[4] = 1;
command[5] = "set-dx/dt(end)";   argument[5] = 1;

command[6] = "set-dy/dt(end)";   argument[6] = 1;
command[7] = "curve-plot-mode";        argument[7] = 0;
/*file[7] = "set_curve_plot_com.h"; */
command[8] = "show";             argument[8] = 0;

command[9] = "help";             argument[9] = 0;
command[10] = "exit";            argument[10] = 0;

brief_help[0] = "Enter a new set of node points for the spline. The old set of node "
"points will be discarded.";

brief_help[1] = "Read a new set of node points from an ASCII file. The old set of node "
"points will be discarded. The file has the following format:\n"
"#xcog-file-type-4\n"
"3\n"
"0.0 1.0\n"
"1.0 1.0\n"
"1.0 2.0\n"
"The first line of the file contains a descriptor which tells xcog that it is a node "
"point file. The second line contains the number of node points, which must be "
"greater than or equal to 2. On the following lines are the x and y coordinates for "
"the node points. The number of x and y coordinates must be at least equal to the "
"given number of node points, and every consecutive pair of node points must be "
"unique.";

brief_help[2] = "Choose one existing node point and change its Cartesian coordinates.";

brief_help[3] = "Specify the boundary condition for the derivative of the "
"x-coordinate of the parametric spline at the first node point of the spline. "
"A value greater than or equal to 1e30 changes the boundary condition to a "
"`free' condition, i.e. d^2 x/ dt^2 = 0.";

brief_help[4] = "Specify the boundary condition for the derivative of the "
"y-coordinate of the parametric spline at the first node point of the spline. "
"A value greater than or equal to 1e30 changes the boundary condition to a "
"`free' condition, i.e. d^2 y/ dt^2 = 0.";

brief_help[5] = "Specify the boundary condition for the derivative of the "
"x-coordinate of the parametric spline at the last node point of the spline. "
"A value greater than or equal to 1e30 changes the boundary condition to a "
"`free' condition, i.e. d^2 x/ dt^2 = 0.";

brief_help[6] = "Specify the boundary condition for the derivative of the "
"y-coordinate of the parametric spline at the last node point of the spline. "
"A value greater than or equal to 1e30 changes the boundary condition to a "
"`free' condition, i.e. d^2 y/ dt^2 = 0.";

brief_help[7] = "Set the graphical appearence in the `Curve' window or save the "
"graphical contents on a postscript file.";

brief_help[8] = "Display the value of the parameters describing the cubic spline.";

brief_help[9] = NULL;

brief_help[10] = "Stop modifying the cubic spline and exit from this command "
"level.";
