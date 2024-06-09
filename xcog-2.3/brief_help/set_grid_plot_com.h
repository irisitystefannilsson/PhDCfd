#define ncom 11
char *command[ncom+1], *brief_help[ncom+1];
int argument[ncom+1], *save_on_copy=NULL;

command[0] ="grid-boundary";           argument[0]=0;
command[1] ="grid-directional-arrows"; argument[1]=0;
command[2] ="grid-tickmarks";          argument[2]=0;

command[3] ="all-grid-lines";          argument[3]=0;
command[4] ="boundary-condition";      argument[4]=0;
command[5] ="curve-value";             argument[5]=0;

command[6] ="title";                   argument[6]=0;
command[7] ="coordinate-axis";         argument[7]=0;
command[8] ="postscript-copy";         argument[8]=1;

command[9] ="show";                    argument[9]=0;
command[10] ="help";                   argument[10]=0;
command[11] ="exit";                   argument[11]=0;

brief_help[0] = "Toggle the plotting of the outline of the grid.";

brief_help[1] = "Toggle the plotting of the arrows and labels that indicate the "
"direction of the r and s parameters of the mapping.";

brief_help[2] = "Toggle the plotting of the tickmarks that show where the grid "
"lines start from the boundary.";

brief_help[3] = "Toggle the plotting of the grid lines.";

brief_help[4] = "Toggle the labeling of the sides with the value of the "
"boundary condition.";

brief_help[5] = "Toggle the labeling of the sides with the value of the curve "
"variable, that is used to assist the overlap algorithm in setting up the "
"interpolation formula close to boundaies.";

brief_help[6] = "Toggle the plotting of a title that shows the name of the mapping.";

brief_help[7] = "Toggle the plotting of coordinate axes.";

brief_help[8] = "Save the present content of the `Mapping' window in a named file "
"in postscript format.";

brief_help[9] = "Show the enabled and disabled plot options.";

brief_help[10] = NULL;

brief_help[11] = "Stop modifying the plot options and exit from this command level.";
