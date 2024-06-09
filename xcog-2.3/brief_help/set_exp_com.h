#define ncom 6
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] = "starting-grid-size";      argument[0] = 1;
command[1] = "reverse-focus";           argument[1] = 0;
command[2] = "grid-lines";              argument[2] = 1;

command[3] = "default-strength";        argument[3] = 0;

command[4] = "show-parameters";         argument[4] = 0;
command[5] = "help";                    argument[5] = 0;
command[6] = "exit";                    argument[6] = 0;

brief_help[0] = "Change the starting grid size. When there is no curve, the grid "
"size is computed under the assumtion that the grid line has unit length and is "
"uniformly parametrized."; 

brief_help[1] = "Focus the stretching around the other end point.";

brief_help[2] = "Change the number of grid lines.";

brief_help[3] = "Use the default parameters for the exponential stretching.";

brief_help[4] = "Show the strength of the exponential stretching.";

brief_help[5] = NULL;

brief_help[6] = "Stop changing the parameters of the exponential stretching "
"and exit from this command level.";

