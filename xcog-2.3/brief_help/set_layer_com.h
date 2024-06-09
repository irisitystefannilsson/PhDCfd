#define ncom 6
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] = "change-strength";        argument[0] = 1;
command[1] = "change-width";           argument[1] = 1;
command[2] = "default-strength-width"; argument[2] = 0;

command[3] = "grid-lines";             argument[3] = 1;

command[4] = "show-parameters";        argument[4] = 0;
command[5] = "help";                   argument[5] = 0;
command[6] = "exit";                   argument[6] = 0;

brief_help[0] = "Change the strength of the layer stretching at one node point.";

brief_help[1] = "Change the width of the layer stretching at one node point.";

brief_help[2] = "Use the default parameters for the layer stretching.";

brief_help[3] = "Change the number of grid lines.";

brief_help[4] = "Show the strength and width of the layer stretching at each node "
"point.";

brief_help[5] = NULL;

brief_help[6] = "Stop changing the parameters of the layer stretching and exit "
"from this command level.";

