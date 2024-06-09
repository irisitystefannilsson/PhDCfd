#define ncom 4
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] ="set-r-lines"; argument[0]=1;
command[1] ="set-s-lines"; argument[1]=1;
command[2] ="show";        argument[2]=0;

command[3] ="help";        argument[3]=0;
command[4] ="exit";        argument[4]=0;

brief_help[0] = "Change the number of grid lines in the r-direction of the grid.";

brief_help[1] = "Change the number of grid lines in the s-direction of the grid.";

brief_help[2] = "Print the current number of grid lines in both directions of "
"the grid.";

brief_help[3] = NULL;

brief_help[4] = "Stop changing the number of grid lines and exit from this "
"command level.";
