#define ncom 6
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] ="low-r"; argument[0] = 1;
command[1] ="high-r"; argument[1] = 1;
command[2] ="low-s"; argument[2] = 1;

command[3] ="high-s"; argument[3] = 1;
command[4] ="show";    argument[4] = 0;
command[5] ="help";    argument[5] = 0;

command[6] ="exit";     argument[6] = 0;

brief_help[0] = "Change the value of the boundary condition on the low-r side of the "
"grid where the parameter r=0. The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is necessary information "
"for a PDE solver that uses the grid.";
brief_help[1] = "Change the value of the boundary condition on the high-r side of the "
"grid where the parameter r=1. The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is necessary information "
"for a PDE solver that uses the grid.";
brief_help[2] = "Change the value of the boundary condition on the low-s side of the "
"grid where the parameter s=0. The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is necessary information "
"for a PDE solver that uses the grid.";
brief_help[3] = "Change the value of the boundary condition on the high-s side of the "
"grid where the parameter s=1. The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is necessary information "
"for a PDE solver that uses the grid.";
brief_help[4] = "Print the boundary condition on all sides of the grid.";

brief_help[5] = NULL;

brief_help[6] = "Stop changing the boundary condition and exit from this "
"command level.";
