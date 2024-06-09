#define ncom 6
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] ="low-r"; 
argument[0] = 1;
brief_help[0] = "Change the curve value on the low-r side of the grid where the "
"parameter r=0.";

command[1] ="high-r"; 
argument[1] = 1;
brief_help[1] = "Change the curve value on the high-r side of the grid where the "
"parameter r=1.";

command[2] ="low-s"; 
argument[2] = 1;
brief_help[2] = "Change the curve value on the low-s side of the grid where the "
"parameter s=0.";

command[3] ="high-s"; 
argument[3] = 1;
brief_help[3] = "Change the curve value on the high-s side of the grid where the "
"parameter s=1.";

command[4] ="show";    
argument[4] = 0;
brief_help[4] = "Print the curve value for each side of the grid.";

command[5] ="help";    
argument[5] = 0;
brief_help[5] = NULL;

command[6] ="exit";    
argument[6] = 0;
brief_help[6] = "Stop changing the curve value and exit from this command level.";
