enum{LOW_R, HIGH_R, LOW_S, HIGH_S, SHOW, HELP, EXIT};
const int LAST_COM = EXIT, LEVEL=1;
char *COMMAND[7], *BRIEF_HELP[7];
int ARGUMENT[7], *SAVE_ON_COPY=NULL;

COMMAND[LOW_R] ="low-r"; 
ARGUMENT[LOW_R] = 1;
BRIEF_HELP[LOW_R] = "Change the external segment for the low-r boundary, that is, "
"the side of the grid where the parameter r=0.";

COMMAND[HIGH_R] ="high-r"; 
ARGUMENT[HIGH_R] = 1;
BRIEF_HELP[HIGH_R] = "Change the external segment for the high-r boundary, that is, "
"the side of the grid where the parameter r=1.";

COMMAND[LOW_S] ="low-s"; 
ARGUMENT[LOW_S] = 1;
BRIEF_HELP[LOW_S] = "Change the external segment for the low-s boundary, that is, "
"the side of the grid where the parameter s=0.";

COMMAND[HIGH_S] ="high-s"; 
ARGUMENT[HIGH_S] = 1;
BRIEF_HELP[HIGH_S] = "Change the external segment for the high-s boundary, that is, "
"the side of the grid where the parameter s=1.";

COMMAND[SHOW] ="show";    
ARGUMENT[SHOW] = 0;
BRIEF_HELP[SHOW] = "Print the external sections on all sides of the grid.";

COMMAND[HELP] ="help";    
ARGUMENT[HELP] = 0;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT] ="exit";     
ARGUMENT[EXIT] = 0;
BRIEF_HELP[EXIT] = "Stop changing the external segments and exit from this "
"command level.";
