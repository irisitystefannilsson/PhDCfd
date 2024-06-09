enum {LOW_R1, HIGH_R1, SHOW, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[5], *BRIEF_HELP[5];
int ARGUMENT[5], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[LOW_R1]  ="low-r";  
ARGUMENT[LOW_R1]  = 1;
BRIEF_HELP[LOW_R1] = 
"Change the boundary condition for the r=0 side.\n"
"The following boundary conditions are available:\n"
"free\n"
"x-constant\n"
"y-constant\n";

COMMAND[HIGH_R1] ="high-r"; 
ARGUMENT[HIGH_R1] = 1;
BRIEF_HELP[HIGH_R1] = 
"Change the boundary condition for the r=1 side.\n"
"The following boundary conditions are available:\n"
"free\n"
"x-constant\n"
"y-constant\n";

COMMAND[SHOW]  ="show";      
ARGUMENT[SHOW]  = 0;
BRIEF_HELP[SHOW] = "Print the boundary condition on all sides of the grid.";

COMMAND[HELP]  ="help";      
ARGUMENT[HELP]  = 0;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT]  ="exit";      
ARGUMENT[EXIT]  = 0;
BRIEF_HELP[EXIT] = "Stop changing the boundary condition and exit from this "
"command level.";
