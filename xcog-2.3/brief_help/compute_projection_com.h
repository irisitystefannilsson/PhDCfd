enum {LEFT, RIGHT, HELP, CANCEL};
const int LAST_COM = CANCEL;

char *COMMAND[4], *BRIEF_HELP[4];
int ARGUMENT[4], *SAVE_ON_COPY=NULL;
const int LEVEL=3;

COMMAND[LEFT]   = "project-low-r";
ARGUMENT[LEFT]   = 0;
BRIEF_HELP[LEFT]  = "Project the grid side r=0 onto the curve.";

COMMAND[RIGHT]  = "project-high-r";       
ARGUMENT[RIGHT]  = 0;
BRIEF_HELP[RIGHT] = "Project the grid side r=1 onto the curve.";

COMMAND[HELP]   = "help";
ARGUMENT[HELP]   = 1;
BRIEF_HELP[HELP]  = NULL;

COMMAND[CANCEL] = "cancel";
ARGUMENT[CANCEL] = 0;
BRIEF_HELP[CANCEL] = "Do not project any side onto the curve. Instead exit to "
"the previous command level.";


