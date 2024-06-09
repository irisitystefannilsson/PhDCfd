#define ncom 4
char *command[ncom+1], *brief_help[ncom+1];
int *argument=NULL, *save_on_copy=NULL;

command[0] ="explicit";
command[1] ="implicit";
command[2] = "show";
command[3] = "help";
command[4] = "cancel";

brief_help[0] = "When the interpolation is explicit, the overlap between the grids "
"is sufficiently wide to ensure that the interpolation formula does NOT include "
"interpolation points from other grids.";

brief_help[1] = "When the interpolation is implicit, the interpolation formula is "
"allowed to include interpolation points from other grids. Therefore, the overlap "
"can usually be narrower than when explicit interpolation is used.";

brief_help[2] = "Print the present interpolation type.";

brief_help[3] = NULL;

brief_help[4] = "Exit from this command level without changing the "
"interpolation type.";
