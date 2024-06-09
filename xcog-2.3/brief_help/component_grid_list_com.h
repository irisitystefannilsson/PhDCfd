#define ncom 8
char *command[ncom+1], *brief_help[ncom+1];
int argument[ncom+1], *save_on_copy=NULL;

command[0] = "insert-all-mappings"; argument[0] = 0;
command[1] = "insert-first";        argument[1] = 1;
command[2] = "remove-from-list";    argument[2] = 1;
     
command[3] = "move-first";          argument[3] = 1;
command[4] = "move-last";           argument[4] = 1;
command[5] = "move-up";             argument[5] = 1;
     
command[6] = "show";                argument[6] = 0;
command[7] = "help";                argument[7] = 0;
command[8] = "exit";                argument[8] = 0;

brief_help[0] = "Insert all defined mappings in order into the list of mappings to "
"be used in the composite grid.";

brief_help[1] = "Select one mapping and insert it first in the list of mappings "
"to be used in the composite grid. The first mapping has the highest priority in "
"the composite grid.";

brief_help[2] = "Remove one mapping from the list of mappings to be used in "
"the composite grid.";

brief_help[3] = "Select one mapping in the composite grid list and move it first in "
"that list. The first mapping has the highest priority in the composite grid.";

brief_help[4] = "Select one mapping in the composite grid list and move it last in "
"that list. The last mapping has the lowest priority in the composite grid.";

brief_help[5] = "Select one mapping in the composite grid list and move it one "
"position up in that list. This increases the priority level of that mapping in "
"the composite grid by one.";

brief_help[6] = "Present the list of mappings in the composite grid and the value "
"of the overlap parameters, i.e. the discretization width, etc.";

brief_help[7] = NULL;

brief_help[8] = "Stop modifying the list of component grids that will be included in the overlapping grid and return to the previous command level.";
