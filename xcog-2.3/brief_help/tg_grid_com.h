enum {TE_RADIUS, S_SCALING, QUALITY, R_STRETCH, NO_R_STRETCH, S_STRETCH, NO_S_STRETCH, 
      R_LINES, S_LINES, PLOT_MODE, SHOW, HELP, EXIT};

const int LAST_COM = EXIT, LEVEL=2;
char *COMMAND[13], *BRIEF_HELP[13];
int ARGUMENT[13], *SAVE_ON_COPY = NULL;

COMMAND[TE_RADIUS] = "trailing-edge-radius";    ARGUMENT[TE_RADIUS] = TRUE;
BRIEF_HELP[TE_RADIUS] = "Parameter related to the trailing edge radius. "
"(the rest of the documentation is sitting in Bjorn Regnstrom's chair).";

COMMAND[S_SCALING] = "s-scaling";    ARGUMENT[S_SCALING] = TRUE;
BRIEF_HELP[S_SCALING] = "Parameter related to the width of the mapping. A larger value "
"makes the mapping thicker. "
"(the rest of the documentation is sitting in Bjorn Regnstrom's chair).";

COMMAND[QUALITY] = "grid-quality";          
ARGUMENT[QUALITY] = TRUE;
BRIEF_HELP[QUALITY] = "Evaluate a number of measures of the quality of a grid "
"associated with a mapping.";

COMMAND[R_STRETCH] = "r-stretching";     ARGUMENT[R_STRETCH] = FALSE;
/*file[R_STRETCH] = "choose_stretching_com.h"; */
BRIEF_HELP[R_STRETCH] = "Enable stretching of the grid points in the r-direction.";

COMMAND[NO_R_STRETCH] = "no-r-stretching"; ARGUMENT[NO_R_STRETCH] = FALSE;
BRIEF_HELP[NO_R_STRETCH] = "Disable the stretching in the r-direction.";

COMMAND[S_STRETCH] = "s-stretching"; ARGUMENT[S_STRETCH] = FALSE;             
/*file[S_STRETCH] = "choose_stretching_com.h"; */
BRIEF_HELP[S_STRETCH] = "Enable stretching of the grid points in the s-direction.";

COMMAND[NO_S_STRETCH] = "no-s-stretching"; ARGUMENT[NO_S_STRETCH] = FALSE;
BRIEF_HELP[NO_S_STRETCH] = "Disable stretching in the s-direction.";

COMMAND[R_LINES] = "set-r-lines";  ARGUMENT[R_LINES] = TRUE;
BRIEF_HELP[R_LINES] = "Change the number of grid lines in the r-direction of the grid.";

COMMAND[S_LINES] = "set-s-lines";  ARGUMENT[S_LINES] = TRUE;
BRIEF_HELP[S_LINES] = "Change the number of grid lines in the s-direction of the grid.";

COMMAND[PLOT_MODE] = "mapping-plot-mode"; ARGUMENT[PLOT_MODE] = FALSE;
/*file[PLOT_MODE] = "set_grid_plot_com.h"; */
BRIEF_HELP[PLOT_MODE] = "Set the graphical appearence in the `Mapping' window or save "
"the graphical contents on a postscript file.";

COMMAND[SHOW] ="show";  ARGUMENT[SHOW] = FALSE;  
BRIEF_HELP[SHOW] = "Display the value of the parameters describing the Thedorsen-"
"Garrick mapping.";

COMMAND[HELP] ="help"; ARGUMENT[HELP] = TRUE;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT] ="exit";  ARGUMENT[EXIT] = FALSE;
BRIEF_HELP[EXIT] = "Stop changing the mapping and exit to the previous command level";
