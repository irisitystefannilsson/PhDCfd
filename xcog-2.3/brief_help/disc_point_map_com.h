enum{R_STRETCH, NO_R_STRETCH, S_STRETCH, NO_S_STRETCH, R_LINES, S_LINES, PLOT_MODE, SHOW,
     ORDER, HELP, EXIT};
#define ncom 10
char *command[ncom+1], *brief_help[ncom+1];
int argument[ncom+1], *save_on_copy=NULL;

for (i=0; i<=ncom; i++) argument[i]=0;

command[ORDER] = "interpolation-order";
brief_help[ORDER] = "Set the interpolation scheme to be linear or cubic."; 

command[R_STRETCH] = "r-stretching";     
/*file[R_STRETCH] = "choose_stretching_com.h"; */
brief_help[R_STRETCH] = "Enable stretching of the grid points in the r-direction.";

command[NO_R_STRETCH] = "no-r-stretching";
brief_help[NO_R_STRETCH] = "Disable the stretching in the r-direction.";

command[S_STRETCH] = "s-stretching";              
/*file[S_STRETCH] = "choose_stretching_com.h"; */
brief_help[S_STRETCH] = "Enable stretching of the grid points in the s-direction.";

command[NO_S_STRETCH] = "no-s-stretching";
brief_help[NO_S_STRETCH] = "Disable stretching in the s-direction.";

command[R_LINES] = "set-r-lines";  argument[R_LINES] = 1;
brief_help[R_LINES] = "Change the number of grid lines in the r-direction of the grid.";

command[S_LINES] = "set-s-lines";  argument[S_LINES] = 1;
brief_help[S_LINES] = "Change the number of grid lines in the s-direction of the grid.";

command[PLOT_MODE] = "mapping-plot-mode";
/*file[PLOT_MODE] = "set_grid_plot_com.h"; */
brief_help[PLOT_MODE] = "Set the graphical appearence in the `Mapping' window or save the "
"graphical contents on a postscript file.";

command[SHOW] = "show";    
brief_help[SHOW] = "Display the value of the parameters describing the discrete point "
"mapping.";

command[HELP] ="help";
brief_help[HELP] = NULL;

command[EXIT] ="exit";     
brief_help[EXIT] = "Stop modifying the discrete point mapping and exit from this "
"command level.";

