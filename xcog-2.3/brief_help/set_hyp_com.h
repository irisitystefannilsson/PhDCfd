enum {BOTH_SIDES, PROJECT_SIDE, COMPUTE, CURV_FACTOR, SMOOTH_FACTOR, V_MIN, WIDTH, 
      C_GRID, N_CUT, BRANCH_LOC, BOUNDARY_CONDITION, R1_LINES, R2_LINES, S_STRETCHING, 
      NO_S_STRETCHING, R_STRETCHING, NO_R_STRETCHING, CHANGE_CURVE, FLIP, QUALITY, 
      SHOW, PLOT_MODE, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[24], *BRIEF_HELP[24];
int ARGUMENT[24], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[C_GRID]  = "c-grid";
ARGUMENT[C_GRID]  = 0;
BRIEF_HELP[C_GRID] = "Toggle the c-grid mode on or off.";

COMMAND[N_CUT]  = "branch-points";
ARGUMENT[N_CUT]  = 1;
BRIEF_HELP[N_CUT] = "Number of grid points on the branch cut in a c-grid.";

COMMAND[BRANCH_LOC]  = "end-location";
ARGUMENT[BRANCH_LOC]  = 0;
BRIEF_HELP[BRANCH_LOC] = "The starting location for the branch-cut, relative "
"to the periodic point on the curve. (Only applicable for c-grids.)";

COMMAND[PROJECT_SIDE]  = "project-side";
/*file[PROJECT_SIDE]  = "compute_projection_com.h"; */
ARGUMENT[PROJECT_SIDE]  = 1;
BRIEF_HELP[PROJECT_SIDE] = "Project one side of the mapping onto a curve.";

COMMAND[BOTH_SIDES] = "both-sided-grid";
ARGUMENT[BOTH_SIDES] = 0;
BRIEF_HELP[BOTH_SIDES] = "Toggle on/off growing the grid in both directions from the "
"curve. This could be useful for constructing a grid around a one-dimensional feature "
"of the solution, or for making a grid around an infinitely thin body, like a sail.";

COMMAND[QUALITY] = "evaluate-grid-quality";
ARGUMENT[QUALITY] = 0;
BRIEF_HELP[QUALITY] = "Toggle the flag that determines wether the quality of "
"the grid should be plotted or not.";

COMMAND[CHANGE_CURVE] = "change-curve";
ARGUMENT[CHANGE_CURVE] = 0;
BRIEF_HELP[CHANGE_CURVE] = "Change the shape of the boundary curve.";

COMMAND[FLIP] = "reverse-curve-parametrization";
ARGUMENT[FLIP] = 0;
BRIEF_HELP[FLIP] = "Flip the curve parametrization to make the grid grow out in the "
"opposite direction.";

COMMAND[BOUNDARY_CONDITION]  = "boundary-condition";
ARGUMENT[BOUNDARY_CONDITION]  = 0;
BRIEF_HELP[BOUNDARY_CONDITION] = "Set the boundary condition for the hyperbolic "
"grid generation.";

COMMAND[COMPUTE]    = "compute-mapping";  
ARGUMENT[COMPUTE]   = 0;
BRIEF_HELP[COMPUTE] = "Solve the hyperbolic PDE with the current parameter values.";

COMMAND[CURV_FACTOR]    = "curvature-coefficient";  
ARGUMENT[CURV_FACTOR]   = 1;
BRIEF_HELP[CURV_FACTOR] = "Set the coefficient e in (1 - e*k), which is "
"the propagation speed of the grid surface. In this expression, k=(k1+k2)/2 "
"is the mean curvature of the grid surface at the grid point. Also see the command "
"`velocity-threshold'.";

COMMAND[SMOOTH_FACTOR]    = "averaging-coefficient";  
ARGUMENT[SMOOTH_FACTOR]   = 1;
BRIEF_HELP[SMOOTH_FACTOR] = "Set the coefficient for the spacial "
"averaging. The averaging is done on the velocity * normal. A factor of 0.0 "
"means no averaging, so the velocity * normal is determined pointwise. A value "
"of 1.0 implies that the point value is replaced by the algebraic average of the "
"nearest neighbors. A value inbetween 0 and 1 implies a linear interpolation "
"between the two extrema.";

COMMAND[V_MIN]    = "velocity-threshold";  
ARGUMENT[V_MIN]   = 1;
BRIEF_HELP[V_MIN] = "Set the minimum propagation speed of the grid surface. This "
"value is used near convex corners, where the propagation speed otherwise could "
"go negative.";

COMMAND[WIDTH]    = "thickness";  
ARGUMENT[WIDTH]   = 1;
BRIEF_HELP[WIDTH] = "Set the approximate thickness of the grid in the r3-direction. "
"The actual thickness will vary and be smaller where the surface is convex and be "
"larger where it is concave.";

COMMAND[R1_LINES]    = "set-r-lines";
ARGUMENT[R1_LINES]   = 1;
BRIEF_HELP[R1_LINES] = "Change the number of grid lines in the r-direction.";

COMMAND[R2_LINES]    = "set-s-lines";
ARGUMENT[R2_LINES]   = 1;
BRIEF_HELP[R2_LINES] = "Change the number of grid lines in the s-direction.";

COMMAND[S_STRETCHING]    = "s-stretching";
/*file[S_STRETCHING] = "choose_stretching_com.h"; */
ARGUMENT[S_STRETCHING]   = 0;
BRIEF_HELP[S_STRETCHING] = "Introduce a stretching function in the normal direction, "
"or, if a stretching function already is in use, modify its parameters.";

COMMAND[NO_S_STRETCHING]    = "no-s-stretching";
ARGUMENT[NO_S_STRETCHING]   = 0;
BRIEF_HELP[NO_S_STRETCHING] = "Disable the stretching function in the normal direction "
"and use a uniform distribution of grid points instead.";

COMMAND[R_STRETCHING]    = "r-stretching";
/*file[R_STRETCHING] = "choose_stretching_com.h"; */
ARGUMENT[R_STRETCHING]   = 0;
BRIEF_HELP[R_STRETCHING] = "Introduce a stretching function in the tangential "
"direction, or, if a stretching function already is in use, modify its parameters.";

COMMAND[NO_R_STRETCHING]    = "no-r-stretching";
ARGUMENT[NO_R_STRETCHING]   = 0;
BRIEF_HELP[NO_R_STRETCHING] = "Disable the stretching function in the tangential "
"direction and use a uniform distribution of grid points instead.";

COMMAND[SHOW]    = "show-parameters";
ARGUMENT[SHOW]   = 0;
BRIEF_HELP[SHOW] = "Display the parameters that control the mapping.";

COMMAND[PLOT_MODE]    = "mapping-plot-mode";
/*file[PLOT_MODE] = "set_grid_plot_mode_com.h";*/
ARGUMENT[PLOT_MODE]   = 0;
BRIEF_HELP[PLOT_MODE] = "Modify the graphical contents of the "
"`Volume mappings' window.";

COMMAND[HELP]    = "help";
ARGUMENT[HELP]   = 0;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT]    = "exit";
ARGUMENT[EXIT]   = 0;
BRIEF_HELP[EXIT] = "Stop modifying the parameters and exit to the previous "
"command level.";


