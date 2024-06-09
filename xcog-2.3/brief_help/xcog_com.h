enum {MAKE_MAPPING, CHANGE_MAPPING, COPY_MAPPING, TRANSFORM_MAPPING, RENAME_MAPPING,
      DELETE_MAPPING, COLOR_MAPPING, MARK_MAPPING_BNDRY,
      MAKE_CURVE, CHANGE_CURVE, COPY_CURVE, TRANSFORM_CURVE,
      RENAME_CURVE, DELETE_CURVE, BOUNDARY_CONDITION, CURVE_LABEL, 
      GRID_LINES, GRID_QUALITY, OVERLAP_PARAMETERS, LIST_OF_COMPONENTS, COMPUTE_OVERLAP,
      RE_COMPUTE, SAVE_COMPONENT, SAVE_OVERLAPPING, SHOW, PLOT_MODE, READ_COMMAND,
      START_SAVING, STOP_SAVING, RESET_XCOG, READ_STATE, SAVE_STATE, HELP, TUTORIAL, 
      QUIT};

const int LAST_COM = QUIT, LEVEL=0;
char *COMMAND[35], *BRIEF_HELP[35];
int ARGUMENT[35], SAVE_ON_COPY[35];

COMMAND[MAKE_MAPPING] = "make-mapping";           
/*file[MAKE_MAPPING] = "choose_mapping_com.h"; */
SAVE_ON_COPY[MAKE_MAPPING] = 1;  
ARGUMENT[MAKE_MAPPING] = 1;  
BRIEF_HELP[MAKE_MAPPING] = "Design a new mapping which either can be cartesian, normals "
"out from a previously defined curve, linear interpolation between two previously "
"defined curves, interpolation between discrete data points, or a Theodorsen-Garrick "
"mapping for wing-profiles.";

COMMAND[CHANGE_MAPPING] = "change-mapping";         
SAVE_ON_COPY[CHANGE_MAPPING] = CHANGE_MAPPING;  
ARGUMENT[CHANGE_MAPPING] = 1;
BRIEF_HELP[CHANGE_MAPPING] = "Change the properties of an existing mapping. For "
"example, the "
"shape of the mapping, the number of grid lines, the graphical presentation of the "
"mapping. It is also possible to save the plot of the mapping on a postscript file.";

COMMAND[COPY_MAPPING] = "copy-mapping";           
SAVE_ON_COPY[COPY_MAPPING] = 1;  
ARGUMENT[COPY_MAPPING] = 1;
BRIEF_HELP[COPY_MAPPING] = "Create a new mapping and copy all values from an existing "
"mapping.";

COMMAND[TRANSFORM_MAPPING] = "transform-mapping";      
/*file[TRANSFORM_MAPPING] = "set_transformation_com.h"; */
SAVE_ON_COPY[TRANSFORM_MAPPING] = 1;  
ARGUMENT[TRANSFORM_MAPPING] = 1;
BRIEF_HELP[TRANSFORM_MAPPING] = "Transform a mapping by first rotating it, then "
"translating it, and finally scaling it.";

COMMAND[RENAME_MAPPING] = "rename-mapping";         
SAVE_ON_COPY[RENAME_MAPPING] = 1;  
ARGUMENT[RENAME_MAPPING] = 1;
BRIEF_HELP[RENAME_MAPPING] = "Give a mapping a new name.";

COMMAND[DELETE_MAPPING] = "delete-mapping";         
SAVE_ON_COPY[DELETE_MAPPING] = 1;  
ARGUMENT[DELETE_MAPPING] = 1;
BRIEF_HELP[DELETE_MAPPING] = "Delete an existing mapping and purge all lists and arrays "
"associated with the mapping.";

COMMAND[COLOR_MAPPING] = "color-mapping";         
SAVE_ON_COPY[COLOR_MAPPING] = 1;  
ARGUMENT[COLOR_MAPPING] = 1;
BRIEF_HELP[COLOR_MAPPING] = "Change the color of a mapping.";

COMMAND[MARK_MAPPING_BNDRY] = "mark-mapping-boundary";         
/*file[MARK_MAPPING_BNDRY] = "mark_mapping_bndry_com.h"; */
SAVE_ON_COPY[MARK_MAPPING_BNDRY] = 1;  
ARGUMENT[MARK_MAPPING_BNDRY] = 1;
BRIEF_HELP[MARK_MAPPING_BNDRY] = "Identify external segments of mappings that in part "
"are aligned with the boundary of the computational domain.";

COMMAND[MAKE_CURVE] = "make-curve";             
/*file[MAKE_CURVE] = "choose_curve_com.h"; */
SAVE_ON_COPY[MAKE_CURVE] = 1;  
ARGUMENT[MAKE_CURVE] = 1; 
BRIEF_HELP[MAKE_CURVE] = "Design a new curve which can be a straight line, a circular "
"arc, a cubic spline or a smoothed polygon. It is also possible to save the plot of "
"the mapping on a postscript file.";

COMMAND[CHANGE_CURVE] = "change-curve";           
SAVE_ON_COPY[CHANGE_CURVE] = 1;  
ARGUMENT[CHANGE_CURVE] = 1;
BRIEF_HELP[CHANGE_CURVE] = "Change an existing curve, i.e. change the shape of the curve or "
"the graphical presentation of the curve. It is also possible to save the plot of "
"the mapping on a postscript file.";

COMMAND[COPY_CURVE] = "copy-curve";             
SAVE_ON_COPY[COPY_CURVE] = 1;  
ARGUMENT[COPY_CURVE] = 1;
BRIEF_HELP[COPY_CURVE] = "Create a new curve and copy all values from an existing curve.";

COMMAND[TRANSFORM_CURVE] = "transform-curve";        
/*file[TRANSFORM_CURVE] = "set_curve_trans_com.h"; */
SAVE_ON_COPY[TRANSFORM_CURVE] = 1;  
ARGUMENT[TRANSFORM_CURVE] = 1;
BRIEF_HELP[TRANSFORM_CURVE] = "Transform a curve by first rotating it, then translating it, "
"and finally scaling it.";

COMMAND[RENAME_CURVE] = "rename-curve";          
SAVE_ON_COPY[RENAME_CURVE] = 1; 
ARGUMENT[RENAME_CURVE] = 1;
BRIEF_HELP[RENAME_CURVE] = "Give a curve a new name.";

COMMAND[DELETE_CURVE] = "delete-curve";          
SAVE_ON_COPY[DELETE_CURVE] = 1; 
ARGUMENT[DELETE_CURVE] = 1;
BRIEF_HELP[DELETE_CURVE] = "Delete an existing curve and purge all lists and arrays "
"associated with the curve.";

COMMAND[BOUNDARY_CONDITION] = "boundary-condition";    
/*file[BOUNDARY_CONDITION] = "set_boundary_condition_com.h"; */
SAVE_ON_COPY[BOUNDARY_CONDITION] = 1; 
ARGUMENT[BOUNDARY_CONDITION] = 1;
BRIEF_HELP[BOUNDARY_CONDITION] = "Change the boundary conditions for a mapping. The "
"value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is necessary information "
"for a PDE solver that uses the grid.";

COMMAND[CURVE_LABEL] = "curve-label";           
/*file[CURVE_LABEL] = "set_curve_com.h"; */
SAVE_ON_COPY[CURVE_LABEL] = 1; 
ARGUMENT[CURVE_LABEL] = 1;
BRIEF_HELP[CURVE_LABEL] = "Change the curve label for the boundaries of a mapping. "
"The curve label is used to identify the different parts of the boundary of the "
"computational domain. All sides of all component grids that are aligned with "
"the boundary must be assigning a curve label. Sides that are completely aligned "
"with the boundary should be "
"given a positive curve-label, and sides that are partly on the boundary, and "
"partly inside of another grid, should be given a negative curve-label. The "
"absolute value of the curve-label should be the same for all sides of all "
"grids that are aligned with the same part of the boundary. Hence, if the "
"domain is simply connected, only one curve-label value should be used. For "
"a doubly connected domain, two values should be used, and so on.";

COMMAND[GRID_LINES] = "grid-lines";            
/*file[GRID_LINES] = "set_gridlines_com.h"; */
SAVE_ON_COPY[GRID_LINES] = 1; 
ARGUMENT[GRID_LINES] = 1;
BRIEF_HELP[GRID_LINES] = "Set the number of grid lines for a mapping.";

COMMAND[GRID_QUALITY] = "grid-quality";          
SAVE_ON_COPY[GRID_QUALITY] = 1; 
ARGUMENT[GRID_QUALITY] = 1;
BRIEF_HELP[GRID_QUALITY] = "Evaluate a number of measures of the quality of a grid "
"associated with a mapping.";

COMMAND[OVERLAP_PARAMETERS] = "overlap-parameters";    
/*file[OVERLAP_PARAMETERS] = "overlap_parameters_com.h"; */
SAVE_ON_COPY[OVERLAP_PARAMETERS] = 1; 
ARGUMENT[OVERLAP_PARAMETERS] = 0;
BRIEF_HELP[OVERLAP_PARAMETERS] = "Change the overlap parameters, i.e. the discretization width, the "
"interpolation width, etc.";

COMMAND[LIST_OF_COMPONENTS] = "list-of-components";    
/*file[LIST_OF_COMPONENTS] = "component_grid_list_com.h"; */
SAVE_ON_COPY[LIST_OF_COMPONENTS] = 1; 
ARGUMENT[LIST_OF_COMPONENTS] = 0;
BRIEF_HELP[LIST_OF_COMPONENTS] = "Change or reorder the list of mappings which will be included "
"in the overlapping grid.";

COMMAND[COMPUTE_OVERLAP] = "compute-overlap";       
/*file[LIST_OF_COMPONENTS] = "compute_overlap_com.h"; */
SAVE_ON_COPY[COMPUTE_OVERLAP] = 1; 
ARGUMENT[COMPUTE_OVERLAP] = 0;
BRIEF_HELP[COMPUTE_OVERLAP] = "Set the parameters for the composite grid and compute the "
"intergrid communication from a selected subset of component grids.\n"
"If the algorithm succeeds, the interpolation points will be marked by CIRCLES.\n"
"If the algorithm fails, the grid points that didn't satisfy the requirements for "
"interpolation point or discretization point are marked in the following way:\n"
"A BOX indicates a bad discretization point.\n"
"A CROSS is a dead interpolation location.\n"
"A FILLED CIRCLE represents a non-explicit interpolation location.\n"
"A TRIANGLE marks a bad interpolation point because of a dead interpolation "
"location.\n"
"A DOWNWARD TRIANGLE indicates a bad interpolation point because of a non explicit "
"interpolation location.";

COMMAND[RE_COMPUTE] = "re-compute-overlap";    
SAVE_ON_COPY[RE_COMPUTE] = 1; 
ARGUMENT[RE_COMPUTE] = 0;
BRIEF_HELP[RE_COMPUTE] = "Re-apply the overlapping grid algorithm to the already existing "
"overlapping grid. This might be necessary to remove orphan points when many grids "
"overlap each other in a complicated way.";

COMMAND[SAVE_COMPONENT] = "save-component-grid";   
SAVE_ON_COPY[SAVE_COMPONENT] = 1; 
ARGUMENT[SAVE_COMPONENT] = 1;
BRIEF_HELP[SAVE_COMPONENT] = "Save one component grid on an ASCII file in PLOT3D format.";

COMMAND[SAVE_OVERLAPPING] = "save-overlapping-grid"; 
SAVE_ON_COPY[SAVE_OVERLAPPING] = 1; 
ARGUMENT[SAVE_OVERLAPPING] = 0;
BRIEF_HELP[SAVE_OVERLAPPING] = "Save the overlapping grid on an ASCII or HDF file. "
"Optionally, the Jacobian of "
"the transformation can be saved. The Jacobian is computed by calling the "
"forward mapping function, if it is available. Otherwise, the Jacobian is "
"approximated by second order accurate divided differences.";

COMMAND[SHOW] = "show-curves-mappings";  
SAVE_ON_COPY[SHOW] = 1; 
ARGUMENT[SHOW] = 0;
BRIEF_HELP[SHOW] = "Display the lists of all curves and mappings presently defined in "
"the program.";

COMMAND[PLOT_MODE] = "composite-plot-mode";             
/*file[PLOT_MODE] = "set_global_plot_com.h"; */
SAVE_ON_COPY[PLOT_MODE] = 1; 
ARGUMENT[PLOT_MODE] = 0;
BRIEF_HELP[PLOT_MODE] = "Modify the graphical contents of the `Composite grid' window "
"or save the graphics on a postscript file.";

COMMAND[READ_COMMAND] = "read-command-file";     
SAVE_ON_COPY[READ_COMMAND] = 0; 
ARGUMENT[READ_COMMAND] = 1;
BRIEF_HELP[READ_COMMAND] = "Start reading commands from a file instead from the keyboard.";

COMMAND[START_SAVING] = "start-saving-commands"; 
SAVE_ON_COPY[START_SAVING] = 1; 
ARGUMENT[START_SAVING] = 1;
BRIEF_HELP[START_SAVING] = "Start saving a copy of all commands on a file that later can be "
"read by `read-command-file'.";

COMMAND[STOP_SAVING] = "stop-saving-commands";  
SAVE_ON_COPY[STOP_SAVING] = 0; 
ARGUMENT[STOP_SAVING] = 1;
BRIEF_HELP[STOP_SAVING] = "Stop saving commands and close the command file that was opened "
"by `start-saving-commands'.";

COMMAND[RESET_XCOG] = "reset-xcog";            
SAVE_ON_COPY[RESET_XCOG] = 1; 
ARGUMENT[RESET_XCOG] = 0;
BRIEF_HELP[RESET_XCOG] = "Clear the lists of curves and mappings, purge any composite grid "
"and reset the graphics.";

COMMAND[READ_STATE] = "read-state";            
SAVE_ON_COPY[READ_STATE] = 1; 
ARGUMENT[READ_STATE] = 1;
BRIEF_HELP[READ_STATE] = "Recreate curves and mappings from a HDF file that previously "
"was saved with the command `save-curves-mappings'.";

COMMAND[SAVE_STATE] = "save-state";            
SAVE_ON_COPY[SAVE_STATE] = 1; 
ARGUMENT[SAVE_STATE] = 1;
BRIEF_HELP[SAVE_STATE] = "Save all curves and mappings on a HDF file that later can be "
"read with the command `read-curves-mappings'.";

COMMAND[HELP] = "help";                  
SAVE_ON_COPY[HELP] = 1; 
ARGUMENT[HELP] = 0;
BRIEF_HELP[HELP] = NULL;

COMMAND[TUTORIAL] = "tutorial";              
/*file[TUTORIAL] = "run_tutorial_com.h"; */
SAVE_ON_COPY[TUTORIAL] = 1; 
ARGUMENT[TUTORIAL] = 0;
BRIEF_HELP[TUTORIAL] = "Introduction to Xcog and composite overlapping grids.";

COMMAND[QUIT] = "quit";                  
SAVE_ON_COPY[QUIT] = 1; 
ARGUMENT[QUIT] = 0;
BRIEF_HELP[QUIT] = "Terminate the Xcog session. To use the results of "
"this session outside of Xcog, you must do `save-overlapping-grid' "
"before you quit. If you want to use the currently defined curves and "
"mappings in a future Xcog session, you should also do `save-state'.";
