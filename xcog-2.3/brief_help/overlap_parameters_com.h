enum{CORNER_WIDTH, DISC_WIDTH, GHOSTPOINTS, INTERP_TYPE, INTERP_WIDTH, NORMAL_WIDTH,
     TANG_WIDTH, ALL_WIDTHS, TOLERANCE, CORRECT_MISMATCH, TRIM_STYLE, PHYSICAL_BNDRY, 
     SHOW_HOLES, SHOW, HELP, EXIT};

const int LAST_COM=EXIT;
int argument[16], *save_on_copy=NULL;
char *command[16], *brief_help[16];

command[SHOW_HOLES] ="show-holes";
argument[SHOW_HOLES] = 1;
brief_help[SHOW_HOLES] = "Toggle the flag that determines if the "
"used grid points in the overlapping grid should be shown after the "
"hole-cutting has been applied. If the overlap algorithm fails, it is often "
"possible to determine if it failed because of an inconsistently defined boundary, "
"or because of too few grid points in the overlap region, "
"by looking at the grid after the holes have been cut. If holes have been cut in "
"the wrong places, it is probably because some sides of some grids have been "
"given an inconsistent curve-label. But if the overlap algorithm fails "
"eventhough the hole are cut correctly, the problem is most likely "
"caused by too few grid points in the overlap region.";

command[PHYSICAL_BNDRY] ="show-physical-boundaries";
argument[PHYSICAL_BNDRY] = 1;
brief_help[PHYSICAL_BNDRY] = "Toggle the flag that determines if the global "
"boundary of the computational domain will be shown during the overlap algorithm.";

command[CORNER_WIDTH] ="corner-width";         
argument[CORNER_WIDTH] = 1;
brief_help[CORNER_WIDTH] = "Change the corner width. The corner width equals the "
"width of "
"the discretization stencil at a corner point in the grid and at points so close "
"to the corner that the regular discretization stencil can't be used.";

command[DISC_WIDTH] ="discretization-width"; 
argument[DISC_WIDTH] = 1;
brief_help[DISC_WIDTH] = "Change the discretization width, which is defined as the "
"maximum number of grid points in any coordinate direction of the discretization "
"stencil away from boundaries.";

command[GHOSTPOINTS] ="ghostpoints";          
argument[GHOSTPOINTS] = 1;
brief_help[GHOSTPOINTS] = "Change the number of extra grid points outside the "
"boundary of "
"the grid. These points have no influence on the composite grid algorithm, but are "
"sometimes useful for discretizing a PDE close to a boundary.";

command[INTERP_TYPE] ="interpolation-type";   
argument[INTERP_TYPE] = 0;
brief_help[INTERP_TYPE] = "Change the interpolation type, which either can be "
"implicit or explicit. When the interpolation is explicit, the overlap between "
"the grids is "
"sufficiently wide to ensure that no interpolation point is present in the "
"interpolation formula for any interpolation point. When the interpolation "
"is implicit, this restriction is removed, so the overlap can usually be narrower.";

command[INTERP_WIDTH] ="interpolation-width";  
argument[INTERP_WIDTH] = 1;
brief_help[INTERP_WIDTH] = "Change the width of the interpolation formula, which "
"is defined as the maximum number number of grid points in any coordinate direction "
"of the interpolation stencil.";

command[NORMAL_WIDTH] ="normal-width";         
argument[NORMAL_WIDTH] = 1;
brief_help[NORMAL_WIDTH] = "Change the width of the discretization stencil at the "
"boundary, normal to the boundary.";

command[TANG_WIDTH] ="tangential-width";     
argument[TANG_WIDTH] = 1;
brief_help[TANG_WIDTH] = "Change the width of the discretization stencil at the "
"boundary, tangent to the boundary.";

command[ALL_WIDTHS] ="set-all-widths";       
argument[ALL_WIDTHS] = 1;
brief_help[ALL_WIDTHS] = "Set all widths at once, i.e. the discretization width, the "
"normal width, the tangential width, the corner width and the interpolation width.";

command[TOLERANCE] ="extra-curve-tolerance";      
argument[TOLERANCE] = 1;           
brief_help[TOLERANCE] = "Set the extra tolerance for curve missmatch at physical "
"boundaries. The basic mismatch tolerance is estimated by the code by looking at "
"the smoothness of the polygon between the boundary grid points which represent "
"the global boundary. The extra tolerance "
"is added to the basic value to yield the total mismatch tolerance, which is used "
"together with matching curve values to help the intergrid communication algorithm "
"find interpolation points close to boundaries.";

command[CORRECT_MISMATCH] ="correct-mismatch";       
argument[CORRECT_MISMATCH] = 1;
brief_help[CORRECT_MISMATCH] = "Set the flag that determines if the interpolation "
"weights and locations should be corrected to account for boundary mismatch. The "
"mismatch occurs when a boundary is represented by the sides of several components "
"when the corresponding mappings have a non-unique representation of the "
"boundary. The mismatch problem becomes more pronounced when the grids are very "
"fine in the direction normal to the boundary, and the grid cells have a large "
"aspect ratio.";

command[TRIM_STYLE] = "trimming-style";     
argument[TRIM_STYLE] = 0;
brief_help[TRIM_STYLE] = "Toggle the trimming style between minimizing and "
"maximizing the overlap in the overlapping grid.";

command[SHOW] ="show";                
argument[SHOW] = 0;
brief_help[SHOW] = "Show the present value of the overlap parameters.";

command[HELP] ="help";                
argument[HELP] = 0;
brief_help[HELP] = NULL;

command[EXIT] ="exit";                
argument[EXIT] = 0;
brief_help[EXIT] = "Stop changing the overlap parameters and exit from this command "
"level.";
