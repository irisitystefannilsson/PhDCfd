#define ncom 5
  char *command[ncom+1], *brief_help[ncom+1];
  int *argument=NULL, *save_on_copy=NULL;

command[0] = "straight-line";
/*file[0] = "straight_line_com.h"; */
brief_help[0] = "Define curve to be a straight line. Proceed to modify the default "
"parameters of the curve.";

command[1] = "circular-arc";      
/*file[1] = "circular_arc_com.h"; */
brief_help[1] = "Define curve to be a circular arc. Proceed to modify the default "
"parameters of the curve.";

command[2] = "cubic-spline";        
/*file[2] = "cubic_spline_com.h"; */
brief_help[2] = "Define curve to be a cubic spline. Proceed to modify the default "
"parameters of the curve.";

command[3] = "smooth-polygon";   
/*file[3] = "smooth_poly_curve_com.h"; */
brief_help[3] = "Define curve to be a smooth polygon. Proceed to modify the default "
"parameters of the curve.";

command[4] = "help";                    
brief_help[4] = NULL;

command[5] = "cancel";                    
brief_help[5] = "Don't define a curve. Instead go back to the previous command "
"level.";

