#define ncom 5
char *command[ncom+1], *brief_help[ncom+1];
int argument[ncom+1], *save_on_copy=NULL;

for (i=0; i<=ncom; i++) argument[i]=0;

command[0] = "arclength-stretching";      
command[1] = "exponential-stretching";
/*file[1] = "set_exp_com.h"; */
command[2] = "layer-stretching"; 
/*file[2] = "set_layer_com.h"; */

command[3] = "hyperbolic-tangent-stretching"; 
/*file[3] = "set_tanh_com.h"; */

command[4] = "help";
command[5] = "cancel";

brief_help[0] = "Reparametrize the curve with respect to arclength. This will only "
"change the distribution of grid points along the curve, and not the shape of the "
"curve itself.";

brief_help[1] = "Introduce an exponential stretching. This stretching is mostly used "
"to concentrate grid points close to the starting or ending point of the curve.";

brief_help[2] = "Attract grid points to the node points of the curve. This "
"stretching is mostly used for smoothed polygon curves.";

brief_help[3] = "Introduce a hyperbolic tangent stretching function. This function "
"can be used to attract grid points to both ends of the curve.";

brief_help[4] = NULL;

brief_help[5] = "Don't pick a stretching. Instead go back to the previous command "
"level.";

