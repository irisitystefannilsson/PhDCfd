enum {CARTESIAN, NORMAL_CURVE, LINEAR_INTERP, DISC_POINT,
#ifdef THE_GAR
      TEO_GAR,
#endif
      HYPERBOLIC, HELP, CANCEL};

const int LAST_COM = CANCEL, LEVEL=1;
#ifdef THE_GAR
char *COMMAND[8], *BRIEF_HELP[8];
int ARGUMENT[8], *SAVE_ON_COPY=NULL;
#else
char *COMMAND[7], *BRIEF_HELP[7];
int ARGUMENT[7], *SAVE_ON_COPY=NULL;
#endif

COMMAND[CARTESIAN] = "cartesian-mapping";      
/*file[CARTESIAN] = "cartesian_grid_com.h"; */
ARGUMENT[CARTESIAN] = 0;   
BRIEF_HELP[CARTESIAN] = "Define the mapping to be Cartesian and proceed to modify the "
"default parameters of the mapping.";

COMMAND[NORMAL_CURVE] = "normal-curve-mapping";  
/*file[NORMAL_CURVE] = "normal_curve_com.h"; */
ARGUMENT[NORMAL_CURVE] = 1;   
BRIEF_HELP[NORMAL_CURVE] = "Define the mapping to be normals sticking out a constant or "
"varying distance from an existing curve.";

COMMAND[LINEAR_INTERP] = "linear-interpolation-mapping"; 
/*file[LINEAR_INTERP] = "linear_interp_grid_com.h"; */
ARGUMENT[LINEAR_INTERP] = 0;   
BRIEF_HELP[LINEAR_INTERP] = "Define the mapping to be linear interpolation between two "
"different existing curves.";

COMMAND[DISC_POINT] = "discrete-point-mapping"; 
/*file[DISC_POINT] = "init_disc_point_com.h"; */
ARGUMENT[DISC_POINT] = 0;   
BRIEF_HELP[DISC_POINT] = "Define the mapping to interpolate inbetween a discrete set of "
"grid points, which are read from a file that for instance could be the output from "
"an external grid generator.";

COMMAND[HYPERBOLIC] = "hyperbolic-mapping";  
/*file[HYPERBOLIC] = "set_hyp_com.h";*/
ARGUMENT[HYPERBOLIC] = 1;
BRIEF_HELP[HYPERBOLIC] = "Define the mapping by solving a Hamilton-Jacobi equation "
"in the direction normal to a curve. (The mapping is called hyperbolic since it "
"belongs to the class of grid generation techniques that commonly is known as "
"hyperbolic). This mapping does NOT have a known analytical inverse. The mapping "
"function is\n"
"\n"
"d/dt (x, y) = ds/dt * v(kappa) * unit_normal\n"
"\n"
"where t and s are the parameter and arclength, respectively, along a grid line "
"normal to the surface. kappa is the mean curvature of the surface and\n"
"\n"
"v(kappa) = max( v_min, 1 - eps * kappa )\n"
"\n"
"is the speed function. Here v_min is the velocity threshold that dictates the "
"propagation speed close to highly curved convex surfaces. "
"The quantity v(kappa) * unit_normal, can be averaged to make the grid smoother. "
"In this case the algebraic average of the nearest neighbors is computed and "
"weighted together with the point value.";

COMMAND[HELP] = "help";
ARGUMENT[HELP] = 0;   
BRIEF_HELP[HELP] = NULL;

COMMAND[CANCEL] = "cancel";
ARGUMENT[CANCEL] = 0;   
BRIEF_HELP[CANCEL] = "Don't define a mapping. Instead go back to the previous command "
"level.";

#ifdef THE_GAR
COMMAND[TEO_GAR] = "theodorsen-garrick-mapping"; 
/*file[TEO_GAR] = "tg_grid_com.h"; */
ARGUMENT[TEO_GAR] = 1;   
BRIEF_HELP[TEO_GAR] = "Define a Theodorsen-Garrick mapping. This conformal mapping is "
"suitable for making grids around wing sections. You will be prompted for a "
"coordinate file where the wing "
"section is defined. If you don't give a valid file-name, you will get the default "
"NACA-0012 profile. The coordinate file should have the same format as a node point "
"file for a cubic spline, i.e.:\n"
"#xcog-file-type-4\n"
"5\n"
"   1.00    0.00\n"
"   0.50    0.05\n"
"   0.00    0.00\n"
"   0.50   -0.05\n"
"   1.00    0.00\n"
"The first line of the file contains a descriptor which tells xcog that it is a node "
"point file. The second line contains the number of node points, which must be "
"greater than or equal to 3. On the following lines are the x and y coordinates for "
"the node points. The number of x and y coordinates must be equal to the "
"given number of node points, and every consecutive pair of node points must be "
"unique. All coordinates must satisfy 0 <= x <= 1, and the first and last "
"coordinates must have x=1, y=0. Furthermore, there must be exactly one coordinate "
"at x=0.\n";
#endif
