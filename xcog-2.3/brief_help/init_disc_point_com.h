#define ncom 2
char *command[ncom+1], *brief_help[ncom+1];
int argument[ncom+1], *save_on_copy=NULL;

command[0] = "plot3d-ascii"; 
argument[0] = 1;
/*file[DISC_POINT] = "disc_point_map_com.h"; */
brief_help[0] = "Read the grid points from a PLOT3D ascii file with the "
"following format:\n" 
" nii njj nkk\n"
" x(1,1)  x(2,1)  ... x(nr,1)\n"
" x(1,2)  x(2,2)  ... x(nr,2)\n"
" ...\n"
" x(1,ns) x(2,ns) ... x(nr,ns)\n"
" y(1,1)  y(2,1)  ... y(nr,1)\n"
" y(1,2)  y(2,2)  ... y(nr,2)\n"
" ...\n"
" y(1,ns) y(2,ns) ... y(nr,ns)\n"
"Because xcog only can handle two-dimensional grids, it is required that at least one of nii, njj and nkk equal ONE. The first non-1 of nii and njj is set to nr and the remaining non-1 is set to ns.\n"
"Example: nii = 3, njj = 1, nkk = 4 yields nr = 3 and ns = 4.";

command[1] = "help";
argument[1] = 0;
brief_help[1] = NULL;

command[2] = "cancel";
argument[2] = 0;
brief_help[2] = "Exit from this command level without reading a grid point file.";
