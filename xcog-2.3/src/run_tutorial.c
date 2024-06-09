#include "xcog.h"

void 
run_tutorial( input_output *io_ptr ){
  FILE *fp;
  char *xcog_home;
  int quit, icom;

#include "run_tutorial_com.h"

/* always present some general hints and a brief outline of the principles */
  printf(
"Welcome to the tutorial of Xcog version "VERSION".\n"
"\n"
"A user's guide to the program can be found at\n"
"http://www.na.chalmers.se/~andersp/xcog/xcog.html\n"
"\n");
  general_help();

  do{

    icom = get_command( io_ptr, "tutorial>", COMMAND, LAST_COM, 
		       LEVEL, SAVE_ON_COPY, ARGUMENT);

    quit = 1;
    switch( icom ){

    case STEPS:
/* overlapping grid steps */
      tutorial_intro();
      quit = 0;
      break;

    case C_INTERP:
/* command interpreter */
      general_help();
      quit = 0;
      break;

    case CYL_IN_SQ:
/* cylinder in square */
      if ((fp = open_this_ascii_file( "cylinder-square.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
      }
      else
	quit = 0;
      break;

    case PAPER_KNIFE:
/* paper knife */
      if ((fp = open_this_ascii_file( "paper-knife.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
      }
      else
	quit = 0;
      break;

    case SHOCK_CYL:
/* cylinder shock */
      if ((fp = open_this_ascii_file( "shock-cylinder.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
      }
      else
	quit = 0;
      break;

    case INT_FLOW:
/* internal flow */
      if ((fp = open_this_ascii_file( "internal-flow.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
      }
      else
	quit = 0;
      break;

    case MEX_GULF:
/* read externally defined grid on plot3d format */
      if ((fp = open_this_ascii_file( "mex-gulf.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
      }
      else
	quit = 0;
      break;

    case HYP_WING:
      if ((fp = open_this_ascii_file( "hyp-wing.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
      }
      else
	quit = 0;
      break;

    case C_GRID:
      if ((fp = open_this_ascii_file( "c-grid.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
      }
      else
	quit = 0;
      break;

#ifdef THE_GAR
    case FOIL_BOX:
      if ((fp = open_this_ascii_file( "foil-box.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
      }
      else
	quit = 0;
      break;
#endif

    case MISC:
/* misc info: postscript, read saved state, change appearance in plot window, etc. */
      if ((fp = open_this_ascii_file( "misc.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
      }
      else
	quit = 0;
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );

      quit = 0;
      break;

    case EXIT:
/* exit tutorial */
      break;

    default:
      quit = 0;
    }
  } while (!quit);

}

void 
tutorial_intro( void ){
  printf(
"Go through the following steps to construct an overlapping grid with Xcog:\n\n"
" 1) Generate curves that describe the geometry.\n"
" 2) Create grids (mappings) that discretize the domain close to each curve.\n"
" 3) Fill the remaining part of the domain with (Cartesian) background grids.\n"
" 4) Identify the different parts of the boundary of the computational domain\n"
"    by assigning a curve_label to the sides of the grids that are aligned with\n"
"    the boundary. Sides that are completely aligned with the boundary should be\n"
"    given a positive curve-label, and sides that are partly on the boundary, and\n"
"    partly inside of another grid, should be given a negative curve-label. The\n"
"    absolute value of the curve-label should be the same for all sides of all\n"
"    grids that are aligned with the same part of the boundary. Hence, if the\n"
"    domain is simply connected, only one curve_label value should be used. For\n"
"    a doubly connected domain, two values should be used, and so on.\n"
" 5) To make the overlapping grid useful for a PDE solver, assign boundary\n"
"    condition values to the sides of the components that are aligned with the\n"
"    physical boundary.\n"
" 6) Order the different mappings in the grid hierarcy. The first grid will\n"
"    be given the highest priority in the overlapping grid.\n"
" 7) Set the discretization and interpolation widths to suit the discretization\n"
"    of the partial differential equation.\n"
" 8) Apply the intergrid communication algorithm to compute the overlapping grid.\n"
" 9) Save the overlapping grid.\n\n");
}
