#include "xcog.h"

static void 
set_interpolation_type(input_output *io_ptr, xcog_data *xcog_data_ptr);

void 
overlap_parameters(input_output *io_ptr, xcog_data *xcog_data_ptr){
  int level=2, no_indent=0;
  char *prompt;
  int icom, quit, all_widths;
  real new_max_n_dist;

#include "overlap_parameters_com.h"

  prompt = "overlap parameters>";
  quit = 0;
  do{
    icom = get_command(io_ptr, prompt, command, LAST_COM, level,
		       save_on_copy, argument);

    switch (icom) {

    case PHYSICAL_BNDRY:
      xcog_data_ptr->show_physical_boundary = 
	get_yes_no(io_ptr, "Display the global boundary during the "
		   "overlap algorithm? ", level);
      break;

    case SHOW_HOLES:
      xcog_data_ptr->show_holes = 
	get_yes_no(io_ptr, "Display the result of the hole-cutting during the "
		   "overlap algorithm? ", level);
      break;

    case CORRECT_MISMATCH:
      xcog_data_ptr->correct_mismatch = 
	get_yes_no(io_ptr, "Correct interpolation locations for boundary mismatch? ",
		   level);
      break;

    case CORNER_WIDTH:
      xcog_data_ptr->corner_width = 
	int_max(1, get_int(io_ptr, "corner width > 0: ", 
			   xcog_data_ptr->corner_width, no_indent));
      break;

    case DISC_WIDTH:
      xcog_data_ptr->disc_width = 
	int_max(1, get_int(io_ptr, 
			   "discretization width > 0: ", 
			   xcog_data_ptr->disc_width, no_indent));
      break;

    case GHOSTPOINTS:
      xcog_data_ptr->extra = 
	int_max(0, get_int(io_ptr, 
			   "number of ghostpoints >= 0: ", 
			   xcog_data_ptr->extra, no_indent));
      break;

    case INTERP_TYPE:
      set_interpolation_type(io_ptr, xcog_data_ptr);
      break;

    case INTERP_WIDTH:
      xcog_data_ptr->interp_width = 
	int_max(1, get_int(io_ptr, 
			   "interpolation width > 0: ", 
			   xcog_data_ptr->interp_width, no_indent));
      break;

    case NORMAL_WIDTH:
      xcog_data_ptr->normal_width = 
	int_max(1, get_int( io_ptr, "normal width > 0: ", 
			   xcog_data_ptr->normal_width, no_indent));
      break;

    case TANG_WIDTH:
      xcog_data_ptr->tangent_width = 
	int_max(1, get_int(io_ptr, "tangential width > 0: ", 
			   xcog_data_ptr->tangent_width, no_indent));
      break;

    case ALL_WIDTHS:
      all_widths = (xcog_data_ptr->disc_width + xcog_data_ptr->interp_width + 
		    xcog_data_ptr->normal_width + xcog_data_ptr->tangent_width +
		    xcog_data_ptr->corner_width)/5.0;
      all_widths = 
	int_max(1, get_int(io_ptr, "all widths > 0: ", all_widths, no_indent));
/* assign */
      xcog_data_ptr->disc_width    =  all_widths;
      xcog_data_ptr->tangent_width =  all_widths;
      xcog_data_ptr->normal_width  =  all_widths;
      xcog_data_ptr->corner_width  =  all_widths;
      xcog_data_ptr->interp_width  =  all_widths;

      break;

    case TOLERANCE:
      do {
	new_max_n_dist = get_real(io_ptr, "Extra curve tolerance >= 0: ", 
				  xcog_data_ptr->max_n_dist, no_indent);
	if ( new_max_n_dist < 0 )
	  printf("Input error: the extra curve tolerance must be non-negative!\n");
      }
      while ( new_max_n_dist < 0 );
      xcog_data_ptr->max_n_dist = new_max_n_dist;
      break;

    case TRIM_STYLE:
/* toggle trimming style */
      xcog_data_ptr->trim_style = !xcog_data_ptr->trim_style;
      printf("The trimming style is to %s the overlap.\n", 
	     (xcog_data_ptr->trim_style)? "maximize": "minimize");
      break;

    case SHOW:
      show_overlap_param( xcog_data_ptr );
      break;

    case HELP:
      while ( (icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				  command, LAST_COM, level+1, NULL, NULL)) == -1 );
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case EXIT:
      quit = 1;
      break;

    default:
      ;
    }
  }
  while(!quit);
}


void 
show_overlap_param( xcog_data *xcog_data_ptr ){
  printf("\n");
  printf("corner width: %i\n",xcog_data_ptr->corner_width);
  printf("discretization width: %i\n",xcog_data_ptr->disc_width);
  printf("ghostpoints: %i\n",xcog_data_ptr->extra);
  printf("interpolation width: %i\n",xcog_data_ptr->interp_width);
  printf("normal width: %i\n",xcog_data_ptr->normal_width);
  printf("tangential width: %i\n",xcog_data_ptr->tangent_width);
  
  printf("\n");
  printf("interpolation type is ");
  if (xcog_data_ptr->interp_type == 'e')
    printf("explicit\n");
  else if (xcog_data_ptr->interp_type == 'i')
    printf("implicit\n");
  else
    printf("unknown\n");

  printf("\n");
  printf("extra curve tolerance: %e\n", xcog_data_ptr->max_n_dist);

  printf("\n");
  printf("The global boundary %s be shown during the overlap algorithm.\n", 
	 (xcog_data_ptr->show_physical_boundary? "WILL" : "WILL NOT") );
  printf("The used grid points %s be shown after the holes have been cut, \n"
	 "during the overlap algorithm.\n", 
	 (xcog_data_ptr->show_holes? "WILL" : "WILL NOT") );
  printf("The interpolation locations %s be corrected for boundary mismatch.\n", 
	 (xcog_data_ptr->correct_mismatch? "WILL" : "WILL NOT") );
  printf("\n");
  printf("The trimming style is to %s the overlap.\n", 
	 (xcog_data_ptr->trim_style)? "MAXIMIZE": "MINIMIZE");
}


static void 
set_interpolation_type(input_output *io_ptr, xcog_data *xcog_data_ptr){
  char *prompt;
  int quit, icom;
  int level=0;

#include "set_interp_type_com.h"

  prompt = "set interpolation type>";


  do {
    icom = get_command(io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    quit = 0;
    switch (icom) {

    case 0: 
      xcog_data_ptr->interp_type = 'e';
      quit = 1;
      break;

    case 1:
      xcog_data_ptr->interp_type = 'i';
      quit = 1;
      break;

    case 2:
      printf("The interpolation type is ");
      if (xcog_data_ptr->interp_type == 'e')
	printf("explicit.\n");
      else if (xcog_data_ptr->interp_type == 'i')
	printf("implicit.\n");
      else
	printf("unknown.\n");
      
      printf("\n");
      break;

    case 3:
/* help */
      while ( (icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				  command, ncom, level+1, NULL, NULL)) == -1 );
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );

      break;

    case 4:
      quit = 1;
      break;

    default:
      ;
    }
  }
  while (!quit);
}


