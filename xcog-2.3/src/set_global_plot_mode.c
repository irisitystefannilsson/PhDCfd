#include "xcog.h"

void set_global_plot_mode(input_output *io_ptr, xcog_data *xcog_data_ptr,
			  overlapping_grid *over_ptr){

  char *prompt, *token, *file_name, curve_string[25];
  int icom, quit, replot, col_flag;
  const int save_command=1;
  generic_mapping *grid_ptr;
  generic_curve *curve_ptr;
  static int coordinate_axis=1;
  FILE *fp;
  boundary_part *boundary_ptr=NULL;
  int boundary_plot_mode = 1|2, curve_label;

#include "set_global_plot_com.h"

  prompt = "global plot-mode>";
  quit = 0;
  replot = FALSE; /* assume that the current plot already is plotted */
  do{

    if (replot){
/* use the main window */
      PL_erase(0);
      PL_start_plot(0);
/* plot the curves */
      for (curve_ptr = xcog_data_ptr->first_c; curve_ptr != NULL;
	   curve_ptr = curve_ptr->next)
	plot_curve( curve_ptr, xcog_data_ptr->curve_plot_mode );
/* plot the mappings */
      for (grid_ptr = xcog_data_ptr->first_m; grid_ptr != NULL;
	   grid_ptr = grid_ptr->next_m)
	plot_generic_mapping( grid_ptr, xcog_data_ptr->grid_plot_mode );
/* plot the composite grid */
      plot_overlapping( over_ptr, xcog_data_ptr->comp_plot_mode );
      plot_global_boundary( boundary_ptr, boundary_plot_mode );
      PL_end_plot();
      replot = FALSE;

    }

    icom = get_command( io_ptr, prompt, command, LAST_COM, 
		       LEVEL, save_on_copy, argument);

    switch (icom) {
    case CURVES:
/* curve */
      if (!(xcog_data_ptr->curve_plot_mode & 1))
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
      else
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode & ~1;
      replot = TRUE;
      break;

    case CURVE_DIR:
/* curve direction */
      if (!(xcog_data_ptr->curve_plot_mode & 2))
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 2;
      else
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode & ~2;
      replot = TRUE;
      break;

    case CURVE_TICK:
/* curve tickmarks*/
      if (!(xcog_data_ptr->curve_plot_mode & 4))
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 4;
      else
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode & ~4;
      replot = TRUE;
      break;

    case CURVE_CONTROL_PTS:
/* curve control points*/
      if (!(xcog_data_ptr->curve_plot_mode & 16))
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 16;
      else
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode & ~16;
      replot = TRUE;
      break;

    case GRID_BOUNDARY:
/* grid boundary */
      if (!(xcog_data_ptr->grid_plot_mode & 1))
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
      else
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode & ~1;
      replot = TRUE;
      break;

    case GRID_ARROWS:
/* grid directional arrows */
      if (!(xcog_data_ptr->grid_plot_mode & 2))
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 2;
      else
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode & ~2;
      replot = TRUE;
      break;

    case GRID_TICKS:
/* grid tickmarks */
      if (!(xcog_data_ptr->grid_plot_mode & 4))
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 4;
      else
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode & ~4;
      replot = TRUE;
      break;

    case ALL_GRID_LINES:
/* all gridlines */
      if (!(xcog_data_ptr->grid_plot_mode & 8))
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 8;
      else
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode & ~8;
      replot = TRUE;
      break;

    case BOUNDARY_COND:
/* boundary condition */
      if (!(xcog_data_ptr->grid_plot_mode & 16))
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 16;
      else
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode & ~16;
      replot = TRUE;
      break;

    case CURVE_LABEL:
/* curve value */
      if (!(xcog_data_ptr->grid_plot_mode & 32))
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 32;
      else
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode & ~32;
      replot = TRUE;
      break;

    case OVERLAPPING_GRID:
/* overlapping grid */
      if (!(xcog_data_ptr->comp_plot_mode & 1))
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode | 1;
      else
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode & ~1;
      replot = TRUE;
      break;

    case INTERP_PTS:
/* interpolation points */
      if (!(xcog_data_ptr->comp_plot_mode & 2))
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode | 2;
      else
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode & ~2;
      replot = TRUE;
      break;

    case BAD_PTS:
/* bad points */
      if (!(xcog_data_ptr->comp_plot_mode & 4))
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode | 4;
      else
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode & ~4;
      replot = TRUE;
      break;

    case PHYS_BNDRY:
      if (over_ptr){
	curve_label = (boundary_ptr)? boundary_ptr->curve_label: 1;
	curve_label = get_int( io_ptr, "Which curve-label> ", curve_label, NO_INDENT );
	if (!(boundary_ptr = find_curve( over_ptr->boundary_list, curve_label ))){
	  printf("Sorry, there is no physical boundary corresponding to "
		 "curve_label = %i.\n", curve_label);
	} 
	else{
	  replot = TRUE;
	}
      }
      else{
	printf("Sorry, there is no overlapping grid to choose a physical "
	       "boundary from.\n");
      }
      break;

    case NO_PHYS_BNDRY:
      if (boundary_ptr)
	replot = TRUE;

      boundary_ptr = NULL;
      break;

    case FLAG:
/* flag array */
      if (!(xcog_data_ptr->comp_plot_mode & 8))
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode | 8;
      else
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode & ~8;
      replot = TRUE;
      break;

    case TITLE:
      if (!(xcog_data_ptr->comp_plot_mode & 16)){
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode | 16;
	token = get_word( io_ptr, "New title: ", "My_xcog_data", 1 );
	if (over_ptr->title != NULL)
	  free( over_ptr->title );
	over_ptr->title = (char *) malloc( (strlen(token)+1) * sizeof(char) );
	over_ptr->title = strcpy(over_ptr->title, token);
      }
      else{
	xcog_data_ptr->comp_plot_mode = xcog_data_ptr->comp_plot_mode & ~16;
      }
      replot = TRUE;
      break;

    case C_AXIS:
      coordinate_axis = !coordinate_axis;
      PL_scale(0, coordinate_axis);
      replot = TRUE;
      break;

    case POSTSCRIPT:
/* save postscript */
      if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
				 "intro.ps", &file_name, 'w', 0,
				 save_command) ) != NULL){
	col_flag = get_yes_no( io_ptr, "Do you want to save color "
				 "information>", LEVEL );
	PL_postscript( 0, fp, col_flag );
	fclose( fp );
      }
      break;

    case SHOW:
/* show */
      printf("Curve Plotmode:\n");
      printf("Curve: %s\n", (xcog_data_ptr->curve_plot_mode & 1)? "on" : "off");
      printf("Curve direction: %s\n", (xcog_data_ptr->curve_plot_mode & 2)? "on" : 
	     "off");
      printf("Curve tickmarks: %s\n", (xcog_data_ptr->curve_plot_mode & 4)? "on" : 
	     "off");
      printf("\n");
      printf("Grid Plotmode:\n");
      printf("Grid boundary: %s\n", (xcog_data_ptr->grid_plot_mode & 1)? "on" : "off");
      printf("Grid directional arrows: %s\n", 
	     (xcog_data_ptr->grid_plot_mode & 2)? "on" : "off");
      printf("Grid tickmarks: %s\n", (xcog_data_ptr->grid_plot_mode & 4)? "on" : "off");
      printf("All gridpoints: %s\n", (xcog_data_ptr->grid_plot_mode & 8)? "on" : "off");
      printf("Boundary condition: %s\n", 
	     (xcog_data_ptr->grid_plot_mode & 16)? "on" : "off");
      printf("Curve value: %s\n", (xcog_data_ptr->grid_plot_mode & 32)? "on" : "off");
      printf("\n");
      printf("Overlapping Grid Plotmode:\n");
      printf("Composite grid: %s\n", (xcog_data_ptr->comp_plot_mode & 1)? "on" : "off");
      printf("Interpolation points: %s\n", (xcog_data_ptr->comp_plot_mode & 2)? 
	     "on" : "off");
      printf("Bad points: %s\n", (xcog_data_ptr->comp_plot_mode & 4)? "on" : "off");
      printf("Flag array: %s\n", (xcog_data_ptr->comp_plot_mode & 8)? "on" : "off");
      if (boundary_ptr)
	sprintf(curve_string, " (curve_label = %i)", boundary_ptr->curve_label);
      else
	curve_string[0] = '\0';
      printf("Physical boundary: %s%s", (boundary_ptr)? "on" : "off", curve_string);
      printf("\n");
      break;

    case CLEAR:
      xcog_data_ptr->curve_plot_mode = 0;
      xcog_data_ptr->grid_plot_mode = 0;
      xcog_data_ptr->comp_plot_mode = 0;
      replot = TRUE;
      break;

    case HELP:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case EXIT:
/* end */
      quit = 1;
      break;

    default:
      ;
    }

  }
  while(!quit);

}
