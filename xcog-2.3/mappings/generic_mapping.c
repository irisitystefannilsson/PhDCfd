#include "mappings.h"

/* static prototypes */
static generic_mapping *
new_generic_mapping( char *name);
/* static prototypes */

void 
mark_mapping_bndry(input_output *io_ptr, generic_mapping *map_ptr){
  char prompt[80];
  int quit, replot, icom;
  const int no_indent=0, plot_mode=1|2|4;

#include "mark_mapping_bndry_com.h"

  sprintf(prompt, "%s: mark boundary>", map_ptr->grid_name);
  quit = FALSE;
  replot = TRUE;
  PL_window(1, map_ptr->xyab);
  do{

    if (replot){
/* plot the mapping */
      PL_erase(1);
      PL_start_plot(1);
      plot_generic_mapping( map_ptr, plot_mode );
      PL_end_plot();
      replot = FALSE;
    }

    replot = TRUE;
    switch (get_command(io_ptr, prompt, COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT)){

    case LOW_R:
      map_ptr->gap_low_r[0] = get_real(io_ptr, "Starting s-value for external gap: ", 
			      map_ptr->gap_low_r[0], no_indent);
      map_ptr->gap_low_r[1] = get_real(io_ptr, "Ending s-value for external gap: ", 
			      map_ptr->gap_low_r[1], no_indent);
      replot = TRUE;
      break;

    case HIGH_R:
      map_ptr->gap_high_r[0] = get_real(io_ptr, "Starting s-value for external gap: ", 
			       map_ptr->gap_high_r[0], no_indent);
      map_ptr->gap_high_r[1] = get_real(io_ptr, "Ending s-value for external gap: ", 
			       map_ptr->gap_high_r[1], no_indent);
      replot = TRUE;
      break;

    case LOW_S:
      map_ptr->gap_low_s[0] = get_real(io_ptr, "Starting r-value for external gap: ", 
			      map_ptr->gap_low_s[0], no_indent);
      map_ptr->gap_low_s[1] = get_real(io_ptr, "Ending r-value for external gap: ", 
			      map_ptr->gap_low_s[1], no_indent);
      replot = TRUE;
      break;

    case HIGH_S:
      map_ptr->gap_high_s[0] = get_real(io_ptr, "Starting r-value for external gap: ", 
					map_ptr->gap_high_s[0], no_indent);
      map_ptr->gap_high_s[1] = get_real(io_ptr, "Ending r-value for external gap: ", 
					map_ptr->gap_high_s[1], no_indent);
      replot = TRUE;
      break;

    case SHOW:
      if (map_ptr->gap_low_s[0] > 1.0 || map_ptr->gap_low_s[1] < 0.0 ||
	  map_ptr->gap_low_s[1] < map_ptr->gap_low_s[0])
	printf("There is no external gap along the s=0 boundary.\n");
      else
	printf("The external gap along the s=0 boundary is at %f <= r <= %f\n",
	       map_ptr->gap_low_s[0], map_ptr->gap_low_s[1]);

      if (map_ptr->gap_high_s[0] > 1.0 || map_ptr->gap_high_s[1] < 0.0 ||
	  map_ptr->gap_high_s[1] < map_ptr->gap_high_s[0])
	printf("There is no external gap along the s=1 boundary.\n");
      else
	printf("The external gap along the s=1 boundary is at %f <= r <= %f\n",
	       map_ptr->gap_high_s[0], map_ptr->gap_high_s[1]);

      if (map_ptr->gap_low_r[0] > 1.0 || map_ptr->gap_low_r[1] < 0.0 ||
	  map_ptr->gap_low_r[1] < map_ptr->gap_low_r[0])
	printf("There is no external gap along the r=0 boundary.\n");
      else
	printf("The external gap along the r=0 boundary is at %f <= s <= %f\n",
	       map_ptr->gap_low_r[0], map_ptr->gap_low_r[1]);

      if (map_ptr->gap_high_r[0] > 1.0 || map_ptr->gap_high_r[1] < 0.0 ||
	  map_ptr->gap_high_r[1] < map_ptr->gap_high_r[0])
	printf("There is no external gap along the r=1 boundary.\n");
      else
	printf("The external gap along the r=1 boundary is at %f <= s <= %f\n",
	       map_ptr->gap_high_r[0], map_ptr->gap_high_r[1]);

      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
      quit = TRUE;
      break;

    default:
      ;
    }
  }
  while(!quit);
}

generic_mapping *
read_mapping( int32 map_dir, generic_curve *first_c ){
  generic_mapping *grid_ptr;

  int32 specific_dir;

/* create a generic mapping structure */
  grid_ptr = (generic_mapping *) malloc( sizeof( generic_mapping ) );

/* important to initialize the pointers */
  grid_ptr->next_m = grid_ptr->prev_m = NULL;

/* read the name */
  grid_ptr->grid_name = hget_string( "grid_name", map_dir );

  hget_int( &(grid_ptr->priority), "priority", map_dir );
  hget_int( &(grid_ptr->inverse_known), "inverse_known", map_dir);
  hget_int( &(grid_ptr->grid_type), "grid_type", map_dir );
  hget_int( &(grid_ptr->r_points), "r_points", map_dir );
  hget_int( &(grid_ptr->s_points), "s_points", map_dir );

  hget_int( &(grid_ptr->r_period), "r_period", map_dir );
  hget_int( &(grid_ptr->s_period), "s_period", map_dir );

  grid_ptr->bc_ptr = hget_int_array_2d( "bc", map_dir );
  grid_ptr->curve_ptr = hget_int_array_2d( "curve", map_dir );

  hget_real( &(grid_ptr->x_min_0), "x_min_0", map_dir );
  hget_real( &(grid_ptr->x_max_0), "x_max_0", map_dir );
  hget_real( &(grid_ptr->y_min_0), "y_min_0", map_dir );
  hget_real( &(grid_ptr->y_max_0), "y_max_0", map_dir );

  hget_real( &(grid_ptr->x_min), "x_min", map_dir );
  hget_real( &(grid_ptr->x_max), "x_max", map_dir );
  hget_real( &(grid_ptr->y_min), "y_min", map_dir );
  hget_real( &(grid_ptr->y_max), "y_max", map_dir );

  hget_int( &(grid_ptr->grid_plot_mode), "grid_plot_mode", map_dir );
  hget_real( &(grid_ptr->xyab[0][0]), "xyab00", map_dir);
  hget_real( &(grid_ptr->xyab[1][0]), "xyab10", map_dir);
  hget_real( &(grid_ptr->xyab[0][1]), "xyab01", map_dir);
  hget_real( &(grid_ptr->xyab[1][1]), "xyab11", map_dir);

/* color info */
  hget_double( &(grid_ptr->grid_color.r), "grid_color.r", map_dir );
  hget_double( &(grid_ptr->grid_color.g), "grid_color.g", map_dir );
  hget_double( &(grid_ptr->grid_color.b), "grid_color.b", map_dir );

/* rotation, translation, scaling */
  hget_real( &(grid_ptr->theta_rot), "theta_rot", map_dir );
  hget_real( &(grid_ptr->cos_theta), "cos_theta", map_dir );
  hget_real( &(grid_ptr->sin_theta), "sin_theta", map_dir );
  hget_real( &(grid_ptr->x_trans), "x_trans", map_dir );
  hget_real( &(grid_ptr->y_trans), "y_trans", map_dir );
  hget_real( &(grid_ptr->scaling), "scaling", map_dir );

/* external gaps */
  if (!hget_real( &(grid_ptr->gap_low_r[0]), "gap_low_r_0", map_dir ))
    grid_ptr->gap_low_r[0] =  2.0;
  if (!hget_real( &(grid_ptr->gap_low_r[1]), "gap_low_r_1", map_dir ))
    grid_ptr->gap_low_r[1] = -1.0;
  if (!hget_real( &(grid_ptr->gap_high_r[0]), "gap_high_r_0", map_dir ))
    grid_ptr->gap_high_r[0] =  2.0;
  if (!hget_real( &(grid_ptr->gap_high_r[1]), "gap_high_r_1", map_dir ))
    grid_ptr->gap_high_r[1] = -1.0;
  if (!hget_real( &(grid_ptr->gap_low_s[0]), "gap_low_s_0", map_dir ))
    grid_ptr->gap_low_s[0] =  2.0;
  if (!hget_real( &(grid_ptr->gap_low_s[1]), "gap_low_s_1", map_dir ))
    grid_ptr->gap_low_s[1] = -1.0;
  if (!hget_real( &(grid_ptr->gap_high_s[0]), "gap_high_s_0", map_dir ))
    grid_ptr->gap_high_s[0] =  2.0;
  if (!hget_real( &(grid_ptr->gap_high_s[1]), "gap_high_s_1", map_dir ))
    grid_ptr->gap_high_s[1] = -1.0;

/* get the directory for the curve-specific parameters */
  if ((specific_dir = locate_dir("specific_params", map_dir)) != -1){

/* assign the appropriate mapping-specific pointers */
    switch( grid_ptr->grid_type ){

    case 1:
/* cartesian */
      if (!read_cartesian_mapping( specific_dir, grid_ptr, first_c )){
/* free the space */
	grid_ptr->cleanup_mapping_data = NULL;
	delete_generic_mapping( grid_ptr );
	grid_ptr = NULL;
      }
      break;

    case 4:
/* normal-curve */
      if (!read_normal_curve( specific_dir, grid_ptr, first_c )){
/* free the space */
	grid_ptr->cleanup_mapping_data = NULL;
	delete_generic_mapping( grid_ptr );
	grid_ptr = NULL;
      }
      break;

    case 5:
/* linear interpolation */
      if (!read_linear_interp( specific_dir, grid_ptr, first_c )){
/* free the space */
	grid_ptr->cleanup_mapping_data = NULL;
	delete_generic_mapping( grid_ptr );
	grid_ptr = NULL;
      }
      break;

    case 6:
/* discrete point mapping */
      if (!read_disc_point_map( specific_dir, grid_ptr, first_c )){
/* free the space */
	grid_ptr->cleanup_mapping_data = NULL;
	delete_generic_mapping( grid_ptr );
	grid_ptr = NULL;
      }
      break;

#ifdef THE_GAR
    case 7:
/* Theodorsen-Garrick mapping */
      if (!read_tg_mapping( specific_dir, grid_ptr, first_c )){
/* free the space */
	grid_ptr->cleanup_mapping_data = NULL;
	delete_generic_mapping( grid_ptr );
	grid_ptr = NULL;
      }
      break;
#endif

    case 8:
/* discrete point mapping */
      if (!read_hyp_map( specific_dir, grid_ptr, first_c )){
/* free the space */
	grid_ptr->cleanup_mapping_data = NULL;
	delete_generic_mapping( grid_ptr );
	grid_ptr = NULL;
      }
      break;

    default:
      printf("Unknown mapping type: %i in read_mapping\n", grid_ptr->grid_type);
/* free the space */
      grid_ptr->cleanup_mapping_data = NULL;
      delete_generic_mapping( grid_ptr );
      grid_ptr = NULL;

    } /* end switch */

/* terminate access to the specific dir */
    Vdetach( specific_dir );
  } /* end if locate dir */

  return grid_ptr;
}

void 
write_mapping( int32 map_dir, generic_mapping *grid_ptr ){

  int32 specific_dir;

/* save the name */
  hput_string( grid_ptr->grid_name, "grid_name", map_dir );

  hput_int( grid_ptr->priority, "priority", map_dir );
  hput_int( grid_ptr->inverse_known, "inverse_known", map_dir );
  hput_int( grid_ptr->grid_type, "grid_type", map_dir );
  hput_int( grid_ptr->r_points, "r_points", map_dir );
  hput_int( grid_ptr->s_points, "s_points", map_dir );

  hput_int( grid_ptr->r_period, "r_period", map_dir );
  hput_int( grid_ptr->s_period, "s_period", map_dir );
  
  hput_int_array_2d( grid_ptr->bc_ptr,    "bc",    map_dir );
  hput_int_array_2d( grid_ptr->curve_ptr, "curve", map_dir );

  hput_real( grid_ptr->x_min_0, "x_min_0", map_dir );
  hput_real( grid_ptr->x_max_0, "x_max_0", map_dir );
  hput_real( grid_ptr->y_min_0, "y_min_0", map_dir );
  hput_real( grid_ptr->y_max_0, "y_max_0", map_dir );

  hput_real( grid_ptr->x_min, "x_min", map_dir );
  hput_real( grid_ptr->x_max, "x_max", map_dir );
  hput_real( grid_ptr->y_min, "y_min", map_dir );
  hput_real( grid_ptr->y_max, "y_max", map_dir );

  hput_int( grid_ptr->grid_plot_mode, "grid_plot_mode", map_dir );
  hput_real( grid_ptr->xyab[0][0], "xyab00", map_dir );
  hput_real( grid_ptr->xyab[0][1], "xyab01", map_dir );
  hput_real( grid_ptr->xyab[1][0], "xyab10", map_dir );
  hput_real( grid_ptr->xyab[1][1], "xyab11", map_dir );

/* color info */
  hput_double(grid_ptr->grid_color.r, "grid_color.r", map_dir );
  hput_double(grid_ptr->grid_color.g, "grid_color.g", map_dir );
  hput_double(grid_ptr->grid_color.b, "grid_color.b", map_dir );

/* rotation, translation, scaling */
  hput_real( grid_ptr->theta_rot, "theta_rot", map_dir );
  hput_real( grid_ptr->cos_theta, "cos_theta", map_dir );
  hput_real( grid_ptr->sin_theta, "sin_theta", map_dir );
  hput_real( grid_ptr->x_trans,   "x_trans",   map_dir );
  hput_real( grid_ptr->y_trans,   "y_trans",   map_dir );
  hput_real( grid_ptr->scaling,   "scaling",   map_dir );

/* external gaps */
  hput_real( grid_ptr->gap_low_r[0],  "gap_low_r_0",  map_dir );
  hput_real( grid_ptr->gap_low_r[1],  "gap_low_r_1",  map_dir );
  hput_real( grid_ptr->gap_high_r[0], "gap_high_r_0", map_dir );
  hput_real( grid_ptr->gap_high_r[1], "gap_high_r_1", map_dir );
  hput_real( grid_ptr->gap_low_s[0],  "gap_low_s_0",  map_dir );
  hput_real( grid_ptr->gap_low_s[1],  "gap_low_s_1",  map_dir );
  hput_real( grid_ptr->gap_high_s[0], "gap_high_s_0", map_dir );
  hput_real( grid_ptr->gap_high_s[1], "gap_high_s_1", map_dir );

/* make a directory for the mapping-specific parameters */
  if ((specific_dir = create_dir("specific_params", "specific_params", 
				 map_dir)) != -1){
/* save the specific part of the data */
    (*grid_ptr->write_specific_mapping)( specific_dir, grid_ptr->mapping_data_ptr );
    Vdetach(specific_dir);
  }
  else
    printf("Error: save_mapping: unable to create a directory for the type-specific "
	   "stuff\n");
}


generic_mapping *
choose_mapping(input_output *io_ptr, char *name, generic_curve *first_c ){
  char prompt[80];
  int icom, quit;
  generic_curve *new_curve_ptr, *new_s0_curve_ptr, *new_s1_curve_ptr;
  generic_mapping *grid_ptr=NULL;

#include "choose_mapping_com.h"

  sprintf(prompt,"set mapping type %s>", name);

  do{

    icom = get_command(io_ptr, prompt, COMMAND, LAST_COM, LEVEL,
		       SAVE_ON_COPY, ARGUMENT);

    quit = 0;
    switch (icom) {
/* mappings */

    case CARTESIAN:
/* make mapping */
      grid_ptr = new_generic_mapping( name );
      init_cartesian_mapping( grid_ptr, 0.0, 1.0, 0.0, 1.0 );
      quit = 1;
      break;

    case NORMAL_CURVE:
      if (first_c == NULL){
	printf("This mapping needs at least one curve, but the curve list is "
	       "empty.\n");
      }
      else if ((new_curve_ptr=get_curve_ptr( io_ptr, first_c )) != NULL){
/* initialize the grid */  
	grid_ptr = new_generic_mapping( name );
	init_normal_curve( io_ptr, grid_ptr, new_curve_ptr, "new" );
	quit = 1;
      }

      break;

    case HYPERBOLIC:
      if (first_c == NULL){
	printf("This mapping needs at least one curve, but the curve list is "
	       "empty.\n");
      }
      else if ((new_curve_ptr=get_curve_ptr( io_ptr, first_c )) != NULL){
/* initialize the grid */  
	grid_ptr = new_generic_mapping( name );
	init_hyp_grid( io_ptr, grid_ptr, new_curve_ptr );
	quit = 1;
      }

      break;

    case LINEAR_INTERP:
/* select boundary curves */ 
      if ( first_c == NULL || first_c->next == NULL ){
	printf("At least 2 curves are needed for this mapping.\n");
      }
      else if ( get_2curve_ptr( io_ptr, first_c, &new_s0_curve_ptr, 
			  &new_s1_curve_ptr, LEVEL+1) ){
/* initialize the grid */
	grid_ptr = new_generic_mapping( name );
	init_linear_interp(io_ptr, grid_ptr, new_s0_curve_ptr, new_s1_curve_ptr);
	quit = 1;
      }

      break;

    case DISC_POINT:
/* interpolate-grid-points */ 
      grid_ptr = new_generic_mapping( name );
      if ( init_disc_point_map( io_ptr, grid_ptr ) != NULL){
	quit = 1;
      }
      else{
/* remove the generic mapping if the initialization was unsuccessful */
	delete_generic_mapping( grid_ptr );
      }

      break;

#ifdef THE_GAR
    case TEO_GAR: /* theodorsen-garrick */
      printf("theodorsen-garrick\n");
      grid_ptr = new_generic_mapping( name );
      if (init_tg_mapping( io_ptr, grid_ptr )){
	quit = 1;
      }
      else{
/* remove the generic mapping if the initialization was unsuccessful */
	delete_generic_mapping( grid_ptr );
      }
      break;
#endif

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );

      break;

    case CANCEL:
/* cancel, return 0 to indicate unsucessful choice of mapping */
      return NULL;
      break;

    default:
      ;
    }
  }
  while(!quit);

/* call the specific setup routine to enter all the mapping specific data */
  set_mapping( io_ptr, grid_ptr, first_c );

  return grid_ptr;
}



void 
delete_generic_mapping( generic_mapping *grid_ptr ){

/* deallocate the name */
  free( grid_ptr->grid_name );

/* free the space of the mapping */
  grid_ptr->bc_ptr = delete_int_array_2d( grid_ptr->bc_ptr );
  grid_ptr->curve_ptr = delete_int_array_2d( grid_ptr->curve_ptr );

/* free the data associated with the mapping */
  if (grid_ptr->cleanup_mapping_data != NULL)
    grid_ptr->mapping_data_ptr = 
      (*grid_ptr->cleanup_mapping_data)(grid_ptr->mapping_data_ptr);

/* free the structure itself */
  free(grid_ptr);
}


generic_mapping *
copy_mapping( char *name, generic_mapping *old_grid_ptr){
  generic_mapping *grid_ptr;
  int i, j;

#define old_bc(i,j)    compute_index_2d(old_grid_ptr->bc_ptr,i,j)
#define old_curve(i,j) compute_index_2d(old_grid_ptr->curve_ptr,i,j)

  grid_ptr = (generic_mapping *)malloc(sizeof(generic_mapping));
  
/* initialize all fields */
  grid_ptr->prev_m = NULL;
  grid_ptr->next_m = NULL;
  grid_ptr->priority = 0;

/* copy the name */
  grid_ptr->grid_name = (char *) malloc( (strlen(name)+1)*sizeof(char) );
  grid_ptr->grid_name = strcpy( grid_ptr->grid_name, name );

/* save the pointers to the data and the functions for
   evaluating the mapping */
  grid_ptr->inverse_known = old_grid_ptr->inverse_known;
  grid_ptr->grid_mapping = old_grid_ptr->grid_mapping;
  grid_ptr->inverse_grid_mapping = old_grid_ptr->inverse_grid_mapping;
  grid_ptr->mapping_data_ptr = 
    (*old_grid_ptr->copy_mapping_data)( old_grid_ptr->mapping_data_ptr );
  grid_ptr->cleanup_mapping_data = old_grid_ptr->cleanup_mapping_data;
  grid_ptr->set_mapping = old_grid_ptr->set_mapping;
  grid_ptr->write_specific_mapping = old_grid_ptr->write_specific_mapping;
  grid_ptr->copy_mapping_data = old_grid_ptr->copy_mapping_data;

/* set the grid type */
  grid_ptr->grid_type = old_grid_ptr->grid_type; 

/* set the number of grid lines */
  grid_ptr->r_points = old_grid_ptr->r_points;
  grid_ptr->s_points = old_grid_ptr->s_points;

  grid_ptr->r_period = old_grid_ptr->r_period;
  grid_ptr->s_period = old_grid_ptr->s_period;

  grid_ptr->bc_ptr    = create_int_array_2d(2,2);
  grid_ptr->curve_ptr = create_int_array_2d(2,2);

  grid_ptr->x_min_0 = old_grid_ptr->x_min_0;
  grid_ptr->x_max_0 = old_grid_ptr->x_max_0;
  grid_ptr->y_min_0 = old_grid_ptr->y_min_0;
  grid_ptr->y_max_0 = old_grid_ptr->y_max_0;

  grid_ptr->x_min = old_grid_ptr->x_min;
  grid_ptr->x_max = old_grid_ptr->x_max;
  grid_ptr->y_min = old_grid_ptr->y_min;
  grid_ptr->y_max = old_grid_ptr->y_max;

  grid_ptr->grid_plot_mode = old_grid_ptr->grid_plot_mode;

  grid_ptr->xyab[0][0] = old_grid_ptr->xyab[0][0];
  grid_ptr->xyab[0][1] = old_grid_ptr->xyab[0][1];
  grid_ptr->xyab[1][0] = old_grid_ptr->xyab[1][0];
  grid_ptr->xyab[1][1] = old_grid_ptr->xyab[1][1];

/* rotation, translation, scaling */
  grid_ptr->theta_rot = old_grid_ptr->theta_rot;
  grid_ptr->cos_theta = old_grid_ptr->cos_theta;
  grid_ptr->sin_theta = old_grid_ptr->sin_theta;
  grid_ptr->x_trans = old_grid_ptr->x_trans;
  grid_ptr->y_trans = old_grid_ptr->y_trans;
  grid_ptr->scaling = old_grid_ptr->scaling;

/* external gaps */
  for (i=0; i<2; i++){
    grid_ptr->gap_low_r[i]  = old_grid_ptr->gap_low_r[i];
    grid_ptr->gap_high_r[i] = old_grid_ptr->gap_high_r[i];
    grid_ptr->gap_low_s[i]  = old_grid_ptr->gap_low_s[i];
    grid_ptr->gap_high_s[i] = old_grid_ptr->gap_high_s[i];
  }

/* initialize the 2-d arrays */
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++){
      bc(i,j) = old_bc(i,j);
      curve(i,j) = old_curve(i,j);
    }

/* get a color */
  get_new_color( &(grid_ptr->grid_color), 'g' ); 

  return grid_ptr;

#undef old_bc
#undef old_curve
}

static generic_mapping *
new_generic_mapping( char *name){
  generic_mapping *grid_ptr;
  int i, j;

  if ((grid_ptr = (generic_mapping *)malloc(sizeof(generic_mapping)) ) == NULL)
    printf("memory error in new_generic_mapping\n");
  
/* initialize all fields */
  grid_ptr->prev_m = NULL;
  grid_ptr->next_m = NULL;
  grid_ptr->priority = 0;
  if ((grid_ptr->grid_name = (char *) malloc( (strlen(name)+1)*sizeof(char) )) == NULL)
    printf("memory error in new_generic_mapping\n");
  grid_ptr->grid_name = strcpy( grid_ptr->grid_name, name );

/* save the pointers to the data and the functions for
   evaluating the mapping */
  grid_ptr->inverse_known = 0;
  grid_ptr->grid_mapping = NULL;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->mapping_data_ptr = NULL;
  grid_ptr->cleanup_mapping_data = NULL;
  grid_ptr->set_mapping = NULL;
  grid_ptr->write_specific_mapping = NULL;
  grid_ptr->copy_mapping_data = NULL;

/* set the grid type and the periodicity flag */
  grid_ptr->grid_type= 0; /* unknown grid type */

/* set the number of grid lines */
  grid_ptr->r_points = 10;
  grid_ptr->s_points = 10;

  grid_ptr->r_period = 0;
  grid_ptr->s_period = 0;

  grid_ptr->bc_ptr        = create_int_array_2d(2,2);
  grid_ptr->curve_ptr     = create_int_array_2d(2,2);

  grid_ptr->x_min_0 = 0.0;
  grid_ptr->x_max_0 = 1.0;
  grid_ptr->y_min_0 = 0.0;
  grid_ptr->y_max_0 = 1.0;

  grid_ptr->x_min = 0.0;
  grid_ptr->x_max = 1.0;
  grid_ptr->y_min = 0.0;
  grid_ptr->y_max = 1.0;

  grid_ptr->grid_plot_mode = 1 + 2;

  grid_ptr->xyab[0][0] = 0.0;
  grid_ptr->xyab[0][1] = 1.0;
  grid_ptr->xyab[1][0] = 0.0;
  grid_ptr->xyab[1][1] = 1.0;

/* rotation, translation, scaling */
  grid_ptr->theta_rot = 0.0;
  grid_ptr->cos_theta = 1.0;
  grid_ptr->sin_theta = 0.0;
  grid_ptr->x_trans = 0.0;
  grid_ptr->y_trans = 0.0;
  grid_ptr->scaling = 1.0;

/* external gaps */
  grid_ptr->gap_low_r[0] =  2.0;
  grid_ptr->gap_low_r[1] = -1.0;
  grid_ptr->gap_high_r[0] =  2.0;
  grid_ptr->gap_high_r[1] = -1.0;
  grid_ptr->gap_low_s[0] =  2.0;
  grid_ptr->gap_low_s[1] = -1.0;
  grid_ptr->gap_high_s[0] =  2.0;
  grid_ptr->gap_high_s[1] = -1.0;

/* initialize the 2-d arrays */
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++){
      bc(i,j) = 0;
      curve(i,j) = 0;
    }

/* get a color */
  get_new_color( &(grid_ptr->grid_color), 'g' ); 

  return grid_ptr;
}

void 
plot_generic_mapping(generic_mapping *grid_ptr, int plot_mode ){
  int i, h_center, v_center, font_number, inside;
  char c_flag[20];
  grid_point gp, gp2;
  real pi, theta, length, sigma, sx, sy, px, py, mx, my, xp, yp, xm, ym
    , r_value, s_value, lr, ls, l_min, r_step, s_step;
  grid_quality quality;
/*
   plot_mode controls what is plotted according to:
   plot_mode = 1   grid boundary (only if plot_mode != 8).
   plot_mode = 2   directional arrows.
   plot_mode = 4   grid tickmarks (only if plot_mode != 8)
   plot_mode = 8   all non-fictitious gridpoints.
   plot_mode = 16  boundary condition.
   plot_mode = 32  curve value.
   plot_mode = 64  title
   plot_mode = 128 corner classification
   plot_mode = 256 grid quality
   Different options are combined by logical or'ing, i.e. plot_mode=1|2 gives both
   option 1 and 2. 
*/

/* get the grid sizes */
  r_step = 1.0/(grid_ptr->r_points - 1);
  s_step = 1.0/(grid_ptr->s_points - 1);

/* set the color */
  PL_color( grid_ptr->grid_color );

/* draw the boundary in solid */
  if ((plot_mode & 1) && !(plot_mode & 8) && grid_ptr->grid_mapping != NULL){
/* r=0 */
    if (curve(1,1) > 0){
      PL_line_width(2);
      PL_stop_dash();
    }
    else if (curve(1,1) < 0){
      PL_line_width(2);
      PL_start_dash();
    }
    else{
      PL_line_width(1);
      PL_stop_dash();
    }
    gp.r = 0.0;
    gp.s = 0.0;
    forward_mapping( &gp, grid_ptr );
    PL_poll();
    PL_move(gp.x, gp.y );
    inside = TRUE;
    for (gp.s = s_step; gp.s < 1.0 + 0.5*s_step;
	 gp.s += s_step){
      forward_mapping( &gp, grid_ptr );
/* check the gap */
      if (grid_ptr->gap_low_r[0] <= gp.s && gp.s <= grid_ptr->gap_low_r[1]){
	if (inside){
	  gp2.r = 0.0;
	  gp2.s = grid_ptr->gap_low_r[0];
	  forward_mapping( &gp2, grid_ptr );
	  PL_draw(gp2.x , gp2.y);
	}
	PL_move(gp.x, gp.y);
	inside = FALSE;
      }
      else{
	if (!inside){
	  gp2.r = 0.0;
	  gp2.s = grid_ptr->gap_low_r[1];
	  forward_mapping( &gp2, grid_ptr );
	  PL_move(gp2.x , gp2.y);
	}
	PL_draw(gp.x , gp.y);
	inside = TRUE;
      }
    }
/* draw markers where the external gap begins and ends */
    if (grid_ptr->gap_low_r[0] <= 1.0 && grid_ptr->gap_low_r[1] >= 0.0 &&
	 grid_ptr->gap_low_r[0] <= grid_ptr->gap_low_r[1]){
      gp.s = grid_ptr->gap_low_r[0];
      forward_mapping( &gp, grid_ptr );
      PL_marker(gp.x, gp.y, CROSS);
      gp.s = grid_ptr->gap_low_r[1];
      forward_mapping( &gp, grid_ptr );
      PL_marker(gp.x, gp.y, CROSS);
    }

/* r=1 */
    if (curve(2,1) > 0){
      PL_line_width(2);
      PL_stop_dash();
    }
    else if (curve(2,1) < 0){
      PL_line_width(2);
      PL_start_dash();
    }
    else{
      PL_line_width(1);
      PL_stop_dash();
    }
    gp.r = 1.0;
    gp.s = 0.0;
    forward_mapping( &gp, grid_ptr );
    PL_poll();
    PL_move(gp.x, gp.y );
    inside = TRUE;
    for (gp.s = s_step; gp.s < 1.0 + 0.5*s_step;
	 gp.s += s_step){
      forward_mapping( &gp, grid_ptr );
/* check the gap */
      if (grid_ptr->gap_high_r[0] <= gp.s && gp.s <= grid_ptr->gap_high_r[1]){
	if (inside){
	  gp2.r = 1.0;
	  gp2.s = grid_ptr->gap_high_r[0];
	  forward_mapping( &gp2, grid_ptr );
	  PL_draw(gp2.x , gp2.y);
	}
	PL_move(gp.x , gp.y);
	inside = FALSE;
      }
      else{
	if (!inside){
	  gp2.r = 1.0;
	  gp2.s = grid_ptr->gap_high_r[1];
	  forward_mapping( &gp2, grid_ptr );
	  PL_move(gp2.x , gp2.y);
	}
	PL_draw(gp.x , gp.y);
	inside = TRUE;
      }
    }
/* draw markers where the external gap begins and ends */
    if (grid_ptr->gap_high_r[0] <= 1.0 && grid_ptr->gap_high_r[1] >= 0.0 &&
	 grid_ptr->gap_high_r[0] <= grid_ptr->gap_high_r[1]){
      gp.s = grid_ptr->gap_high_r[0];
      forward_mapping( &gp, grid_ptr );
      PL_marker(gp.x, gp.y, CROSS);
      gp.s = grid_ptr->gap_high_r[1];
      forward_mapping( &gp, grid_ptr );
      PL_marker(gp.x, gp.y, CROSS);
    }

/* s=0 */
    if (curve(1,2) > 0){
      PL_line_width(2);
      PL_stop_dash();
    }
    else if (curve(1,2) < 0){
      PL_line_width(2);
      PL_start_dash();
    }
    else{
      PL_line_width(1);
      PL_stop_dash();
    }
    gp.s = 0.0;
    gp.r = 0.0;
    forward_mapping( &gp, grid_ptr );
    PL_poll();
    PL_move(gp.x, gp.y );
    inside = TRUE;
    for (gp.r = r_step; gp.r < 1.0 + 0.5*r_step;
	 gp.r += r_step){
      forward_mapping( &gp, grid_ptr );
/* check the gap */
      if (grid_ptr->gap_low_s[0] <= gp.r && gp.r <= grid_ptr->gap_low_s[1]){
	if (inside){
	  gp2.s = 0.0;
	  gp2.r = grid_ptr->gap_low_s[0];
	  forward_mapping( &gp2, grid_ptr );
	  PL_draw(gp2.x , gp2.y);
	}
	PL_move(gp.x, gp.y);
	inside = FALSE;
      }
      else{
	if (!inside){
	  gp2.s = 0.0;
	  gp2.r = grid_ptr->gap_low_s[1];
	  forward_mapping( &gp2, grid_ptr );
	  PL_move(gp2.x , gp2.y);
	}
	PL_draw(gp.x , gp.y);
	inside = TRUE;
      }
    } /* end for gp.r ... */
/* draw markers where the external gap begins and ends */
    if (grid_ptr->gap_low_s[0] <= 1.0 && grid_ptr->gap_low_s[1] >= 0.0 &&
	 grid_ptr->gap_low_s[0] <= grid_ptr->gap_low_s[1]){
      gp.r = grid_ptr->gap_low_s[0];
      forward_mapping( &gp, grid_ptr );
      PL_marker(gp.x, gp.y, CROSS);
      gp.r = grid_ptr->gap_low_s[1];
      forward_mapping( &gp, grid_ptr );
      PL_marker(gp.x, gp.y, CROSS);
    }

/* s=1 */
    if (curve(2,2) > 0){
      PL_line_width(2);
      PL_stop_dash();
    }
    else if (curve(2,2) < 0){
      PL_line_width(2);
      PL_start_dash();
    }
    else{
      PL_line_width(1);
      PL_stop_dash();
    }
    gp.s = 1.0;
    gp.r = 0.0;
    forward_mapping( &gp, grid_ptr );
    PL_poll();
    PL_move(gp.x, gp.y );
    inside = TRUE;
    for (gp.r = r_step; gp.r < 1.0 + 0.5*r_step;
	 gp.r += r_step){
      forward_mapping( &gp, grid_ptr );
/* check the gap */
      if (grid_ptr->gap_high_s[0] <= gp.r && gp.r <= grid_ptr->gap_high_s[1]){
	if (inside){
	  gp2.s = 1.0;
	  gp2.r = grid_ptr->gap_high_s[0];
	  forward_mapping( &gp2, grid_ptr );
	  PL_draw(gp2.x , gp2.y);
	}
	PL_move(gp.x, gp.y);
	inside = FALSE;
      }
      else{
	if (!inside){
	  gp2.s = 1.0;
	  gp2.r = grid_ptr->gap_high_s[1];
	  forward_mapping( &gp2, grid_ptr );
	  PL_move(gp2.x , gp2.y);
	}
	PL_draw(gp.x , gp.y);
	inside = TRUE;
      }
    }
/* draw markers where the external gap begins and ends */
    if (grid_ptr->gap_high_s[0] <= 1.0 && grid_ptr->gap_high_s[1] >= 0.0 &&
	 grid_ptr->gap_high_s[0] <= grid_ptr->gap_high_s[1]){
      gp.r = grid_ptr->gap_high_s[0];
      forward_mapping( &gp, grid_ptr );
      PL_marker(gp.x, gp.y, CROSS);
      gp.r = grid_ptr->gap_high_s[1];
      forward_mapping( &gp, grid_ptr );
      PL_marker(gp.x, gp.y, CROSS);
    }
/* reset the line width and stop any dashing */
    PL_line_width(1);
    PL_stop_dash();
  } /* end if plot_mode & 1 */

/* draw arrows to indicate the directions of the parameters */
/* estimate the relative lengths of the parameter lines r=0 and s=0 */
/* draw the non-fictitious points in solid */
  if ((plot_mode & 2) && grid_ptr->grid_mapping != NULL){
    gp.r = 0.0;
    gp.s = 0.0;
    forward_mapping( &gp, grid_ptr );
    lr = sqrt( gp.xr*gp.xr + gp.yr*gp.yr );
    ls = sqrt( gp.xs*gp.xs + gp.ys*gp.ys );
    if ( ls > lr ){
      r_value = 0.4;
      s_value = r_value * lr / ls;
      l_min = lr;
    }
    else{
      s_value = 0.4;
      r_value = s_value * ls / lr;
      l_min = ls;
    }
/* length and angle */
    pi = 4.0*atan(1.0);
    theta = 20.0*pi/180.0;
    length = 0.05*real_max(grid_ptr->x_max-grid_ptr->x_min, 
			   grid_ptr->y_max-grid_ptr->y_min);
    length = real_min( 0.25*l_min, length );
/* get the point on the boundary s=0 */
    gp.s = 0.0;
    gp.r = r_value;
    forward_mapping( &gp, grid_ptr );
/* find the tangent */
    sigma = 1.0/sqrt(gp.xr*gp.xr + gp.yr*gp.yr);
    sx = sigma*gp.xr;
    sy = sigma*gp.yr;
/* direction of the first leg */
    px = cos(theta)*sx - sin(theta)*sy;
    py = sin(theta)*sx + cos(theta)*sy;
/* direction of the second leg */
    mx = cos(-theta)*sx - sin(-theta)*sy;
    my = sin(-theta)*sx + cos(-theta)*sy;
/* all the way out on the first leg */
    xp = gp.x - length*px;
    yp = gp.y - length*py;
/* all the way out on the second leg */
    xm = gp.x - length*mx;
    ym = gp.y - length*my;
/* draw it */
    PL_move(xp, yp);
    PL_draw(gp.x, gp.y);
    PL_draw(xm, ym);
/* draw an 'r' next to the arrow */
    xp = gp.x + 0.5*length * sy;
    yp = gp.y - 0.5*length * sx;
    PL_plot_string("r", xp, yp, 0, 0, 3);


/* get the point on the boundary r=0 */
    gp.s = s_value;
    gp.r = 0.0;
    forward_mapping( &gp, grid_ptr );
/* find the tangent */
    sigma = 1.0/sqrt(gp.xs*gp.xs + gp.ys*gp.ys);
    sx = sigma*gp.xs;
    sy = sigma*gp.ys;
/* direction of the first leg */
    px = cos(theta)*sx - sin(theta)*sy;
    py = sin(theta)*sx + cos(theta)*sy;
/* direction of the second leg */
    mx = cos(-theta)*sx - sin(-theta)*sy;
    my = sin(-theta)*sx + cos(-theta)*sy;
/* all the way out on the first leg */
    xp = gp.x - length*px;
    yp = gp.y - length*py;
/* all the way out on the second leg */
    xm = gp.x - length*mx;
    ym = gp.y - length*my;
/* draw it */
    PL_move(xp, yp);
    PL_draw(gp.x, gp.y);
    PL_draw(xm, ym);
/* draw an 's' next to the arrow */
    xp = gp.x - 0.5*length * sy;
    yp = gp.y + 0.5*length * sx;
    PL_plot_string("s", xp, yp, 0, 0, 3);

  }

/* draw tickmarks on the grid boundaries */
  if ((plot_mode & 4) && !(plot_mode & 8) && grid_ptr->grid_mapping != NULL){
    length = 0.02*real_max(grid_ptr->x_max-grid_ptr->x_min, 
			   grid_ptr->y_max-grid_ptr->y_min);
/* s=0 */
    gp.s = 0.0;
    PL_poll();
    for (i=1; i<= grid_ptr->r_points; i++){
      gp.r = (i-1) * r_step;
      forward_mapping( &gp, grid_ptr );
      sigma = 1.0/sqrt(gp.xr*gp.xr + gp.yr*gp.yr);
/* tangent */
      sx = sigma*gp.xr;
      sy = sigma*gp.yr;
/* start and end points */
      xp = gp.x - length*sy;
      yp = gp.y + length*sx;
      PL_move(gp.x, gp.y);
      PL_draw(xp, yp);
    }

/* s=1 */
    gp.s = 1.0;
    PL_poll();
    for (i=1; i<= grid_ptr->r_points; i++){
      gp.r = (i-1) * r_step;
      forward_mapping( &gp, grid_ptr );
      sigma = 1.0/sqrt(gp.xr*gp.xr + gp.yr*gp.yr);
/* tangent */
      sx = sigma*gp.xr;
      sy = sigma*gp.yr;
/* start and end points */
      xm = gp.x + length*sy;
      ym = gp.y - length*sx;
      PL_move(gp.x, gp.y);
      PL_draw(xm, ym);
    }

/* r=0 */
    gp.r = 0.0;
    PL_poll();
    for (i=1; i<= grid_ptr->s_points; i++){
      gp.s = (i-1) * s_step;
      forward_mapping( &gp, grid_ptr );
      sigma = 1.0/sqrt(gp.xs*gp.xs + gp.ys*gp.ys);
/* tangent */
      sx = sigma*gp.xs;
      sy = sigma*gp.ys;
/* start and end points */
      xm = gp.x + length*sy;
      ym = gp.y - length*sx;
      PL_move(gp.x, gp.y);
      PL_draw(xm, ym);
    }

/* r=1 */
    gp.r = 1.0;
    PL_poll();
    for (i=1; i<= grid_ptr->s_points; i++){
      gp.s = (i-1) * s_step;
      forward_mapping( &gp, grid_ptr );
      sigma = 1.0/sqrt(gp.xs*gp.xs + gp.ys*gp.ys);
/* tangent */
      sx = sigma*gp.xs;
      sy = sigma*gp.ys;
/* start and end points */
      xp = gp.x - length*sy;
      yp = gp.y + length*sx;
      PL_move(gp.x, gp.y);
      PL_draw(xp, yp);
    }
  }

/* draw the non-fictitious points in solid */
  if ((plot_mode & 8) && grid_ptr->grid_mapping != NULL){
/*    PL_shade(1.0); */
/* iso i-lines */
    for (gp.r = 0.0; gp.r < 1.0+0.5*r_step; 
	 gp.r += r_step){
      gp.s = 0.0;
      forward_mapping( &gp, grid_ptr );
      PL_poll();
      PL_move(gp.x, gp.y );
      for (gp.s = s_step; gp.s < 1.0 + 0.5*s_step;
	   gp.s += s_step){
	forward_mapping( &gp, grid_ptr );
	PL_draw(gp.x , gp.y );
      }
    }

/* iso j-lines */
    for (gp.s = 0.0; gp.s < 1.0 + 0.5*s_step;
	 gp.s += s_step){
      gp.r = 0.0;
      forward_mapping( &gp, grid_ptr );
      PL_poll();
      PL_move(gp.x, gp.y );
      for (gp.r = r_step; gp.r < 1.0+0.5*r_step; 
	   gp.r += r_step){
	forward_mapping( &gp, grid_ptr );
	PL_draw(gp.x , gp.y );
      }
    }
  }

/* draw the value of bc close to each boundary */
  if ((plot_mode & 16) && grid_ptr->grid_mapping != NULL ){
    h_center = 0; v_center = 1; 
    font_number=3;
/* r=0 */
    sprintf(c_flag, "bc=%1d", bc(1,1));
    gp.r = 0.0;
    gp.s = 0.5;
    forward_mapping( &gp, grid_ptr );
    PL_plot_string(c_flag, gp.x, gp.y, h_center, v_center, font_number);
/* r=1 */
    sprintf(c_flag, "bc=%1d", bc(2,1));
    gp.r = 1.0;
    gp.s = 0.5;
    forward_mapping( &gp, grid_ptr );
    PL_plot_string(c_flag, gp.x, gp.y, h_center, v_center, font_number);
/* s=0 */
    sprintf(c_flag, "bc=%1d", bc(1,2));
    gp.r = 0.5;
    gp.s = 0.0;
    forward_mapping( &gp, grid_ptr );
    PL_plot_string(c_flag, gp.x, gp.y, h_center, v_center, font_number);
/* s=1 */
    sprintf(c_flag, "bc=%1d", bc(2,2));
    gp.r = 0.5;
    gp.s = 1.0;
    forward_mapping( &gp, grid_ptr );
    PL_plot_string(c_flag, gp.x, gp.y, h_center, v_center, font_number);

  }

/* draw the value of curve close to each boundary */
  if ((plot_mode & 32) && grid_ptr->grid_mapping != NULL ){
    h_center = 0; v_center = -1; 
    font_number=3;
/* r=0 */
    sprintf(c_flag, "curve=%1d", curve(1,1));
    gp.r = 0.0;
    gp.s = 0.5;
    forward_mapping( &gp, grid_ptr );
    PL_plot_string(c_flag, gp.x, gp.y, h_center, v_center, font_number);
/* r=1 */
    sprintf(c_flag, "curve=%1d", curve(2,1));
    gp.r = 1.0;
    gp.s = 0.5;
    forward_mapping( &gp, grid_ptr );
    PL_plot_string(c_flag, gp.x, gp.y, h_center, v_center, font_number);
/* s=0 */
    sprintf(c_flag, "curve=%1d", curve(1,2));
    gp.r = 0.5;
    gp.s = 0.0;
    forward_mapping( &gp, grid_ptr );
    PL_plot_string(c_flag, gp.x, gp.y, h_center, v_center, font_number);
/* s=1 */
    sprintf(c_flag, "curve=%1d", curve(2,2));
    gp.r = 0.5;
    gp.s = 1.0;
    forward_mapping( &gp, grid_ptr );
    PL_plot_string(c_flag, gp.x, gp.y, h_center, v_center, font_number);

  }

/* title */
  if ( plot_mode & 64 ){
    PL_title( 1, grid_ptr->grid_name );
  }


/* title */
  if ( plot_mode & 256 && grid_ptr->grid_mapping != NULL ){
    get_grid_quality( &quality, grid_ptr );

/* min cell area */
    gp.r = (quality.i_area - 1)*r_step;
    gp.s = (quality.j_area - 1)*s_step;
    forward_mapping( &gp, grid_ptr );
    PL_marker( gp.x, gp.y, BOX );

/* min grid size in the r-direction */
    gp.r = (quality.i_msr - 1)*r_step;
    gp.s = (quality.j_msr - 1)*s_step;
    forward_mapping( &gp, grid_ptr );
    PL_marker( gp.x, gp.y, ASTERISK );

/* min grid size in the s-direction */
    gp.r = (quality.i_mss - 1)*r_step;
    gp.s = (quality.j_mss - 1)*s_step;
    forward_mapping( &gp, grid_ptr );
    PL_marker( gp.x, gp.y, CROSS );

/* max aspect ratio */
    gp.r = (quality.i_ar - 1)*r_step;
    gp.s = (quality.j_ar - 1)*s_step;
    forward_mapping( &gp, grid_ptr );
    PL_marker( gp.x, gp.y, FILLED_CIRCLE );

/* min orthogonality */
    gp.r = (quality.i_ortho - 1)*r_step;
    gp.s = (quality.j_ortho - 1)*s_step;
    forward_mapping( &gp, grid_ptr );
    PL_marker( gp.x, gp.y, FILLED_BOX );

/* max change in direction along constant r lines */
    gp.r = (quality.i_theta_r - 1)*r_step;
    gp.s = (quality.j_theta_r - 1)*s_step;
    forward_mapping( &gp, grid_ptr );
    PL_marker( gp.x, gp.y, TRIANGLE );

/* max relaive change in grid size along constant r lines */
    gp.r = (quality.i_change_r - 1)*r_step;
    gp.s = (quality.j_change_r - 1)*s_step;
    forward_mapping( &gp, grid_ptr );
    PL_marker( gp.x, gp.y, FILLED_TRIANGLE );

/* max change in direction along constant s lines */
    gp.r = (quality.i_theta_s - 1)*r_step;
    gp.s = (quality.j_theta_s - 1)*s_step;
    forward_mapping( &gp, grid_ptr );
    PL_marker( gp.x, gp.y, INV_TRIANGLE );

/* max relaive change in grid size along constant r lines */
    gp.r = (quality.i_change_s - 1)*r_step;
    gp.s = (quality.j_change_s - 1)*s_step;
    forward_mapping( &gp, grid_ptr );
    PL_marker( gp.x, gp.y, FILLED_INV_TRIANGLE );
  }

} /* end grid_plot */
#undef ncom


void 
set_boundary_condition(input_output *io_ptr, generic_mapping *grid_ptr){
  char prompt[80];
  int icom, quit, replot;
  const int level=1, no_indent=0, plot_mode=19;

#include "set_boundary_condition_com.h"

  sprintf(prompt, "%s: set boundary condition>", grid_ptr->grid_name);
  quit = 0;
  replot = 1;
  PL_window(1, grid_ptr->xyab);
  do{

    if (replot){
/* plot the mapping */
      PL_erase(1);
      PL_start_plot(1);
      plot_generic_mapping( grid_ptr, plot_mode );
      PL_end_plot();
    }

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {

    case 0: 
      bc(1,1) = get_int(io_ptr, "Enter boundary condition: ", 
			bc(1,1), no_indent);
      break;

    case 1: 
      bc(2,1) = get_int(io_ptr, "Enter boundary condition: ", 
			bc(2,1), no_indent);
      break;

    case 2: 
      bc(1,2) = get_int(io_ptr, "Enter boundary condition: ", 
			bc(1,2), no_indent);
      break;

    case 3:
      bc(2,2) = get_int(io_ptr, "Enter boundary condition: ", 
			bc(2,2), no_indent);
      break;

    case 4:
      printf("Boundary condition at r=0: %i, r=1: %i, s=0: %i, s=1: %i.\n",
	     bc(1,1), bc(2,1), bc(1,2), bc(2,2));
      replot = 0;
      break;

    case 5:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case 6:
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);
#undef ncom
}

void 
set_curve_label(input_output *io_ptr, generic_mapping *grid_ptr){
  int icom, quit, replot, plot_mode;
  const int level=1, no_indent=0;
  char prompt[80];

#include "set_curve_com.h"

  sprintf(prompt, "%s: set curve-value>", grid_ptr->grid_name);
  quit = 0;
  replot = 1;
  plot_mode = 35;
  PL_window(1, grid_ptr->xyab);
  do{

    if (replot){
/* plot the mapping */
      PL_erase(1);
      PL_start_plot(1);
      plot_generic_mapping( grid_ptr, plot_mode );
      PL_end_plot();
    }

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {

    case 0: 
      curve(1,1) = get_int( io_ptr, "Enter curve-value: ", curve(1,1), no_indent);
      break;

    case 1: 
      curve(2,1) = get_int( io_ptr, "Enter curve-value: ", curve(2,1), no_indent);
      break;

    case 2: 
      curve(1,2) = get_int( io_ptr, "Enter curve-value: ", curve(1,2), no_indent);
      break;

    case 3:
      curve(2,2) = get_int( io_ptr, "Enter curve-value: ", curve(2,2), no_indent);
      break;

    case 4:
      printf("curve at r=0: %i, r=1: %i, s=0: %i, s=1: %i.\n",
	     curve(1,1), curve(2,1), curve(1,2), curve(2,2));
      replot = 0;
      break;

    case 5:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case 6:
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);

#undef ncom
}



void 
set_grid_plot_mode( input_output *io_ptr, generic_mapping *grid_ptr ){

  char *prompt, *file_name;
  int level=3;
  const int save_command=1;
  int icom, quit, replot, col_flag;
  static int coordinate_axis=1;
  FILE *fp;

#include "set_grid_plot_com.h"

  prompt = "grid plot-mode>";

  quit = 0;
  replot = 0;
  PL_window(1, grid_ptr->xyab);
  do{

    if (replot){
/* plot the mapping */
      PL_erase(1);
      PL_start_plot(1);
      plot_generic_mapping( grid_ptr, grid_ptr->grid_plot_mode );
      PL_end_plot();
      replot = 0;
    }

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {
    case 0:
      if (!(grid_ptr->grid_plot_mode & 1))
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 1;
      else
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode & ~1;
      break;

    case 1:
      if (!(grid_ptr->grid_plot_mode & 2))
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 2;
      else
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode & ~2;
      break;

    case 2:
      if (!(grid_ptr->grid_plot_mode & 4))
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 4;
      else
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode & ~4;
      break;

    case 3:
      if (!(grid_ptr->grid_plot_mode & 8))
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 8;
      else
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode & ~8;
      break;

    case 4:
      if (!(grid_ptr->grid_plot_mode & 16))
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 16;
      else
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode & ~16;
      break;

    case 5:
      if (!(grid_ptr->grid_plot_mode & 32))
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 32;
      else
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode & ~32;
      break;

    case 6:
      if (!(grid_ptr->grid_plot_mode & 64))
/* turn the title on */
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 64;
      else
/* turn the title off */
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode & ~64;
      break;

    case 7:
/* toggle coordinate axis */
      coordinate_axis = !coordinate_axis;
      PL_scale( 1, coordinate_axis );
      break;

    case 8:
/* save postscript */
      if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
				 "mapping.ps", &file_name, 'w', 0, 
				 save_command) ) != NULL){
	col_flag = get_yes_no( io_ptr, "Do you want to save color "
				 "information>", level );
	PL_postscript( 1, fp, col_flag );
	fclose( fp );
      }
      replot = 0;
      break;


    case 9:
      printf("Grid Plotmode:\n");
      printf("Grid boundary: %s\n", (grid_ptr->grid_plot_mode & 1)? "on" : "off");
      printf("Grid directional arrows: %s\n", 
	     (grid_ptr->grid_plot_mode & 2)? "on" : "off");
      printf("Grid tickmarks: %s\n", (grid_ptr->grid_plot_mode & 4)? "on" : "off");
      printf("All gridpoints: %s\n", (grid_ptr->grid_plot_mode & 8)? "on" : "off");
      printf("Boundary condition: %s\n", 
	     (grid_ptr->grid_plot_mode & 16)? "on" : "off");
      printf("Curve label: %s\n", (grid_ptr->grid_plot_mode & 32)? "on" : "off");
      printf("\n");
/* don't redraw the plot */
      replot = 0;
      break;

    case 10:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case 11:
      quit = 1;
/* don't redraw the plot */
      replot = 0;
      break;

    default:
/* don't redraw the plot */
      replot = 0;
    }

  }
  while(!quit);
#undef ncom
}


void 
set_gridlines(input_output *io_ptr, generic_mapping *grid_ptr){
  char prompt[80];
  int icom, quit, replot;
  const int level=1, no_indent=0, plot_mode=7;

#include "set_gridlines_com.h"

  sprintf(prompt, "%s: number of gridlines>", grid_ptr->grid_name);
  quit = 0;
  replot = 1;
  PL_window(1, grid_ptr->xyab);
  do{

    if (replot){
/* plot the mapping */
      PL_erase(1);
      PL_start_plot(1);
      plot_generic_mapping( grid_ptr, plot_mode );
      PL_end_plot();
    }

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {

    case 0: 
      grid_ptr->r_points = 
	int_max(2, get_int(io_ptr, 
			   "Enter number of gridlines in r >= 2: ", 
			   grid_ptr->r_points,
			   no_indent));
      break;

    case 1:
      grid_ptr->s_points = 
	int_max(2, get_int(io_ptr, 
			   "Enter number of gridlines in s>= 2: ", grid_ptr->s_points,
			   no_indent));
      break;

    case 2:
      printf("Number of gridlines in r: %i and in s: %i.\n",
	     grid_ptr->r_points, grid_ptr->s_points);
      replot = 0;
      break;

    case 3:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case 4:
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);
#undef ncom
}


void 
reset_view_port( generic_mapping *map_ptr ){
  real x_pad, y_pad, x1, x2, x3, x4, y1, y2, y3, y4;

/* compute the transformed coordinates of the corners of the bounding box */
  x1 = map_ptr->scaling* (map_ptr->x_trans + 
			  map_ptr->x_min_0 * map_ptr->cos_theta -
			  map_ptr->y_min_0 * map_ptr->sin_theta);
  y1 = map_ptr->scaling* (map_ptr->y_trans + 
			  map_ptr->x_min_0 * map_ptr->sin_theta +
			  map_ptr->y_min_0 * map_ptr->cos_theta);

  x2 = map_ptr->scaling* (map_ptr->x_trans + 
			  map_ptr->x_max_0 * map_ptr->cos_theta -
			  map_ptr->y_min_0 * map_ptr->sin_theta);
  y2 = map_ptr->scaling* (map_ptr->y_trans + 
			  map_ptr->x_max_0 * map_ptr->sin_theta +
			  map_ptr->y_min_0 * map_ptr->cos_theta);

  x3 = map_ptr->scaling* (map_ptr->x_trans + 
			  map_ptr->x_max_0 * map_ptr->cos_theta -
			  map_ptr->y_max_0 * map_ptr->sin_theta);
  y3 = map_ptr->scaling* (map_ptr->y_trans + 
			  map_ptr->x_max_0 * map_ptr->sin_theta +
			  map_ptr->y_max_0 * map_ptr->cos_theta);

  x4 = map_ptr->scaling* (map_ptr->x_trans + 
			  map_ptr->x_min_0 * map_ptr->cos_theta -
			  map_ptr->y_max_0 * map_ptr->sin_theta);
  y4 = map_ptr->scaling* (map_ptr->y_trans + 
			  map_ptr->x_min_0 * map_ptr->sin_theta +
			  map_ptr->y_max_0 * map_ptr->cos_theta);

/* get the min and max from the transformed corners */
  map_ptr->x_min = real_min( x1, real_min( x2, real_min( x3, x4 ) ) );
  map_ptr->y_min = real_min( y1, real_min( y2, real_min( y3, y4 ) ) );
  
  map_ptr->x_max = real_max( x1, real_max( x2, real_max( x3, x4 ) ) );
  map_ptr->y_max = real_max( y1, real_max( y2, real_max( y3, y4 ) ) );
  
/* compute padding */  
  x_pad = 0.05*(map_ptr->x_max - map_ptr->x_min);
  y_pad = 0.05*(map_ptr->y_max - map_ptr->y_min);
  map_ptr->xyab[0][0] = map_ptr->x_min - x_pad;
  map_ptr->xyab[0][1] = map_ptr->x_max + x_pad;
  map_ptr->xyab[1][0] = map_ptr->y_min - y_pad;
  map_ptr->xyab[1][1] = map_ptr->y_max + y_pad;

}


void 
show_mapping_parameters( generic_mapping *grid_ptr ){

/* grid name */
  printf("Grid name: `%s'\n",grid_ptr->grid_name);

/* mapping type */
  printf("\t");
  switch (grid_ptr->grid_type){

  case 1:
    printf("cartesian mapping\n");
    break;
  case 4:
    printf("normal curve mapping\n");
    break;
  case 5:
    printf("linear interpolation mapping\n");
    break;
  case 6:
    printf("discrete point mapping\n");
    break;
#ifdef THE_GAR
  case 7:
    printf("Theodorsen-Garrick mapping\n");
    break;
#endif
  default:
    printf("unknown mapping\n");
  }

/* boundary conditions */
  printf("\tBoundary condition at r=0: %i, r=1: %i, s=0: %i, s=1: %i.\n",
	 bc(1,1), bc(2,1), bc(1,2), bc(2,2));

/* curve */
  printf("\tCurve label at r=0: %i, r=1: %i, s=0: %i, s=1: %i.\n",
	 curve(1,1), curve(2,1), curve(1,2), curve(2,2));

/* number of gridlines */
  printf("\tNumber of gridlines in r: %i and in s: %i.\n",
	 grid_ptr->r_points, grid_ptr->s_points);

/* one separator line */
  printf("\n");
}


generic_mapping *
get_grid_ptr( input_output *io_ptr, generic_mapping *first_m ){
  generic_mapping *grid_ptr;
  char **command;
  int i, icom, ncom, *save_on_copy=NULL, *argument=NULL, level=0, quit, n_grids;

  if (first_m == NULL) return NULL;

/* count the number of curves */
  n_grids = 0;
  for (grid_ptr = first_m; grid_ptr != NULL; grid_ptr = grid_ptr->next_m)
    n_grids++;

  ncom = n_grids + 1;
  if ((command = (char **) malloc( (ncom+1)*sizeof(char*) )) == NULL)
    printf("memory error in get_grid_ptr\n");

  i=0;
  for (grid_ptr=first_m; grid_ptr != NULL; grid_ptr = grid_ptr->next_m){
    command[i] = grid_ptr->grid_name;
    i++;
  }
  command[ncom-1] = "help";
  command[ncom] = "cancel";

  quit = 0;
  do{
    icom = get_command(io_ptr, "mapping name: ", command, ncom, level, 
		       save_on_copy, argument);

/* default is to quit */
    quit = 1;

    if (icom == -1)
/* not a valid command */
      quit = 0;
    else if (icom == ncom-1){
/* help */
      printf(
" o  Enter one of the mapping names if you wish to proceed.\n"
" o  Enter cancel if you do not wish to proceed. You will then return to the\n"
"    previous command level.\n");
      quit = 0;
    }

  }
  while (!quit);

/* return the pointer to the selected mapping */
  i=0;
  for (grid_ptr=first_m; grid_ptr != NULL; grid_ptr = grid_ptr->next_m){
    if (i==icom) break;
    i++;
  }

  free(command);

  return grid_ptr;

}

int 
forward_mapping( grid_point *gp_ptr, generic_mapping *map_ptr ){
  grid_point gp_old;
  int msg;

  gp_old.r = gp_ptr->r;
  gp_old.s = gp_ptr->s;

/* take care of any errors in the specific mapping routine and pass them on */
  if ( (msg=(*map_ptr->grid_mapping)( &gp_old, map_ptr->mapping_data_ptr )) != 0)
    return msg;

/* the transformation is done in the order: rotation - translation - scaling */
  gp_ptr->x = map_ptr->scaling* (map_ptr->x_trans + 
				 gp_old.x * map_ptr->cos_theta -
				 gp_old.y * map_ptr->sin_theta);

  gp_ptr->y = map_ptr->scaling* (map_ptr->y_trans + 
				 gp_old.x * map_ptr->sin_theta +
				 gp_old.y * map_ptr->cos_theta);

/* jacobian */
  gp_ptr->xr = map_ptr->scaling* (gp_old.xr * map_ptr->cos_theta -
				  gp_old.yr * map_ptr->sin_theta);

  gp_ptr->xs = map_ptr->scaling* (gp_old.xs * map_ptr->cos_theta -
				  gp_old.ys * map_ptr->sin_theta);
				 
  gp_ptr->yr = map_ptr->scaling* (gp_old.xr * map_ptr->sin_theta +
				  gp_old.yr * map_ptr->cos_theta);

  gp_ptr->ys = map_ptr->scaling* (gp_old.xs * map_ptr->sin_theta +
				  gp_old.ys * map_ptr->cos_theta);
				 
  return 0;
}


int 
backward_mapping( grid_point *gp_ptr, generic_mapping *map_ptr ){
  real x_new, y_new;

  x_new = gp_ptr->x;
  y_new = gp_ptr->y;

/* invert the rotation - translation - scaling */
  gp_ptr->x = (x_new/map_ptr->scaling - map_ptr->x_trans) * map_ptr->cos_theta +
    (y_new/map_ptr->scaling - map_ptr->y_trans) * map_ptr->sin_theta;

  gp_ptr->y =-(x_new/map_ptr->scaling - map_ptr->x_trans) * map_ptr->sin_theta +
    (y_new/map_ptr->scaling - map_ptr->y_trans) * map_ptr->cos_theta;

/* invert the un-transformed (x,y) coordinates */
  return (*map_ptr->inverse_grid_mapping)( gp_ptr, map_ptr->mapping_data_ptr );

}


void 
set_transformation(input_output *io_ptr, generic_mapping *map_ptr){
  char prompt[80];
  int icom, quit, replot;
  const int level=1, no_indent=0, plot_mode=1;
  const real eps=1.e-7;

#include "set_transformation_com.h"

  sprintf(prompt, "%s: transformation>", map_ptr->grid_name);
  quit = 0;
  replot = 1;
  do{

    if (replot){
/* update the view port */
      reset_view_port( map_ptr );
      PL_window(1, map_ptr->xyab);
/* plot the mapping */
      PL_erase(1);
      PL_start_plot(1);
      plot_generic_mapping( map_ptr, plot_mode );
      PL_end_plot();
    }

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {

    case 0: 
      map_ptr->theta_rot = get_real( io_ptr, "Angle of rotation (degrees): ", 
				    map_ptr->theta_rot, no_indent);
      map_ptr->cos_theta = cos( ((double)M_PI)/180.0 * map_ptr->theta_rot );
      map_ptr->sin_theta = sin( ((double)M_PI)/180.0 * map_ptr->theta_rot );
      break;

    case 1: 
      map_ptr->x_trans = get_real( io_ptr, "Horizontal translation: ", 
				    map_ptr->x_trans, no_indent);
      break;

    case 2: 
      map_ptr->y_trans = get_real( io_ptr, "Vertical translation: ", 
				    map_ptr->y_trans, no_indent);
      break;

    case 3:
      map_ptr->scaling = real_max(eps, 
				  get_real( io_ptr, "Scaling (>0): ", 
				    map_ptr->scaling, no_indent) );
      break;

    case 4:
      printf("Rotation angle = %f\n"
	     "Horizontal translation = %f\n"
	     "Vertical translation = %f\n"
	     "Scaling factor = %f\n", map_ptr->theta_rot, map_ptr->x_trans, 
	     map_ptr->y_trans, map_ptr->scaling);
      replot = 0;
      break;

    case 5:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case 6:
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);
#undef ncom
}



void 
map_transform(real x_old, real y_old, real *x_new, real *y_new, 
	      generic_mapping *map_ptr ){

/* the transformation is done in the order: rotation - translation - scaling */
  *x_new = map_ptr->scaling* (map_ptr->x_trans + 
				x_old * map_ptr->cos_theta -
				y_old * map_ptr->sin_theta);

  *y_new = map_ptr->scaling* (map_ptr->y_trans + 
				x_old * map_ptr->sin_theta +
				y_old * map_ptr->cos_theta);

}


void 
inv_map_transform(real x_old, real y_old, real *x_new, real *y_new, 
		  generic_mapping *map_ptr ){

/* invert the rotation - translation - scaling */
  *x_new = (x_old/map_ptr->scaling - map_ptr->x_trans) * map_ptr->cos_theta +
    (y_old/map_ptr->scaling - map_ptr->y_trans) * map_ptr->sin_theta;

  *y_new =-(x_old/map_ptr->scaling - map_ptr->x_trans) * map_ptr->sin_theta +
    (y_old/map_ptr->scaling - map_ptr->y_trans) * map_ptr->cos_theta;

}


void 
set_mapping( input_output *io_ptr, generic_mapping *map_ptr, generic_curve *first_c){
  if (map_ptr == NULL || map_ptr->set_mapping == NULL){
    printf("ERROR: set_mapping was called with inconsistent pointers.\n");
  }
  else
    (*map_ptr->set_mapping)( io_ptr, map_ptr, first_c );
}


void 
get_grid_quality( grid_quality *qp, generic_mapping *map_ptr ){
  grid_point gp1, gp2;
  real ar, ortho, theta, dsdr, dsds, xrxs, dsds1, dsds2, xs1xs2
    , dsdr1, dsdr2, xr1xr2, size_r, size_s, area, theta_r, theta_s
      , change_r, change_s, r_step, s_step;
  int i, j;

/* grid sizes */
  r_step = 1.0/( map_ptr->r_points - 1);
  s_step = 1.0/( map_ptr->s_points - 1);

/* aspect ratio and orthogonality */
  qp->max_ar = 0.0;       qp->i_ar = 0;    qp->j_ar = 0; 
  qp->min_ortho = 1.e10;  qp->i_ortho = 0; qp->j_ortho = 0;
  qp->min_size_r = 1.e10; qp->i_msr = 0;   qp->j_msr = 0;
  qp->min_size_s = 1.e10; qp->i_mss = 0;   qp->j_mss = 0; 
  qp->min_area = 1.e10;   qp->i_area = 0;  qp->j_area = 0;  

  for (i=1; i<= map_ptr->r_points; i++)
    for (j=1; j<= map_ptr->s_points; j++){
      gp1.r = (i-1) * r_step;
      gp1.s = (j-1) * s_step;
      forward_mapping( &gp1, map_ptr );

/* find the angle between the r and s lines */
      dsdr = sqrt( gp1.xr * gp1.xr + gp1.yr * gp1.yr );
      dsds = sqrt( gp1.xs * gp1.xs + gp1.ys * gp1.ys );
      xrxs = (gp1.xr * gp1.xs + gp1.yr * gp1.ys) / (dsdr * dsds);
      theta = acos( xrxs );
      ortho = fabs( sin(theta) );
      ar = dsds * fabs( sin(theta) ) * s_step / ( dsdr * r_step );
      if (ar < 1.0) ar = 1.0/real_max( NEWTON_EPS, ar );

/* size_r and size_s */
      size_r = dsdr * r_step;
      size_s = dsds * s_step;

/* find the area of the cell */
      area = dsds * fabs( sin(theta) ) * s_step * dsdr * r_step;

/* check if this cell is extreme */
      if (ar > qp->max_ar){
	qp->max_ar = ar;
	qp->i_ar = i; qp->j_ar = j;
      }

      if (ortho < qp->min_ortho){
	qp->min_ortho = ortho;
	qp->i_ortho = i; qp->j_ortho = j;
      }

      if (area < qp->min_area){
	qp->min_area = area;
	qp->i_area = i;
	qp->j_area = j;
      }

      if (size_r < qp->min_size_r){
	qp->min_size_r = size_r;
	qp->i_msr = i;
	qp->j_msr = j;
      }

      if (size_s < qp->min_size_s){
	qp->min_size_s = size_s;
	qp->i_mss = i;
	qp->j_mss = j;
      }
      
    }

  printf("\n");
  printf("Grid quality for grid `%s':\n\n", map_ptr->grid_name);

  printf("BOX: Min cell area: %f occured at i=%i, j=%i.\n", qp->min_area, 
	 qp->i_area, qp->j_area);
  printf("ASTERIX: Min grid size in the r-direction: %f occured at i=%i, j=%i.\n", 
	 qp->min_size_r, qp->i_msr, qp->j_msr);
  printf("CROSS: Min grid size in the s-direction: %f occured at i=%i, j=%i.\n\n", 
	 qp->min_size_s, qp->i_mss, qp->j_mss);


  printf("FILLED CIRCLE: Max aspect ratio: %f occured at i=%i, j=%i.\n", 
	 qp->max_ar, qp->i_ar, qp->j_ar);
  printf("FILLED BOX: Min orthogonality: %f degrees occured at i=%i, j=%i.\n\n", 
	 asin(qp->min_ortho) * 180.0/((double)M_PI), qp->i_ortho, qp->j_ortho);

/* change in grid size and direction along constant-r-lines */
  qp->max_theta_r = -1.0; qp->i_theta_r=0;  qp->j_theta_r=0;
  qp->max_change_r = 0.5; qp->i_change_r=0; qp->j_change_r=0;
  for (i=1; i<= map_ptr->r_points; i++){
    j = 1;
    gp1.r = (i-1) * r_step;
    gp1.s = (j-1) * s_step;
    forward_mapping( &gp1, map_ptr );
    dsds1 = sqrt( gp1.xs * gp1.xs + gp1.ys * gp1.ys );

    for (j=2; j<= map_ptr->s_points; j++){
      gp2.r = (i-1) * r_step;
      gp2.s = (j-1) * s_step;
      forward_mapping( &gp2, map_ptr );
      dsds2 = sqrt( gp2.xs * gp2.xs + gp2.ys * gp2.ys );

/* find the angle between the s lines */
      xs1xs2 = real_min( 1.0, (gp1.xs * gp2.xs + gp1.ys * gp2.ys) / (dsds1 * dsds2) );
      theta_r = acos( xs1xs2 );

/* change in size */
      change_r = dsds2 / dsds1;
      if (change_r < 1.0) change_r = 1.0/real_max(NEWTON_EPS, change_r);

/* check if this cell is extreme */
      if (theta_r > qp->max_theta_r){
	qp->max_theta_r = theta_r;
	qp->i_theta_r = i; qp->j_theta_r = j;
      }

      if (change_r > qp->max_change_r){
	qp->max_change_r = change_r;
	qp->i_change_r = i; qp->j_change_r = j;
      }

/* copy gp2 to gp1 and dsds2 to dsds1 */
      gp1 = gp2;
      dsds1 = dsds2;
    }
  }

  printf("TRIANGLE: Max change in direction along constant-r-lines: %f degrees "
	 "occured at i=%i, j=%i.\n"
	 , qp->max_theta_r * 180.0/((double)M_PI), qp->i_theta_r, qp->j_theta_r);
  printf("FILLED TRIANGLE: Max relative increase in grid size along constant-r-lines: "
	 "%f occured at i=%i, j=%i.\n\n", 
	 qp->max_change_r - 1.0, qp->i_change_r, qp->j_change_r);

/* change in grid size and direction along constant-s-lines */
  qp->max_theta_s = -1.0; qp->i_theta_s=0;  qp->j_theta_s=0;
  qp->max_change_s = 0.5; qp->i_change_s=0; qp->j_change_s=0;
  for (j=1; j<= map_ptr->s_points; j++){
    i = 1;
    gp1.r = (i-1) * r_step;
    gp1.s = (j-1) * s_step;
    forward_mapping( &gp1, map_ptr );
    dsdr1 = sqrt( gp1.xr * gp1.xr + gp1.yr * gp1.yr );

    for (i=2; i<= map_ptr->r_points; i++){
      gp2.r = (i-1) * r_step;
      gp2.s = (j-1) * s_step;
      forward_mapping( &gp2, map_ptr );
      dsdr2 = sqrt( gp2.xr * gp2.xr + gp2.yr * gp2.yr );

/* find the angle between the s lines */
      xr1xr2 = real_min( 1.0, (gp1.xr * gp2.xr + gp1.yr * gp2.yr) / (dsdr1 * dsdr2) );
      theta_s = acos( xr1xr2 );

/* change_s in size */
      change_s = dsdr2 / dsdr1;
      if (change_s < 1.0) change_s = 1.0/real_max(NEWTON_EPS, change_s);

/* check if this cell is extreme */
      if (theta_s > qp->max_theta_s){
	qp->max_theta_s = theta_s;
	qp->i_theta_s = i; qp->j_theta_s = j;
      }

      if (change_s > qp->max_change_s){
	qp->max_change_s = change_s;
	qp->i_change_s = i; qp->j_change_s = j;
      }

/* copy gp2 to gp1 and dsdr2 to dsdr1 */
      gp1 = gp2;
      dsdr1 = dsdr2;
    }
  }

  printf("INVERTED TRIANGLE: Max change in direction along constant-s-lines: %f "
	 "degrees occured at i=%i, j=%i.\n"
	 , qp->max_theta_s * 180.0/((double)M_PI), qp->i_theta_s, qp->j_theta_s);
  printf("FILLED INVERTED TRIANGLE: Max relative increase in grid size along "
	 "constant-s-lines: %f occured at i=%i, j=%i.\n", 
	 qp->max_change_s - 1.0, qp->i_change_s, qp->j_change_s);

  printf("\n");

}
