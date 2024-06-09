#include "mappings.h"

/* linear interpolation between two curves */

/* private member function declarations */
static void compute_linear_interp(linear_interp_data *info, 
				  generic_mapping *grid_ptr );
static void 
write_linear_interp( int32 dir, void *linear_interp_data_ptr);
static void set_linear_interp_mapping(input_output *io_ptr, 
				      generic_mapping *grid_ptr, 
				  generic_curve *first_c);
static int linear_interp_mapping( grid_point *gp_ptr, void *data_ptr);
static void *cleanup_linear_interp( void *data_ptr );
static void *copy_linear_interp( void *mapping_data_ptr );
/* end */

static void set_linear_interp_mapping(input_output *io_ptr, 
				      generic_mapping *grid_ptr, 
				      generic_curve *first_c){
  char prompt[80];
  int icom, quit, replot;
  int level=2;
  const int no_indent=0;
  linear_interp_data *info;

#include "linear_interp_grid_com.h"

  sprintf(prompt, "%s: linear-interpolation mapping>", grid_ptr->grid_name);

/* cast the data pointer to the right type */
  info = (linear_interp_data *) grid_ptr->mapping_data_ptr;

  quit = 0;
  replot = 1;
/* recompute the view port because it is possible that the curve has been changed */
  compute_linear_interp( info, grid_ptr );
  PL_window(1, grid_ptr->xyab);
  do{

    if (replot){
/* plot the mapping */
      PL_erase(1);
      PL_start_plot(1);
      plot_generic_mapping( grid_ptr, grid_ptr->grid_plot_mode );
      PL_end_plot();
    }

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {

    case 0:
/* modify the s=0 curve */
      set_curve( io_ptr, info->curve_s0_ptr );
      compute_linear_interp( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 1:
/* flip the s=0 curve parametrization */
      info->flip_s0 = !info->flip_s0;
      compute_linear_interp( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 2:
/* stretch-s=0-curve */
      if (info->curve_s0_stretch == NULL){
	info->curve_s0_stretch = choose_stretching( io_ptr, info->curve_s0_ptr, 
						   &grid_ptr->r_points );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	if ( set_stretching( io_ptr, info->curve_s0_stretch, info->curve_s0_ptr, 
			     &grid_ptr->r_points ) ){
	  grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->curve_s0_stretch ) );
	  replot = 0;
	}
      }
      break;

    case 3:
/* don't stretch-s=0-curve */
      if (info->curve_s0_stretch != NULL){	
	info->curve_s0_stretch = 
	  delete_generic_stretching( info->curve_s0_stretch );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	printf("There is no stretching in the r-direction.\n");
	replot = 0;
      }
      break;

    case 4:
/* modify s=1 curve */
      set_curve( io_ptr, info->curve_s1_ptr );
      compute_linear_interp( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 5:
/* flip the s=1 curve parametrization */
      info->flip_s1 = !info->flip_s1;
      compute_linear_interp( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 6:
/* stretch s=1 curve */
      if (info->curve_s1_stretch == NULL){
	info->curve_s1_stretch = choose_stretching( io_ptr, info->curve_s1_ptr, 
						   &grid_ptr->r_points );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	if ( set_stretching( io_ptr, info->curve_s1_stretch, info->curve_s1_ptr, 
			     &grid_ptr->r_points) ){
	  grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->curve_s1_stretch ) );
	  replot = 0;
	}
      }
      break;

    case 7:
/* don't stretch s=1 curve */
      if (info->curve_s1_stretch != NULL){
	info->curve_s1_stretch = delete_generic_stretching( info->curve_s1_stretch );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	printf("There is no stretching in the r-direction.\n");
	replot = 0;
      }
      break;

    case 8:
/* s-stretching */
      if (info->normal_stretch == NULL){
	info->normal_stretch = choose_stretching( io_ptr, NULL, 
						 &grid_ptr->s_points );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	if ( set_stretching( io_ptr, info->normal_stretch, NULL, 
			    &grid_ptr->s_points) ){
	  grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->normal_stretch ) );
	  replot = 0;
	}
      }
      break;

    case 9:
/* no s-stretching */
      if (info->normal_stretch != NULL){
	info->normal_stretch = delete_generic_stretching( info->normal_stretch );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	printf("There is no stretching in the s-direction.\n");
	replot = 0;
      }
      break;

    case 10:
/* r-grid lines */
      grid_ptr->r_points = 
	int_max(2, get_int(io_ptr, 
			   "Enter number of gridlines in r >= 2: ", 
			   grid_ptr->r_points,
			   no_indent));

      grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      break;

    case 11:
/* s-grid lines */
      grid_ptr->s_points = 
	int_max(2, get_int(io_ptr, 
			   "Enter number of gridlines in s>= 2: ", grid_ptr->s_points,
			   no_indent));

      grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      break;

    case 12:
/* plot mode */
      set_grid_plot_mode( io_ptr, grid_ptr );
      replot = 0;
      break;

    case 13:
/* show parameters */
      printf("Properties if the s=0 curve:\n");
      show_curve_parameters( info->curve_s0_ptr );
      printf("%s.\n", stretching_name( info->curve_s0_stretch ) );
      printf("\n");
      printf("Properties if the s=1 curve:\n");
      show_curve_parameters( info->curve_s1_ptr );
      printf("%s.\n", stretching_name( info->curve_s1_stretch ) );
      printf("\n");
      show_mapping_parameters( grid_ptr );
      printf("%s between the boundary curves.\n", 
	     stretching_name( info->normal_stretch ) );
      printf("\n");

      replot = 0;
      break;

    case 14:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case 15:
/* quit */
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);

}


static int linear_interp_mapping( grid_point *gp_ptr, void *data_ptr){
  linear_interp_data *info;
  curve_point cp0, cp1;
  real not_used, s_stretch, s_prime, r0_prime, r1_prime, r0, r1;

/* cast the data to the right type */
  info = (linear_interp_data *) data_ptr;

/* flip the parameter? */
  if (info->flip_s0)
    r0 = 1.0 - gp_ptr->r;
  else
    r0 = gp_ptr->r;

  if (info->flip_s1)
    r1 = 1.0 - gp_ptr->r;
  else
    r1 = gp_ptr->r;

/* stretching */
  uniform_to_stretch( r0, info->curve_s0_stretch, &(cp0.r), &r0_prime, &not_used );

  uniform_to_stretch( r1, info->curve_s1_stretch, &(cp1.r), &r1_prime, &not_used );

  uniform_to_stretch( gp_ptr->s, info->normal_stretch, &s_stretch, &s_prime, 
		     &not_used );

/* call the curve routines */
  curve_function( &cp0, info->curve_s0_ptr );
  if (info->flip_s0){
    cp0.xr = - cp0.xr;
    cp0.yr = - cp0.yr;
  }

  curve_function( &cp1, info->curve_s1_ptr );
  if (info->flip_s1){
    cp1.xr = - cp1.xr;
    cp1.yr = - cp1.yr;
  }

/* Cartesian coordinates */
  gp_ptr->x = cp0.x * (1.0 - s_stretch) + cp1.x * s_stretch;
  gp_ptr->y = cp0.y * (1.0 - s_stretch) + cp1.y * s_stretch;

/* Jacobian */
  gp_ptr->xr = cp0.xr * r0_prime * (1.0 - s_stretch) + cp1.xr * r1_prime * s_stretch;
  gp_ptr->xs = (-cp0.x + cp1.x) * s_prime;
  gp_ptr->yr = cp0.yr * r0_prime * (1.0 - s_stretch) + cp1.yr * r1_prime * s_stretch;
  gp_ptr->ys = (-cp0.y + cp1.y) * s_prime;

  return 0;
}


static void compute_linear_interp(linear_interp_data *info, 
				  generic_mapping *grid_ptr ){
  grid_point gp;

/* the grid is periodic if both curves are periodic */
  grid_ptr->r_period = (info->curve_s0_ptr->periodic && 
			info->curve_s1_ptr->periodic);
  grid_ptr->s_period = 0;

/* compute bounding box */
  grid_ptr->x_min_0 = 1.e10;
  grid_ptr->x_max_0 = -1.e10;
  grid_ptr->y_min_0 = 1.e10;
  grid_ptr->y_max_0 = -1.e10;
/* s=0 */
  gp.s = 0.0;
  for (gp.r=0.0; gp.r<=1.05; gp.r+=0.1){
    linear_interp_mapping( &gp, grid_ptr->mapping_data_ptr );
    grid_ptr->x_min_0 = real_min( gp.x, grid_ptr->x_min_0 );
    grid_ptr->x_max_0 = real_max( gp.x, grid_ptr->x_max_0 );
    grid_ptr->y_min_0 = real_min( gp.y, grid_ptr->y_min_0 );
    grid_ptr->y_max_0 = real_max( gp.y, grid_ptr->y_max_0 );
  }
/* s=1 */
  gp.s = 1.0;
  for (gp.r=0.0; gp.r<=1.05; gp.r+=0.1){
    linear_interp_mapping( &gp, grid_ptr->mapping_data_ptr );
    grid_ptr->x_min_0 = real_min( gp.x, grid_ptr->x_min_0 );
    grid_ptr->x_max_0 = real_max( gp.x, grid_ptr->x_max_0 );
    grid_ptr->y_min_0 = real_min( gp.y, grid_ptr->y_min_0 );
    grid_ptr->y_max_0 = real_max( gp.y, grid_ptr->y_max_0 );
  }

  reset_view_port( grid_ptr );
}


static void *cleanup_linear_interp( void *data_ptr ){
  linear_interp_data *info;

/* cast the data to the right type */
  info = (linear_interp_data *) data_ptr;

/* flag that the curves are no longer used */
  info->curve_s0_ptr->used_by_grid--;
  info->curve_s1_ptr->used_by_grid--;

  delete_generic_stretching( info->curve_s0_stretch );
  delete_generic_stretching( info->curve_s1_stretch );
  delete_generic_stretching( info->normal_stretch );

/* cleanup the struct itself */
  free(info);

  return NULL;
}


static void 
write_linear_interp( int32 dir, void *linear_interp_data_ptr){
  linear_interp_data *info;

  int32 stretch_dir;

  info = (linear_interp_data*) linear_interp_data_ptr;

/* save the names of the curves */
  hput_string( info->curve_s0_ptr->curve_name, "curve_0_name", dir );
  hput_string( info->curve_s1_ptr->curve_name, "curve_1_name", dir );

  hput_int( info->flip_s0, "flip_s0", dir );
  hput_int( info->flip_s1, "flip_s1", dir );

/* stretching */
  if ((stretch_dir = create_dir("s0-stretching", 
				"s0-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( stretch_dir, info->curve_s0_stretch );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }
  else
    printf("Error: write_linear_interp: unable to create a directory for "
	   "`s0-stretching'.\n");

  if ((stretch_dir = create_dir("s1-stretching", 
				"s1-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( stretch_dir, info->curve_s1_stretch );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }
  else
    printf("Error: write_linear_interp: unable to create a directory for "
	   "`s1-stretching'.\n");

  if ((stretch_dir = create_dir("normal-stretching", 
				"normal-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( stretch_dir, info->normal_stretch );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }
  else
    printf("Error: write_linear_interp: unable to create a directory for "
	   "`normal-stretching'.\n");
}


linear_interp_data *
init_linear_interp(input_output *io_ptr, generic_mapping *grid_ptr,
		   generic_curve *curve_s0_ptr, generic_curve *curve_s1_ptr){
  linear_interp_data *info;
  
/* allocate space for the linear interpolation grid */
  if ((info = (linear_interp_data *) malloc( sizeof(linear_interp_data) )) == NULL)
    printf("memory error in init_linear_interp\n");

/* assign the pointers to the default curves */
  info->curve_s0_ptr = curve_s0_ptr;
  info->curve_s1_ptr = curve_s1_ptr;
  info->curve_s0_ptr->used_by_grid++;
  info->curve_s1_ptr->used_by_grid++;

/* initially not flipped */
  info->flip_s0 = 0;
  info->flip_s1 = 0;

/* no stretching */
  info->curve_s0_stretch = NULL;
  info->curve_s1_stretch = NULL;
  info->normal_stretch = NULL;

/* set the grid type and the periodicity flag */
  grid_ptr->grid_type= 5; /* linear interpolation grid */

/* set the number of grid lines */
  grid_ptr->r_points = 1.0/real_min(info->curve_s0_ptr->r_step,
				    info->curve_s1_ptr->r_step );
  grid_ptr->s_points = 10;

/* save the pointers to the data and the functions for
   evaluating the mapping */
  grid_ptr->grid_mapping = linear_interp_mapping;
  grid_ptr->inverse_known = 0;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->mapping_data_ptr = info;
  grid_ptr->cleanup_mapping_data = cleanup_linear_interp;
  grid_ptr->set_mapping = set_linear_interp_mapping;
  grid_ptr->write_specific_mapping = write_linear_interp;
  grid_ptr->copy_mapping_data = copy_linear_interp;

  compute_linear_interp( info, grid_ptr );

  return info;
}

int 
read_linear_interp( int32 dir, generic_mapping *grid_ptr, generic_curve *first_c ){
  linear_interp_data *info;

  char *curve_name;
  int32 stretch_dir;

  info = (linear_interp_data *) malloc( sizeof(linear_interp_data) );

/* read the name of the s=0 curve*/
  curve_name = hget_string( "curve_0_name", dir );

/* look for the s=0 curve in the curve list */
  for (info->curve_s0_ptr = first_c;
       info->curve_s0_ptr != NULL &&
       strcmp(info->curve_s0_ptr->curve_name, curve_name) != 0;
       info->curve_s0_ptr = info->curve_s0_ptr->next);
/* was the curve found? */
  if (info->curve_s0_ptr == NULL){
    printf("ERROR: read_linear_interp. Could not find curve with name: %s\n",
           curve_name);
    free(curve_name);
    free(info);
    return FALSE;
  }
  free( curve_name );

/* read the name of the s=1 curve*/
  curve_name = hget_string( "curve_1_name", dir );

/* look for the s=1 curve in the curve list */
  for (info->curve_s1_ptr = first_c;
       info->curve_s1_ptr != NULL &&
       strcmp(info->curve_s1_ptr->curve_name, curve_name) != 0;
       info->curve_s1_ptr = info->curve_s1_ptr->next);
/* was the curve found? */
  if (info->curve_s1_ptr == NULL){
    printf("ERROR: read_linear_interp. Could not find curve with name: %s\n",
           curve_name);
    free(curve_name);
    free(info);
    return FALSE;
  }
  free( curve_name );

  hget_int( &(info->flip_s0), "flip_s0", dir );
  hget_int( &(info->flip_s1), "flip_s1", dir);

/* stretching */
  if ((stretch_dir = locate_dir("s0-stretching", dir)) != -1){
/* save the generic and specific curve data */
    info->curve_s0_stretch = read_stretching( stretch_dir );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }

  if ((stretch_dir = locate_dir("s1-stretching", dir)) != -1){
/* save the generic and specific curve data */
    info->curve_s1_stretch = read_stretching( stretch_dir );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }

  if ((stretch_dir = locate_dir("normal-stretching", dir)) != -1){
/* save the generic and specific curve data */
    info->normal_stretch = read_stretching( stretch_dir );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }

  grid_ptr->mapping_data_ptr = (void *) info;
  grid_ptr->grid_mapping = linear_interp_mapping;
  grid_ptr->inverse_known = 0;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->cleanup_mapping_data = cleanup_linear_interp;
  grid_ptr->set_mapping = set_linear_interp_mapping;
  grid_ptr->write_specific_mapping = write_linear_interp;
  grid_ptr->copy_mapping_data = copy_linear_interp;

  return OK;
}


static void *copy_linear_interp( void *mapping_data_ptr ){
  linear_interp_data *info, *old_info;
  
  info = (linear_interp_data *) malloc( sizeof(linear_interp_data) );

  old_info = (linear_interp_data *)mapping_data_ptr;

/* only make a shallow copy of the curves. */
  info->curve_s0_ptr = old_info->curve_s0_ptr;
  info->curve_s1_ptr = old_info->curve_s1_ptr;

/* remember that these curves are now used by both mappings */
  info->curve_s0_ptr->used_by_grid++;
  info->curve_s1_ptr->used_by_grid++;

  info->flip_s0 = old_info->flip_s0;
  info->flip_s1 = old_info->flip_s1;

/* make a (deep) copy of the stretching functions */
  info->curve_s0_stretch = copy_stretching( old_info->curve_s0_stretch );
  info->curve_s1_stretch = copy_stretching( old_info->curve_s1_stretch );
  info->normal_stretch   = copy_stretching( old_info->normal_stretch   );

  return (void *)info;
}
