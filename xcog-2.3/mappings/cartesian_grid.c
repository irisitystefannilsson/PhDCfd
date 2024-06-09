#include "mappings.h"

/* mappings for cartesian grids */

/* private member functions */
static void compute_cartesian_mapping( cartesian_grid_info *info, 
				      generic_mapping *grid_ptr );
static void set_cartesian_mapping(input_output *io_ptr, generic_mapping *grid_ptr, 
				  generic_curve *first_c);
static int cartesian_mapping( grid_point *gp_ptr, void *info);
static int inverse_cartesian_mapping( grid_point *gp_ptr, void *info);
static void *cleanup_cartesian_mapping( void *data_ptr );
static void 
write_cartesian_mapping( int32 dir, void *cartesian_data_ptr );
static void *copy_cartesian_mapping( void *mapping_data_ptr );
/* end private member function declaration */

static void set_cartesian_mapping(input_output *io_ptr, generic_mapping *grid_ptr, 
				  generic_curve *first_c){
  char prompt[80];
  int icom, quit, i, replot;
  cartesian_grid_info *info;
  int level=2, no_indent=0;
  real xmin, ymin, xmax, ymax;

#include "cartesian_grid_com.h"

  sprintf(prompt, "%s: cartesian mapping>", grid_ptr->grid_name);

/* cast the data pointer to the right type */
  info = (cartesian_grid_info *) grid_ptr->mapping_data_ptr;

  quit = 0;
  replot = 1;
  PL_window(1, grid_ptr->xyab);
  do{

    if (replot){
/* plot the mapping */
      PL_erase(1);
      PL_start_plot(1);
      plot_generic_mapping( grid_ptr, grid_ptr->grid_plot_mode );
      PL_end_plot();
    }

    icom = get_command(io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {

    case 0: 
/* global mapping */
      map_transform( info->xmin, info->ymin, &xmin, &ymin, grid_ptr );
      info->xmin = get_real( io_ptr, "xmin: ", xmin, no_indent);
/* inverse global mapping */
      inv_map_transform( info->xmin, ymin, &(info->xmin), &(info->ymin), grid_ptr );

      compute_cartesian_mapping( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 1:
/* global mapping */
      map_transform( info->xmax, info->ymax, &xmax, &ymax, grid_ptr );
      info->xmax = get_real( io_ptr, "xmax: ", xmax, no_indent);
/* inverse global mapping */
      inv_map_transform( info->xmax, ymax, &(info->xmax), &(info->ymax), grid_ptr );

      compute_cartesian_mapping( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 2:
/* global mapping */
      map_transform( info->xmin, info->ymin, &xmin, &ymin, grid_ptr );
      info->ymin = get_real( io_ptr, "ymin: ", ymin, no_indent);
/* inverse global mapping */
      inv_map_transform( xmin, info->ymin, &(info->xmin), &(info->ymin), grid_ptr );

      compute_cartesian_mapping( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 3:
/* global mapping */
      map_transform( info->xmax, info->ymax, &xmax, &ymax, grid_ptr );
      info->ymax = get_real( io_ptr, "ymax: ", ymax, no_indent);
/* inverse global mapping */
      inv_map_transform( xmax, info->ymax, &(info->xmax), &(info->ymax), grid_ptr );

      compute_cartesian_mapping( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 4:
/* r-stretching */
      if (info->r_stretch == NULL){
	info->r_stretch = choose_stretching( io_ptr, NULL, &grid_ptr->r_points );
/* the inverse is no longer valid */
	grid_ptr->inverse_known = 0; grid_ptr->inverse_grid_mapping = NULL;
/* update the plot mode */
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	if ( set_stretching( io_ptr, info->r_stretch, NULL, &grid_ptr->r_points) ){
	  grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->r_stretch ) );
	  replot = 0;
	}
      }
      break;

    case 5:
/* no-r-stretching */
      if (info->r_stretch != NULL){
	info->r_stretch = delete_generic_stretching( info->r_stretch );
/* the inverse is valid if also the s-stretching is off */
	if (info->s_stretch == NULL){
	  grid_ptr->inverse_known = 1;
	  grid_ptr->inverse_grid_mapping = inverse_cartesian_mapping;
	}
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	printf("There is no stretching in the s-direction.\n");
	replot = 0;
      }
      break;

    case 6:
/* s-stretching */
      if (info->s_stretch == NULL){
	info->s_stretch = choose_stretching( io_ptr, NULL, &grid_ptr->s_points );
/* the inverse is no longer valid */
	grid_ptr->inverse_known = 0; grid_ptr->inverse_grid_mapping = NULL;
/* update the plot mode */
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	if ( set_stretching( io_ptr, info->s_stretch, NULL, &grid_ptr->s_points) ){
	  grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->s_stretch ) );
	  replot = 0;
	}
      }
      break;

    case 7:
/* no-s-stretching */
      if (info->s_stretch != NULL){
	info->s_stretch = delete_generic_stretching( info->s_stretch );
/* the inverse is valid if also the r-stretching is off */
	if (info->r_stretch == NULL){
	  grid_ptr->inverse_known = 1;
	  grid_ptr->inverse_grid_mapping = inverse_cartesian_mapping;
	}
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	printf("There is no stretching in the s-direction.\n");
	replot = 0;
      }
      break;

    case 8:
/* r-lines */
      grid_ptr->r_points = 
	int_max(2, get_int(io_ptr, 
			   "Enter number of gridlines in r >= 2: ", 
			   grid_ptr->r_points,
			   no_indent));

      grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      break;

    case 9:
/* s-lines */
      grid_ptr->s_points = 
	int_max(2, get_int(io_ptr, 
			   "Enter number of gridlines in s>= 2: ", grid_ptr->s_points,
			   no_indent));

      grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      break;

    case 10:
/* plot-mode */
      set_grid_plot_mode( io_ptr, grid_ptr );
      replot = 0;
      break;

    case 11:
/* show */
      map_transform( info->xmin, info->ymin, &xmin, &ymin, grid_ptr );
      map_transform( info->xmax, info->ymax, &xmax, &ymax, grid_ptr );
      printf("xmin: %e, ymin: %e,\n"
	     "xmax: %e, ymax: %e.\n", xmin, ymin, xmax, ymax);
      show_mapping_parameters( grid_ptr );
      printf("%s in the r-direction and %s in the s-direction.\n", 
	     stretching_name( info->r_stretch ), 
	     stretching_name( info->s_stretch ) );
      replot = 0;
      break;

    case 12:
/* help */
      while ( (icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				  command, ncom, level+1, NULL, NULL)) == -1 );
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case 13:
/* exit */
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);

}

static int cartesian_mapping( grid_point *gp_ptr, void *data_ptr){
  cartesian_grid_info *info;
  real r_stretch, r_ratio, s_stretch, s_ratio, dummy;

/* cast the data to the right type */
  info = (cartesian_grid_info *) data_ptr;

/* r-stretching */
  uniform_to_stretch( gp_ptr->r, info->r_stretch, &r_stretch, 
		     &r_ratio, &dummy );

/* s-stretching */
  uniform_to_stretch( gp_ptr->s, info->s_stretch, &s_stretch, 
		     &s_ratio, &dummy );

  gp_ptr->x= info->xmin +  r_stretch * (info->xmax - info->xmin);
  gp_ptr->y= info->ymin +  s_stretch * (info->ymax - info->ymin);

  gp_ptr->xr= (info->xmax - info->xmin) * r_ratio;
  gp_ptr->xs= 0.;
  gp_ptr->yr= 0.;
  gp_ptr->ys= (info->ymax - info->ymin) * s_ratio;

  return 0;
}


static int inverse_cartesian_mapping( grid_point *gp_ptr, void *data_ptr){
  cartesian_grid_info *info;

/* cast the data to the right type */
  info = (cartesian_grid_info *) data_ptr;

/* we need to take stretching into consideration here */

  gp_ptr->r = (gp_ptr->x - info->xmin) / (info->xmax - info->xmin);
  gp_ptr->s = (gp_ptr->y - info->ymin) / (info->ymax - info->ymin);
  
  return 0;
}


static void compute_cartesian_mapping( cartesian_grid_info *info, 
			       generic_mapping *grid_ptr ){
/* set the bounding box */
  grid_ptr->x_min_0 = info->xmin;
  grid_ptr->x_max_0 = info->xmax;
  grid_ptr->y_min_0 = info->ymin;
  grid_ptr->y_max_0 = info->ymax;

  reset_view_port( grid_ptr );
}


static void 
write_cartesian_mapping( int32 dir, void *cartesian_data_ptr ){
  cartesian_grid_info *info;
  int32 r_stretch_dir, s_stretch_dir;

  info = (cartesian_grid_info *) cartesian_data_ptr;

  hput_real( info->xmin, "xmin", dir );
  hput_real( info->xmax, "xmax", dir );
  hput_real( info->ymin, "ymin", dir );
  hput_real( info->ymax, "ymax", dir );

/* r-stretching */
  if ((r_stretch_dir = create_dir("r-stretching", "r-stretching", dir)) != -1){
/* save the generic and specific stretch data */
    write_stretch( r_stretch_dir, info->r_stretch );
/* release the sub-directory */
    Vdetach(r_stretch_dir);
  }
  else
    printf("Error: save_curves_mappings: unable to create a directory for "
	   "`r-stretching'.\n");

/* s-stretching */
  if ((s_stretch_dir = create_dir("s-stretching", "s-stretching", dir)) != -1){
/* save the generic and specific stretch data */
    write_stretch( s_stretch_dir, info->s_stretch );
/* release the sub-directory */
    Vdetach(s_stretch_dir);
  }
  else
    printf("Error: save_curves_mappings: unable to create a directory for "
	   "`s-stretching'.\n");
}

static void *
cleanup_cartesian_mapping( void *data_ptr ){
  cartesian_grid_info *info;

/* cast the data to the right type */
  info = (cartesian_grid_info *) data_ptr;

  info->r_stretch = delete_generic_stretching( info->r_stretch );

  info->s_stretch = delete_generic_stretching( info->s_stretch );

  if (info != NULL)
    free(info);

  return NULL;
}

cartesian_grid_info *init_cartesian_mapping(generic_mapping *grid_ptr,
			     real xmin, real xmax, 
			     real ymin, real ymax){
  cartesian_grid_info *info;
  
/* this must be deallocated when the mapping is deleted */
  if ((info = (cartesian_grid_info *) 
    malloc (sizeof (cartesian_grid_info))) == NULL)
    printf("memory error in init_cartesian_mapping\n");

  info->xmin= xmin;
  info->xmax= xmax;
  info->ymin= ymin;
  info->ymax= ymax;

/* no stretching */
  info->r_stretch = NULL;
  info->s_stretch = NULL;

/* set the grid type and the periodicity flag */
  grid_ptr->grid_type= 1; /* Cartesian grid */
  grid_ptr->r_period = 0;
  grid_ptr->s_period = 0;

/* set the number of grid lines */
  grid_ptr->r_points = 10;
  grid_ptr->s_points = 10;

/* save the pointers to the data and the functions for
   evaluating the mapping */
  grid_ptr->grid_mapping = cartesian_mapping;
  grid_ptr->inverse_known = 1;
  grid_ptr->inverse_grid_mapping = inverse_cartesian_mapping;
  grid_ptr->mapping_data_ptr = info;
  grid_ptr->cleanup_mapping_data = cleanup_cartesian_mapping;
  grid_ptr->set_mapping = set_cartesian_mapping;
  grid_ptr->write_specific_mapping = write_cartesian_mapping;
  grid_ptr->copy_mapping_data = copy_cartesian_mapping;

  compute_cartesian_mapping( info, grid_ptr );

  return info;
}

int 
read_cartesian_mapping( int32 dir, generic_mapping *grid_ptr, generic_curve *first_c ){
  cartesian_grid_info *info;
  int32 r_stretch_dir, s_stretch_dir;

  info = (cartesian_grid_info *) malloc( sizeof(cartesian_grid_info) );

  hget_real( &(info->xmin), "xmin", dir );
  hget_real( &(info->xmax), "xmax", dir );
  hget_real( &(info->ymin), "ymin", dir );
  hget_real( &(info->ymax), "ymax", dir );

/* r-stretching */
  if ((r_stretch_dir = locate_dir("r-stretching", dir)) != -1){
/* read the generic and specific stretch data */
    info->r_stretch = read_stretching( r_stretch_dir );
/* release the sub-directory */
    Vdetach(r_stretch_dir);
  }

  if ((s_stretch_dir = locate_dir("s-stretching", dir)) != -1){
/* read the generic and specific stretch data */
    info->s_stretch = read_stretching( s_stretch_dir );
/* release the sub-directory */
    Vdetach(s_stretch_dir);
  }

  grid_ptr->mapping_data_ptr = (void *) info;
  grid_ptr->grid_mapping = cartesian_mapping;
  grid_ptr->inverse_known = 1;
  grid_ptr->inverse_grid_mapping = inverse_cartesian_mapping;
  grid_ptr->cleanup_mapping_data = cleanup_cartesian_mapping;
  grid_ptr->set_mapping = set_cartesian_mapping;
  grid_ptr->write_specific_mapping = write_cartesian_mapping;
  grid_ptr->copy_mapping_data = copy_cartesian_mapping;

  return OK;
}


static void *
copy_cartesian_mapping( void *mapping_data_ptr ){
  cartesian_grid_info *info, *old_info;
  
  info = (cartesian_grid_info *) malloc( sizeof(cartesian_grid_info) );

  old_info = (cartesian_grid_info *)mapping_data_ptr;

  info->xmin = old_info->xmin;
  info->xmax = old_info->xmax;
  info->ymin = old_info->ymin;
  info->ymax = old_info->ymax;

/* make a copy of the stretching functions */
  info->r_stretch = copy_stretching( old_info->r_stretch );
  info->s_stretch = copy_stretching( old_info->s_stretch );

  return (void *)info;
}
