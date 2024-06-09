#include "mappings.h"

/* private memeber functions */
static void set_disc_point_map(input_output *io_ptr, generic_mapping *grid_ptr, 
			       generic_curve *first_c);
static int disc_point_mapping( grid_point *gp_ptr, void *info);
static void *cleanup_disc_point_map( void *data_ptr );
static void 
write_disc_point_map( int32 dir, void *disc_point_data_ptr );
static void alloc_disc_point_arrays( disc_point_info *info );
static void free_disc_point_arrays( disc_point_info *info );
static disc_point_info *read_plot3d_file( generic_mapping *grid_ptr, 
					 FILE *fp );
static void *copy_disc_point_map( void *mapping_data_ptr );
/* end private member functions */

static void alloc_disc_point_arrays( disc_point_info *info ){
  info->x_ptr = create_real_array_2d(info->nr,info->ns);
  info->y_ptr = create_real_array_2d(info->nr,info->ns);
}

static void free_disc_point_arrays( disc_point_info *info ){
  info->x_ptr = delete_real_array_2d(info->x_ptr);
  info->y_ptr = delete_real_array_2d(info->y_ptr);
}

static disc_point_info *read_plot3d_file( generic_mapping *grid_ptr, 
					 FILE *fp ){
  int nii, njj, nkk, nr, ns, i, j, left_handed;
  disc_point_info *info;
/* formatted plot3d file with one grid in the file */
  fscanf( fp, "%i%i%i", &nii, &njj, &nkk );
/* check if it is 3d */
  if( nii > 1 && njj > 1 && nkk > 1 ){
    printf("The input file contains a truely 3-d grid, which xcog can't handle "
	   "right now.\n");
    fclose( fp );
    return NULL;
  }
  else if( nkk == 1 ){          
    nr = nii;
    ns = njj;
    left_handed = 0;
  } 
  else if( njj == 1 ){
    nr = nii;
    ns = nkk;
    left_handed = 1;
  }
  else if( nkk == 1 ){
    nr = njj;
    ns = nkk;
    left_handed = 0;
  }
  else{
    printf("The dimensions in the file make no sense. nii = %i, njj = %i, "
	   "nkk = %i.\n", nii, njj, nkk);
    fclose( fp );
    return NULL;
  }

/* create a disc_point structure */
  if ((info = (disc_point_info *) 
    malloc (sizeof (disc_point_info))) == NULL)
    printf("memory error in read_plot3d_file\n");
  
  info->nr = nr;
  info->ns = ns;
/* lets dimension the arrays */
  alloc_disc_point_arrays( info );

/* read the x coordinates */
  for (j = 1; j <= info->ns; j++){
    if (left_handed){
      for (i = info->nr; i >= 1; i--)
#ifdef SINGLE
	if ((fscanf( fp, "%g", &(x_disc(i,j)) )) == EOF){
#else
	if ((fscanf( fp, "%lg", &(x_disc(i,j)) )) == EOF){
#endif
	  printf("Input error in read_plot3d while reading x-coordinates at "
		 "i = %i, j= %i.\n", info->nr-i+1, j);
	  free_disc_point_arrays( info );
	  free( info );
	  fclose( fp );
	  return NULL;
	}
    }
    else{
      for (i = 1; i <= info->nr; i++)
#ifdef SINGLE
	if ((fscanf( fp, "%g", &(x_disc(i,j)) )) == EOF){
#else
	if ((fscanf( fp, "%lg", &(x_disc(i,j)) )) == EOF){
#endif
	  printf("Input error in read_plot3d while reading x-coordinates at "
		 "i = %i, j= %i.\n", i, j);
	  free_disc_point_arrays( info );
	  free( info );
	  fclose( fp );
	  return NULL;
	}
    }  
  }
/* read the y coordinates */
  for (j = 1; j <= info->ns; j++){
    if (left_handed){
      for (i = info->nr; i >= 1; i--)
#ifdef SINGLE
	if ((fscanf( fp, "%g", &(y_disc(i,j)) )) == EOF){
#else
	if ((fscanf( fp, "%lg", &(y_disc(i,j)) )) == EOF){
#endif
	  printf("Input error in read_plot3d while reading y-coordinates at "
		 "i = %i, j= %i.\n", info->nr-i+1, j);
	  free_disc_point_arrays( info );
	  free( info );
	  fclose( fp );
	  return NULL;
	}
    }
    else{
      for (i = 1; i <= info->nr; i++)
#ifdef SINGLE
	if ((fscanf( fp, "%g", &(y_disc(i,j)) )) == EOF){
#else
	if ((fscanf( fp, "%lg", &(y_disc(i,j)) )) == EOF){
#endif
	  printf("Input error in read_plot3d while reading y-coordinates at "
		 "i = %i, j= %i.\n", i, j);
	  free_disc_point_arrays( info );
	  free( info );
	  fclose( fp );
	  return NULL;
	}
    }
  }
/* successful completion */
  fclose( fp );

  return info;
}


static void compute_disc_point_map( disc_point_info *info, 
				   generic_mapping *grid_ptr ){
  int i,j;

/* compute bounding box */
  grid_ptr->x_min_0 = 1.e10;
  grid_ptr->x_max_0 = -1.e10;
  grid_ptr->y_min_0 = 1.e10;
  grid_ptr->y_max_0 = -1.e10;
/* i=1 */
  i = 1;
  for (j=1; j <= info->ns; j++){
    grid_ptr->x_min_0 = real_min( x_disc(i,j), grid_ptr->x_min_0 );
    grid_ptr->x_max_0 = real_max( x_disc(i,j), grid_ptr->x_max_0 );
    grid_ptr->y_min_0 = real_min( y_disc(i,j), grid_ptr->y_min_0 );
    grid_ptr->y_max_0 = real_max( y_disc(i,j), grid_ptr->y_max_0 );
  }
/* i=info->nr */
  i = info->nr;
  for (j=1; j <= info->ns; j++){
    grid_ptr->x_min_0 = real_min( x_disc(i,j), grid_ptr->x_min_0 );
    grid_ptr->x_max_0 = real_max( x_disc(i,j), grid_ptr->x_max_0 );
    grid_ptr->y_min_0 = real_min( y_disc(i,j), grid_ptr->y_min_0 );
    grid_ptr->y_max_0 = real_max( y_disc(i,j), grid_ptr->y_max_0 );
  }
/* j=1 */
  j = 1;
  for (i=1; i <= info->nr; i++){
    grid_ptr->x_min_0 = real_min( x_disc(i,j), grid_ptr->x_min_0 );
    grid_ptr->x_max_0 = real_max( x_disc(i,j), grid_ptr->x_max_0 );
    grid_ptr->y_min_0 = real_min( y_disc(i,j), grid_ptr->y_min_0 );
    grid_ptr->y_max_0 = real_max( y_disc(i,j), grid_ptr->y_max_0 );
  }
/* j=info->ns */
  j = info->ns;
  for (i=1; i <= info->nr; i++){
    grid_ptr->x_min_0 = real_min( x_disc(i,j), grid_ptr->x_min_0 );
    grid_ptr->x_max_0 = real_max( x_disc(i,j), grid_ptr->x_max_0 );
    grid_ptr->y_min_0 = real_min( y_disc(i,j), grid_ptr->y_min_0 );
    grid_ptr->y_max_0 = real_max( y_disc(i,j), grid_ptr->y_max_0 );
  }

  reset_view_port( grid_ptr );
}


static void set_disc_point_map(input_output *io_ptr, generic_mapping *grid_ptr, 
			       generic_curve *first_c){
  char prompt[80];
  int icom, quit, i, replot, ip, jp;
  disc_point_info *info;
  grid_point gp;
  int level=2, no_indent=0;

#include "disc_point_map_com.h"

  sprintf(prompt, "%s: discrete point mapping>", grid_ptr->grid_name);

/* cast the data pointer to the right type */
  info = (disc_point_info *) grid_ptr->mapping_data_ptr;

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

    case R_STRETCH:
/* r-stretching */
      if (info->r_stretch == NULL){
	info->r_stretch = choose_stretching( io_ptr, NULL, &grid_ptr->r_points );
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

    case NO_R_STRETCH:
/* no-r-stretching */
      if (info->r_stretch != NULL){
	info->r_stretch = delete_generic_stretching( info->r_stretch );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	printf("There is no stretching in the s-direction.\n");
	replot = 0;
      }
      break;

    case S_STRETCH:
/* s-stretching */
      if (info->s_stretch == NULL){
	info->s_stretch = choose_stretching( io_ptr, NULL, &grid_ptr->s_points );
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

    case NO_S_STRETCH:
/* no-s-stretching */
      if (info->s_stretch != NULL){
	info->s_stretch = delete_generic_stretching( info->s_stretch );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	printf("There is no stretching in the s-direction.\n");
	replot = 0;
      }
      break;

    case R_LINES:
/* r-lines */
      grid_ptr->r_points = 
	int_max(2, get_int(io_ptr, 
			   "Enter number of gridlines in r >= 2: ", 
			   grid_ptr->r_points,
			   no_indent));

      grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      break;

    case S_LINES:
/* s-lines */
      grid_ptr->s_points = 
	int_max(2, get_int(io_ptr, 
			   "Enter number of gridlines in s>= 2: ", grid_ptr->s_points,
			   no_indent));

      grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      break;

/* compare point values with interpolated values */
/*     case TEST: */
/*       ip = int_max(1, get_int(io_ptr, "Enter i-point: ", 1, no_indent)); */
/*       ip = int_min(ip, info->nr); */
/*       jp = int_max(1, get_int(io_ptr, "Enter j-point: ", 1, no_indent)); */
/*       jp = int_min(jp, info->ns); */
/*       gp.r = real_max( 0.0, get_real(io_ptr, "Enter r: ", 0.5, no_indent)); */
/*       gp.r = real_min( 1.0, gp.r ); */
/*       gp.s = real_max( 0.0, get_real(io_ptr, "Enter s: ", 0.5, no_indent)); */
/*       gp.s = real_min( 1.0, gp.s ); */
/*       disc_point_mapping(&gp, info); */
/*       printf("Point coordinate: (%e, %e)\n", x_disc(ip,jp), y_disc(ip,jp)); */
/*       printf("Mapping coordinate: (%e, %e)\n", gp.x, gp.y); */
/*       break; */

    case ORDER:
      do{
	info->order = get_int(io_ptr, 
			      "Enter interpolation order (2 for linear, 4 for cubic): ", 
			      info->order, no_indent);
	if (!(info->order == 2 || info->order == 4))
	  printf("The interpolation order must be 2 or 4!!!\n");
      }	while (!(info->order == 2 || info->order == 4));
      printf("Using %s interpolation.\n", (info->order == 2)? "LINEAR" : "CUBIC");

      break;

    case PLOT_MODE:
/* plot mode */
      set_grid_plot_mode( io_ptr, grid_ptr );
      replot = 0;
      break;

    case SHOW:
/* show */
      printf("Using %s interpolation.\n\n", (info->order == 2)? "LINEAR" : "CUBIC");
      printf("The mapping was defined with %i grid point in the r-direction\n"
	     "and %i grid points in the s-direction.\n\n", info->nr, info->ns);
      show_mapping_parameters( grid_ptr );
      replot = 0;
      break;

    case HELP:
/* help */
      while ( (icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				  command, ncom, level+1, NULL, NULL)) == -1 );
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case EXIT:
/* exit */
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);
#undef ncom
}


static int 
disc_point_mapping( grid_point *gp_ptr, void *data_ptr){
  disc_point_info *info;
  int i,j,ip,jp;
  int i_loc, j_loc, i_close, j_close, j_m, j_p, i_m, i_p;
  const int iw_low=1, iw_high=2, iw1=3;
  real xi, eta, d_xi, d_eta, r_stretch, r_ratio, s_stretch, s_ratio, dummy, dr, ds;

#include "cubic_interp.h"

/* cast the pointer to the right type */
  info = (disc_point_info *) data_ptr;

  dr = 1.0/((real) info->nr - 1);
  ds = 1.0/((real) info->ns - 1);

/* r-stretching */
  uniform_to_stretch( gp_ptr->r, info->r_stretch, &r_stretch, &r_ratio, &dummy );

/* s-stretching */
  uniform_to_stretch( gp_ptr->s, info->s_stretch, &s_stretch, &s_ratio, &dummy );

  if (info->r_period){
    if (r_stretch < 0.0)
      r_stretch += 1.0;
    else if (r_stretch > 1.0)
      r_stretch -= 1.0;
  }

  if (info->s_period){
    if (s_stretch < 0.0)
      s_stretch += 1.0;
    else if (s_stretch > 1.0)
      s_stretch -= 1.0;
  }

  if (info->order == 2)
    {
/* bi-linear interpolation */
/* get the nearest grid point below and to the left of (r,s) */
      i = int_max( 1, 1 + (int) ( (info->nr-1)*r_stretch ) ); 
      j = int_max( 1, 1 + (int) ( (info->ns-1)*s_stretch ) ); 

/* prevent overflow */
      i = int_min( i, info->nr-1 ); 
      j = int_min( j, info->ns-1 ); 

/* step sizes */
      d_xi  = 1.0/ (info->nr - 1); 
      d_eta = 1.0/ (info->ns - 1); 

/* local coordinates */
      xi = (r_stretch - (i-1)*d_xi)/d_xi; 
      eta = (s_stretch - (j-1)*d_eta)/d_eta; 

/* bi-linear interpolation */
      gp_ptr->x =  
	(1.0 - xi) * (1.0 - eta) * x_disc(i,j) + 
	  xi * (1.0 - eta) * x_disc(i+1,j) + 
	    (1.0 - xi) * eta * x_disc(i,j+1) + 
	      xi * eta * x_disc(i+1,j+1); 
      gp_ptr->y =  
	(1.0 - xi) * (1.0 - eta) * y_disc(i,j) + 
	  xi * (1.0 - eta) * y_disc(i+1,j) + 
	    (1.0 - xi) * eta * y_disc(i,j+1) + 
	      xi * eta * y_disc(i+1,j+1); 

/* jacobian */
      gp_ptr->xr =  
	(- (1.0 - eta) * x_disc(i,j) + (1.0 - eta) * x_disc(i+1,j) - 
	 eta * x_disc(i,j+1) + eta * x_disc(i+1,j+1) )/d_xi * r_ratio; 
      gp_ptr->xs = 
	( - (1.0 - xi) * x_disc(i,j) - xi * x_disc(i+1,j) + 
	 (1.0 - xi) * x_disc(i,j+1) + xi * x_disc(i+1,j+1) )/d_eta * s_ratio; 
      gp_ptr->yr =  
	(- (1.0 - eta) * y_disc(i,j) + (1.0 - eta) * y_disc(i+1,j) - 
	 eta * y_disc(i,j+1) + eta * y_disc(i+1,j+1) )/d_xi * r_ratio; 
      gp_ptr->ys = 
	( - (1.0 - xi) * y_disc(i,j) - xi * y_disc(i+1,j) + 
	 (1.0 - xi) * y_disc(i,j+1) + xi * y_disc(i+1,j+1) )/d_eta * s_ratio; 

    }
  else if (info->order == 4)
    {

/* bi-cubic interpolation */
/* get the nearest grid point below and to the left of (r,s) */
      i_close = 1 + (int) ( (info->nr-1)*r_stretch ); 
      j_close = 1 + (int) ( (info->ns-1)*s_stretch ); 

      if (i_close-iw_low < 1) 
	{
	  if (info->r_period)
	    {
	      i_m = info->nr-1; 
	      i_loc = 0;
	      i_p = 3; 
	    }
	  else
	    {
	      i_m = i_loc = 1; 
	      i_p = 4;
	    }
	}
      else if (i_close+iw_high > info->nr) 
	{
	  if (info->r_period)
	    {
	      i_p = 2; 
	      i_m = i_loc = info->nr-2; 
	    }
	  else
	    {
	      i_m = i_loc = info->nr-iw1; 
	      i_p = info->nr;
	    }
	}
      else 
	{
	  i_m = i_loc = i_close-iw_low;
	  i_p = i_m + 3;
	}
      
      if (j_close-iw_low < 1) 
	{ 
	  if (info->s_period) 
	    { 
	      j_m = info->ns-1; 
	      j_loc = 0;
	      j_p = 3; 
	    } 
	  else
	    { 
	      j_m = j_loc = 1; 
	      j_p = 4; 
	    } 
	} 
      else if (j_close+iw_high > info->ns) 
	{ 
	  if (info->s_period) 
	    { 
	      j_p = 2; 
	      j_m = j_loc = info->ns-2; 
	    } 
	  else
	    { 
	      j_m = j_loc = info->ns-iw1; 
	      j_p = info->ns; 
	    } 
	} 
      else
	{ 
	  j_m = j_loc = j_close-iw_low; 
	  j_p = j_m + 3; 
	} 

      gp_ptr->x  = CUB(r_stretch, s_stretch, x_disc); 
      gp_ptr->y  = CUB(r_stretch, s_stretch, y_disc); 
      gp_ptr->xr = CUB_R(r_stretch, s_stretch, x_disc); 
      gp_ptr->yr = CUB_R(r_stretch, s_stretch, y_disc); 
      gp_ptr->xs = CUB_S(r_stretch, s_stretch, x_disc); 
      gp_ptr->ys = CUB_S(r_stretch, s_stretch, y_disc); 

    }
  else
    {
      printf("ERROR: order in disc_point_mapping must be 2 or 4 and not %i\n", 
	     info->order);
      gp_ptr->x  = 0.;
      gp_ptr->y  = 0.;
      gp_ptr->xr = 0.;
      gp_ptr->yr = 0.;
      gp_ptr->xs = 0.;
      gp_ptr->ys = 0.;
      return 1;
    }

  return 0;
}

static void *cleanup_disc_point_map( void *data_ptr ){
  disc_point_info *info;

/* cast the data pointer to the right type */
  info = (disc_point_info *) data_ptr;

/* free the memory */
  free_disc_point_arrays( info );

  delete_generic_stretching( info->r_stretch );
  delete_generic_stretching( info->s_stretch );

  free( info );

  return NULL;
}

static void 
write_disc_point_map( int32 dir, void *data_ptr ){
  disc_point_info *info;
  int32 stretch_dir;

/* cast the data pointer to the right type */
  info = (disc_point_info *) data_ptr;

/* save the sizes */
  hput_int( info->nr, "nr", dir );
  hput_int( info->ns, "ns", dir );

/* save the arrays */
  hput_real_array_2d( info->x_ptr, "x", dir );
  hput_real_array_2d( info->y_ptr, "y", dir );

/* r-stretching */
  if ((stretch_dir = create_dir("r-stretching", 
				"r-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( stretch_dir, info->r_stretch );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }
  else
    printf("Error: write_linear_interp: unable to create a directory for "
	   "`r-stretching'.\n");

/* s-stretching */
  if ((stretch_dir = create_dir("s-stretching", 
				"s-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( stretch_dir, info->s_stretch );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }
  else
    printf("Error: write_linear_interp: unable to create a directory for "
	   "`s-stretching'.\n");
}


disc_point_info *init_disc_point_map( input_output *io_ptr, 
				     generic_mapping *grid_ptr ){
  disc_point_info *info=NULL;
  FILE *fp;
  char *file_name, *prompt;
  int icom, level=1, quit, i, j;
  real dist, box_size;
  const int save_command=1;

#include "init_disc_point_com.h"

  prompt = "select file format for reading grid points>";

  do{

    icom = get_command(io_ptr, prompt, command, ncom, level,
		       save_on_copy, argument);

    quit = 0;
    switch (icom) {

    case 0:
/* plot3d ascii format */
      if ( (fp = open_ascii_file( io_ptr, "Enter plot3d ascii file: ", 
				 "test.plot3d", &file_name, 'r', 0, 
				 save_command)) == NULL ||
	  (info = read_plot3d_file( grid_ptr, fp )) == NULL )
	return NULL;
      quit = 1;
      break;

    case 1:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case 2:
/* cancel */
      return NULL;
      break;

    default:
      ;

    }
  } while( !quit );


/* no stretching */
  info->r_stretch = NULL;
  info->s_stretch = NULL;

/* set the grid type */
  grid_ptr->grid_type= 6;
  
/* set the number of grid lines */
  grid_ptr->r_points = info->nr;
  grid_ptr->s_points = info->ns;
  
/* save the pointers to the data and the functions for
   evaluating the mapping */
  grid_ptr->grid_mapping = disc_point_mapping;
  grid_ptr->inverse_known = 0;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->mapping_data_ptr = (void *) info;
  grid_ptr->cleanup_mapping_data = cleanup_disc_point_map;
  grid_ptr->set_mapping = set_disc_point_map;
  grid_ptr->write_specific_mapping = write_disc_point_map;
  grid_ptr->copy_mapping_data = copy_disc_point_map;
  
/* cubic interpolation by default */
  info->order = 4;

  compute_disc_point_map( info, grid_ptr );

/* get the size of the bounding box */
  box_size = real_max(grid_ptr->x_max_0 - grid_ptr->x_min_0, 
		      grid_ptr->y_max_0 - grid_ptr->y_min_0);

/* check for periodicity in the r-direction */
  dist = 0;
  for (j=1; j<=grid_ptr->s_points; j++)
    dist += fabs(x_disc(1,j) - x_disc(grid_ptr->r_points,j)) +
	fabs(y_disc(1,j) - y_disc(grid_ptr->r_points,j));

  info->r_period = grid_ptr->r_period = (dist <= 1.e-4 * box_size)? 1 : 0;

/* check for periodicity in the s-direction */
  dist = 0;
  for (i=1; i<=grid_ptr->r_points; i++)
    dist += fabs(x_disc(i,1) - x_disc(i,grid_ptr->s_points)) +
	fabs(y_disc(i,1) - y_disc(i,grid_ptr->s_points));

  info->s_period = grid_ptr->s_period = (dist <= 1.e-4 * box_size)? 1 : 0;

/* tmp */
  if (grid_ptr->r_period) printf("The grid seems to be periodic in r\n");
  if (grid_ptr->s_period) printf("The grid seems to be periodic in s\n");
/* end tmp */

  return info;
}


int
read_disc_point_map( int32 dir, generic_mapping *grid_ptr, generic_curve *first_c ){
  disc_point_info *info;
  int32 stretch_dir;

/* cast the data pointer to the right type */
  info = (disc_point_info *) malloc( sizeof( disc_point_info ) );

/* retrieve the sizes */
  hget_int( &(info->nr), "nr", dir );
  hget_int( &(info->ns), "ns", dir );

/* retrieve the arrays */
  info->x_ptr = hget_real_array_2d( "x", dir );
  info->y_ptr = hget_real_array_2d( "y", dir );

/* r-stretching */
  if ((stretch_dir = locate_dir("r-stretching", dir)) != -1){
/* read the generic and specific curve data */
    info->r_stretch = read_stretching( stretch_dir );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }

/* s-stretching */
  if ((stretch_dir = locate_dir("s-stretching", dir)) != -1){
/* read the generic and specific curve data */
    info->s_stretch = read_stretching( stretch_dir );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }

  grid_ptr->mapping_data_ptr = (void *) info;
  grid_ptr->grid_mapping = disc_point_mapping;
  grid_ptr->inverse_known = 0;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->cleanup_mapping_data = cleanup_disc_point_map;
  grid_ptr->set_mapping = set_disc_point_map;
  grid_ptr->write_specific_mapping = write_disc_point_map;
  grid_ptr->copy_mapping_data = copy_disc_point_map;

/* copy the periodicity flag */
  info->r_period = grid_ptr->r_period;
  info->s_period = grid_ptr->s_period;

  return OK;
}

static void *
copy_disc_point_map( void *mapping_data_ptr ){
  disc_point_info *info, *old_info;
  int i,j;
  
#define old_x_disc(i,j) compute_index_2d(old_info->x_ptr,i,j)
#define old_y_disc(i,j) compute_index_2d(old_info->y_ptr,i,j)

  info = (disc_point_info *) malloc( sizeof(disc_point_info) );

  old_info = (disc_point_info *)mapping_data_ptr;

  info->nr = old_info->nr;
  info->ns = old_info->ns;

/* make a (deep) copy of the stretching functions */
  info->r_stretch = copy_stretching( old_info->r_stretch );
  info->s_stretch = copy_stretching( old_info->s_stretch );

/* allocate space for the arrays */
  alloc_disc_point_arrays( info );

/* copy the array elements */
  for (i=1; i<=info->nr; i++)
    for (j=1; j<=info->ns; j++){
      x_disc(i,j) = old_x_disc(i,j);
      y_disc(i,j) = old_y_disc(i,j);
    }

  return (void *)info;
#undef old_x_disc
#undef old_y_disc
}

