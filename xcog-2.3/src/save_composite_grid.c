#include "xcog.h"

#include "hdf_stuff.h"

/* static prototypes */
static void 
ascii_save_xcog_data( FILE *fp, overlapping_grid *xcog_data_ptr, 
		     int save_jacobian );
static void 
hdf_save_xcog_data(int32 root, overlapping_grid *xcog_data_ptr, 
		   int save_jacobian );
static void 
hdf_save_component_grid( int32 comp_dir, component_grid *grid_ptr, 
			int save_jacobian );
static void 
compute_jacobian( component_grid *grid_ptr, 
		 real_array_2d *xr_ptr, real_array_2d *xs_ptr, 
		 real_array_2d *yr_ptr, real_array_2d *ys_ptr );
/* end static prototypes */


void 
save_xcog_data( input_output *io_ptr, overlapping_grid *xcog_data_ptr ){
  FILE *fp;
  int32 root;

  char prompt[80], *file_name, *name;
  int icom, quit, save_jacobian;
  int level=1;
  const int save_command=1;

#include "save_composite_grid_com.h"

  sprintf(prompt, "save-overlapping-grid>");

  do{

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    quit = 0;
    switch (icom) {

    case 0:
/* ascii file */
      if ((fp = open_ascii_file( io_ptr, "File name for composite grid: ",
				"test.acg", &file_name, 'w', 5, 
				save_command) ) != NULL ){
	save_jacobian = get_yes_no( io_ptr, "Do you want to save the Jacobian "
				   "of the transformation in every grid point>"
				   , level );
	ascii_save_xcog_data( fp, xcog_data_ptr, save_jacobian );
	fclose( fp );
	quit = 1;
      }
      break;

    case 1:
/* hdf file */
      name = get_word( io_ptr, "File name for overlapping grid: ", "test.hdf", 1);
  
      if ((root = open_hdf_file(name, 'i')) > 0){
	save_jacobian = get_yes_no( io_ptr, "Do you want to save the Jacobian "
				   "of the transformation in every grid point>"
				   , level );
/* save the xcog data */
	hdf_save_xcog_data( root, xcog_data_ptr, save_jacobian );
/* close the database file */
	close_hdf_file(root);
	quit = 1;
      }
      else{
	printf("Unable to open the database file %s\n", name);
      }

      break;

    case 2:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>",
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case 3:
/* cancel */
      quit = 1;
      break;
      
    default:
      ;
    }
  } while (!quit);

}

static void 
ascii_save_xcog_data( FILE *fp, overlapping_grid *xcog_data_ptr, int save_jacobian ){
  component_grid *grid_ptr;
  interp_point *interp_ptr;
  char *time_string;
  time_t present_time;
  real_array_2d *xr_ptr, *xs_ptr, *yr_ptr, *ys_ptr;

/* start composite grid */
  fprintf( fp, "@begin_xcog_data\n" );

/* get the present time */
  present_time = time( NULL );
  time_string = ctime( &present_time );

/* save an info string */
  fprintf( fp, "%s%s", "@info_string Composite grid created by Xcog version " 
	  VERSION " on ", time_string );

  fprintf( fp, "#\n# Composite grid parameters:\n#\n" );  
/* save the composite grid stuff */
  fprintf( fp, "@n_grids %i\n", xcog_data_ptr->n_grids );
  fprintf( fp, "@extra %i\n", xcog_data_ptr->extra );
  fprintf( fp, "@extra_period %i\n", xcog_data_ptr->extra_period );
  fprintf( fp, "@disc_width %i\n", xcog_data_ptr->disc_width );
  fprintf( fp, "@normal_width %i\n", xcog_data_ptr->normal_width );
  fprintf( fp, "@tangent_width %i\n", xcog_data_ptr->tangent_width );
  fprintf( fp, "@corner_width %i\n", xcog_data_ptr->corner_width );
  fprintf( fp, "@interp_width %i\n", xcog_data_ptr->interp_width );
  fprintf( fp, "@interp_type %c\n", xcog_data_ptr->interp_type );
  fprintf( fp, "@x_min %.18G\n", xcog_data_ptr->x_min );
  fprintf( fp, "@x_max %.18G\n", xcog_data_ptr->x_max );
  fprintf( fp, "@y_min %.18G\n", xcog_data_ptr->y_min );
  fprintf( fp, "@y_max %.18G\n", xcog_data_ptr->y_max );

  fprintf( fp, "#\n# Size information for the component grids:\n#\n" );  
  fprintf( fp, "@begin_size_information\n" );  
  for (grid_ptr = xcog_data_ptr->first; grid_ptr != NULL;
       grid_ptr = grid_ptr->next){
    fprintf( fp, "#\n# A new component grid starts on next line.\n#\n" );
    fprintf( fp, "@begin_component_grid\n" );
    fprintf( fp, "@grid_name %s\n", grid_ptr->grid_name );
    fprintf( fp, "@r_dim %i\n", grid_ptr->r_dim );
    fprintf( fp, "@s_dim %i\n", grid_ptr->s_dim );
    fprintf( fp, "@r_period %i\n", grid_ptr->r_period );
    fprintf( fp, "@s_period %i\n", grid_ptr->s_period );
    fprintf( fp, "# Number of interpolation points:\n" );
    fprintf( fp, "@n_interp %i\n", grid_ptr->n_interp );
    fprintf( fp, 
"# The flag save_jacobian == 1 if the Jacobian of the mapping is saved for\n"
"# all grid points.\n" );
    fprintf( fp, "@save_jacobian %i\n", save_jacobian );
    fprintf( fp, "@number_of_real_arrays %i\n", 2 + save_jacobian*4 );
    fprintf( fp, "@number_of_int_arrays %i\n", 1 );

    fprintf( fp, "@end_component_grid\n" );
  }
  fprintf( fp, "@end_size_information\n" );  

  fprintf( fp, "#\n# Complete information for the component grids:\n#\n" );  
  fprintf( fp, "@begin_component_grid_information\n" );  
/* save the mapping stuff */
  for (grid_ptr = xcog_data_ptr->first; grid_ptr != NULL;  
       grid_ptr = grid_ptr->next){ 
/* first save all scalar quantities */
    fprintf( fp, "#\n# A new component grid starts on next line.\n#\n" );
    fprintf( fp, "@begin_component_grid\n" );
    fprintf( fp, "@grid_name %s\n", grid_ptr->grid_name );
    fprintf( fp, "@priority %i\n", grid_ptr->priority );
    fprintf( fp, "@grid_type %i\n", grid_ptr->grid_type );
    fprintf( fp, "@r_points %i\n", grid_ptr->r_points );
    fprintf( fp, "@s_points %i\n", grid_ptr->s_points );
    fprintf( fp, "@r_dim %i\n", grid_ptr->r_dim );
    fprintf( fp, "@s_dim %i\n", grid_ptr->s_dim );
    fprintf( fp, "@r_period %i\n", grid_ptr->r_period );
    fprintf( fp, "@s_period %i\n", grid_ptr->s_period );
    fprintf( fp, "# Number of interpolation points:\n" );
    fprintf( fp, "@n_interp %i\n", grid_ptr->n_interp );
    fprintf( fp, "@r_step %.18G\n", grid_ptr->r_step );
    fprintf( fp, "@s_step %.18G\n", grid_ptr->s_step );
    fprintf( fp, "@x_min %.18G\n", grid_ptr->x_min );
    fprintf( fp, "@x_max %.18G\n", grid_ptr->x_max );
    fprintf( fp, "@y_min %.18G\n", grid_ptr->y_min );
    fprintf( fp, "@y_max %.18G\n", grid_ptr->y_max );

/* variables to distinguish C- and H-topologies from normal grids */
    fprintf( fp, "# Topology information:\n" );
    fprintf( fp, "@topology %i\n", grid_ptr->topology );
    fprintf( fp, "@branch_direction %i\n", grid_ptr->branch_direction );
    fprintf( fp, "@first_point %i\n", grid_ptr->first_point );
    fprintf( fp, "@last_point %i\n",  grid_ptr->last_point );

/* save all array quantities */
    fprintf( fp,
"# range(is,id): array range on side is=1,2 in direction id=1,2 (r=1, s=2)\n" );
    ascii_int_array_2d( grid_ptr->range_ptr, fp, "range" );
    fprintf( fp,
"# bc(is,id): boundary condition on side is=1,2 in direction id=1,2 (r=1, s=2)\n" );
    ascii_int_array_2d( grid_ptr->bc_ptr, fp, "bc" );
    fprintf( fp,
"# curve(is,id): curve value on side is=1,2 in direction id=1,2 (r=1, s=2)\n" );
    ascii_int_array_2d( grid_ptr->curve_ptr, fp, "curve" );
    ascii_int_array_2d( grid_ptr->flag_ptr, fp, "flag" );
    ascii_real_array_2d( grid_ptr->x_ptr, fp, "x" );
    ascii_real_array_2d( grid_ptr->y_ptr, fp, "y" );

/* save info if the jacobian is saved or not */
    fprintf( fp, "@save_jacobian %i\n", save_jacobian );

/* only save the jacobian if save_jacobian == 1 */
    if (save_jacobian == 1){
      xr_ptr = create_real_array_2d(grid_ptr->r_dim, grid_ptr->s_dim);
      xs_ptr = create_real_array_2d(grid_ptr->r_dim, grid_ptr->s_dim);
      yr_ptr = create_real_array_2d(grid_ptr->r_dim, grid_ptr->s_dim);
      ys_ptr = create_real_array_2d(grid_ptr->r_dim, grid_ptr->s_dim);

      compute_jacobian( grid_ptr, xr_ptr, xs_ptr, yr_ptr, ys_ptr );

      ascii_real_array_2d( xr_ptr, fp, "xr" ); 
      ascii_real_array_2d( xs_ptr, fp, "xs" ); 
      ascii_real_array_2d( yr_ptr, fp, "yr" ); 
      ascii_real_array_2d( ys_ptr, fp, "ys" ); 

      xr_ptr = delete_real_array_2d( xr_ptr );
      xs_ptr = delete_real_array_2d( xs_ptr );
      yr_ptr = delete_real_array_2d( yr_ptr );
      ys_ptr = delete_real_array_2d( ys_ptr );
    }
/* save the interpolation points in array form */
    fprintf( fp, "#\n# Interpolation information.\n#\n" );
/* i_point */
    fprintf( fp, "@int_array_1d %s\n", "i_point" );
/* output the size */
    fprintf( fp, "&dim1 %i\n", grid_ptr->n_interp );
/* output the data */
    fprintf( fp, "# The first index changes the fastest.\n" );
    fprintf( fp, "&int_array_1d_data\n" );
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev)
      fprintf( fp, "%i\n", interp_ptr->i_point );

/* j_point */
    fprintf( fp, "@int_array_1d %s\n", "j_point" );
/* output the size */
    fprintf( fp, "&dim1 %i\n", grid_ptr->n_interp );
/* output the data */
    fprintf( fp, "# The first index changes the fastest.\n" );
    fprintf( fp, "&int_array_1d_data\n" );
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev)
      fprintf( fp, "%i\n", interp_ptr->j_point );

/* grid_loc->priority */
    fprintf( fp, "@int_array_1d %s\n", "grid_loc->priority" );
/* output the size */
    fprintf( fp, "&dim1 %i\n", grid_ptr->n_interp );
/* output the data */
    fprintf( fp, "# The first index changes the fastest.\n" );
    fprintf( fp, "&int_array_1d_data\n" );
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev)
      fprintf( fp, "%i\n", interp_ptr->grid_loc->priority );

/* i_loc */
    fprintf( fp, "@int_array_1d %s\n", "i_loc" );
/* output the size */
    fprintf( fp, "&dim1 %i\n", grid_ptr->n_interp );
/* output the data */
    fprintf( fp, "# The first index changes the fastest.\n" );
    fprintf( fp, "&int_array_1d_data\n" );
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev)
      fprintf( fp, "%i\n", interp_ptr->i_loc );

/* j_loc */
    fprintf( fp, "@int_array_1d %s\n", "j_loc" );
/* output the size */
    fprintf( fp, "&dim1 %i\n", grid_ptr->n_interp );
/* output the data */
    fprintf( fp, "# The first index changes the fastest.\n" );
    fprintf( fp, "&int_array_1d_data\n" );
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev)
      fprintf( fp, "%i\n", interp_ptr->j_loc );

/* r_loc */
    fprintf( fp, "@real_array_1d %s\n", "r_loc" );
/* output the size */
    fprintf( fp, "&dim1 %i\n", grid_ptr->n_interp );
/* output the data */
    fprintf( fp, "# The first index changes the fastest.\n" );
    fprintf( fp, "&real_array_1d_data\n" );
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev)
      fprintf( fp, "%.18G\n", interp_ptr->r_loc );

/* s_loc */
    fprintf( fp, "@real_array_1d %s\n", "s_loc" );
/* output the size */
    fprintf( fp, "&dim1 %i\n", grid_ptr->n_interp );
/* output the data */
    fprintf( fp, "# The first index changes the fastest.\n" );
    fprintf( fp, "&real_array_1d_data\n" );
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev)
      fprintf( fp, "%.18G\n", interp_ptr->s_loc );

    fprintf( fp, "@end_component_grid\n" );
  }

/* end component grid information */
  fprintf( fp, "@end_component_grid_information\n" );  

/* end composite grid */
  fprintf( fp, "@end_xcog_data\n" );
}


static void 
hdf_save_xcog_data(int32 root, overlapping_grid *xcog_data_ptr, 
		   int save_jacobian ){
  int32 over_dir, comp_dir;
  char comp_dir_name[80], creator[120];
  component_grid *grid_ptr;
  real_array_2d *bbox_;
  time_t today;

#define bbox(i,j) compute_index_2d(bbox_, i, j)

  today = time(NULL);
  sprintf(creator, "xcog version "VERSION" generated this file on %s", ctime(&today));
  hput_string(creator, "creator", root);

/* make a directory for the overlapping grid */
  if ((over_dir = create_dir("overlapping grid", "overlapping grid", root)) != -1){

/* save all data for the overlapping grid */
    hput_string("The one and only overlapping grid in this file!", 
		"overlapping grid name", over_dir);

#ifdef SINGLE
    hput_int(32, "precision", over_dir);
#else
    hput_int(64, "precision", over_dir);
#endif

/* this is a 2-D grid */
    hput_int(2, "number of dimensions", over_dir);

    hput_int(xcog_data_ptr->n_grids, "n_components", over_dir);
    hput_string((xcog_data_ptr->interp_type == 'e')? "explicit" : "implicit", 
		"interpolation type", over_dir);
    hput_int(xcog_data_ptr->interp_width, "interpolation width", over_dir);
    hput_int(xcog_data_ptr->disc_width, "discretization width", over_dir);
    hput_int(xcog_data_ptr->normal_width, "normal width", over_dir);
    hput_int(xcog_data_ptr->tangent_width, "tangent width", over_dir);
    hput_int(xcog_data_ptr->corner_width, "corner width", over_dir);
    hput_int(xcog_data_ptr->extra, "ghost points", over_dir);
    hput_int(xcog_data_ptr->extra_period, "periodic overlap", over_dir);
/*    hput_int(xcog_data_ptr->trim_style, "trim style", over_dir); */
/*    hput_real(xcog_data_ptr->max_n_dist, "max curve missfit", over_dir); */

/* make a real array 2d to put the bounding box in */
    bbox_ = create_real_array_2d(2,2);
    bbox(1,1) = xcog_data_ptr->x_min;
    bbox(2,1) = xcog_data_ptr->x_max;
    bbox(1,2) = xcog_data_ptr->y_min;
    bbox(2,2) = xcog_data_ptr->y_max;
/* save the bounding box */
    hput_real_array_2d(bbox_, "bounding box", over_dir);
    bbox_ = delete_real_array_2d( bbox_ );

/* save all component grids */
    for (grid_ptr = xcog_data_ptr->first; grid_ptr != NULL; 
	 grid_ptr = grid_ptr->next){
    
      sprintf(comp_dir_name, "component grid %i", grid_ptr->priority);
      if ((comp_dir = create_dir(comp_dir_name, "component grid directory", 
				 over_dir)) != -1){
	hdf_save_component_grid(comp_dir, grid_ptr, save_jacobian);
/* release the sub-directory */
	Vdetach(comp_dir);
      }
      else
	printf("Error: save_overlapping_grid: unable to create a directory for %s\n",
	       comp_dir_name);
    }
/* release the overlapping grid-directory */
    Vdetach(over_dir);
  }
  else
    printf("Error: save_overlapping_grid: unable to create a directory for "
	   "the overlapping grid\n");
}

static void 
hdf_save_component_grid(int32 comp_dir, component_grid *grid_ptr, 
			int save_jacobian ){
  int i;
  int_array_1d *dimension_, *periodicity_;
  real_array_1d *step_;
  real_array_2d *xr_ptr, *xs_ptr, *yr_ptr, *ys_ptr, *donor_par_, *bbox_;
  int_array_2d *inter_point_, *donor_point_;
  interp_point *interp_ptr;
#define donor_par(i,j)   compute_index_2d(donor_par_, i, j)
#define donor_point(i,j) compute_index_2d(donor_point_, i, j)
#define inter_point(i,j) compute_index_2d(inter_point_, i, j)
#define bbox(i,j) compute_index_2d(bbox_, i, j)
#define dimension(i) compute_index_1d(dimension_, i)
#define periodicity(i) compute_index_1d(periodicity_, i)
#define step(i) compute_index_1d(step_, i)

  hput_int(grid_ptr->priority, "priority", comp_dir);
  hput_string(grid_ptr->grid_name, "component grid name", comp_dir);

  hput_string((grid_ptr->grid_type == 1)? "cartesian" : "curvilinear", 
	      "grid type", comp_dir);

  dimension_ = create_int_array_1d(2);
  dimension(1) = grid_ptr->r_dim;
  dimension(2) = grid_ptr->s_dim;
  hput_int_array_1d(dimension_, "dimension", comp_dir);
  dimension_ = delete_int_array_1d(dimension_);

  periodicity_ = create_int_array_1d(2);
  periodicity(1) = grid_ptr->r_period;
  periodicity(2) = grid_ptr->s_period;
  hput_int_array_1d(periodicity_, "periodicity", comp_dir);
  periodicity_ = delete_int_array_1d(periodicity_);

  step_ = create_real_array_1d(2);
  step(1) = grid_ptr->r_step;
  step(2) = grid_ptr->s_step;
  hput_real_array_1d(step_, "step", comp_dir);
  step_ = delete_real_array_1d(step_);

  hput_int_array_2d(grid_ptr->range_ptr, "range", comp_dir);
  hput_int_array_2d(grid_ptr->bc_ptr, "boundary condition", comp_dir);
  hput_int_array_2d(grid_ptr->curve_ptr, "curve label", comp_dir);

/* make a real array 2d to put the bounding box in */
  bbox_ = create_real_array_2d(2,2);
  bbox(1,1) = grid_ptr->x_min;
  bbox(2,1) = grid_ptr->x_max;
  bbox(1,2) = grid_ptr->y_min;
  bbox(2,2) = grid_ptr->y_max;
/* save the bounding box */
  hput_real_array_2d(bbox_, "bounding box", comp_dir);
  bbox_ = delete_real_array_2d( bbox_ );

/* Topology information */
  hput_int(grid_ptr->topology, "topology", comp_dir);
  hput_int(grid_ptr->branch_direction, "branch_direction", comp_dir);
  hput_int(grid_ptr->first_point, "first_point", comp_dir);
  hput_int(grid_ptr->last_point,  "last_point",  comp_dir);

/* large arrays */
  hput_int_array_2d(grid_ptr->flag_ptr, "flag", comp_dir);
  hput_real_array_2d(grid_ptr->x_ptr, "x", comp_dir);
  hput_real_array_2d(grid_ptr->y_ptr, "y", comp_dir);

/* only save the jacobian if save_jacobian == 1 */
  if (save_jacobian == 1){
    hput_int(1, "save jacobian", comp_dir);
    xr_ptr = create_real_array_2d(grid_ptr->r_dim, grid_ptr->s_dim);
    xs_ptr = create_real_array_2d(grid_ptr->r_dim, grid_ptr->s_dim);
    yr_ptr = create_real_array_2d(grid_ptr->r_dim, grid_ptr->s_dim);
    ys_ptr = create_real_array_2d(grid_ptr->r_dim, grid_ptr->s_dim);

    compute_jacobian( grid_ptr, xr_ptr, xs_ptr, yr_ptr, ys_ptr );

    hput_real_array_2d( xr_ptr, "xr", comp_dir ); 
    hput_real_array_2d( xs_ptr, "xs", comp_dir ); 
    hput_real_array_2d( yr_ptr, "yr", comp_dir ); 
    hput_real_array_2d( ys_ptr, "ys", comp_dir ); 

    xr_ptr = delete_real_array_2d( xr_ptr );
    xs_ptr = delete_real_array_2d( xs_ptr );
    yr_ptr = delete_real_array_2d( yr_ptr );
    ys_ptr = delete_real_array_2d( ys_ptr );
  }
  else{
    hput_int(0, "save jacobian", comp_dir);
  }
/* # interpolation points */
  hput_int(grid_ptr->n_interp, "n_interp", comp_dir);

  if (grid_ptr->n_interp > 0){
/* make 2-d arrays of the interpolation list */
    donor_par_ = create_real_array_2d(2, grid_ptr->n_interp);
    donor_point_ = create_int_array_2d(3, grid_ptr->n_interp);
    inter_point_ = create_int_array_2d(2, grid_ptr->n_interp);

    for (interp_ptr = grid_ptr->last_interp, i=1; interp_ptr != NULL; 
	 interp_ptr = interp_ptr->prev, i++){ /* cool syntax for i ! */
      donor_par(1, i)   = interp_ptr->r_loc;
      donor_par(2, i)   = interp_ptr->s_loc;

      donor_point(1, i) = interp_ptr->i_loc;
      donor_point(2, i) = interp_ptr->j_loc;
      donor_point(3, i) = interp_ptr->grid_loc->priority;

      inter_point(1, i) = interp_ptr->i_point;
      inter_point(2, i) = interp_ptr->j_point;
    }

/* save the interpolation information */
    hput_real_array_2d( donor_par_, "donor parameter", comp_dir );
    hput_int_array_2d( donor_point_, "donor point", comp_dir );
    hput_int_array_2d( inter_point_, "interpolation point", comp_dir );

/* release the temporary storage */
    delete_real_array_2d( donor_par_ );
    delete_int_array_2d( donor_point_ );
    delete_int_array_2d( inter_point_ );
  }

#undef dimension
}


#define xr(i,j)     compute_index_2d(xr_ptr,i,j)
#define xs(i,j)     compute_index_2d(xs_ptr,i,j)
#define yr(i,j)     compute_index_2d(yr_ptr,i,j)
#define ys(i,j)     compute_index_2d(ys_ptr,i,j)

static void 
compute_jacobian(component_grid *grid_ptr, 
		 real_array_2d *xr_ptr, real_array_2d *xs_ptr, 
		 real_array_2d *yr_ptr, real_array_2d *ys_ptr){
  int i,j;
  grid_point gp;

/* compute the jacobian by calling the mapping function, if it is available */
/*   if (grid_ptr->grid_mapping != NULL){ */
/*     for (i=1; i<=grid_ptr->r_dim; i++){ */
/*       for (j=1; j<=grid_ptr->s_dim; j++){ */
/* 	gp.r = (i-range(1,1)) * grid_ptr->r_step; */
/* 	gp.s = (j-range(1,2)) * grid_ptr->s_step; */
/* 	forward_grid_mapping( &gp, grid_ptr ); */
/* 	xr(i,j) = gp.xr; */
/* 	xs(i,j) = gp.xs; */
/* 	yr(i,j) = gp.yr; */
/* 	ys(i,j) = gp.ys; */
/*       } */
/*     } */
/*   } */

/* otherwise, approximate the jacobian by divided differences */
/*  else{*/

/* always approximate the jacobian by 2nd order divided differences */
/* left side */
    i=1;
    j=1;
    xr(i,j) = 0.5*(-3.0*x(i,j) + 4.0*x(i+1,j) - x(i+2,j))/grid_ptr->r_step;
    yr(i,j) = 0.5*(-3.0*y(i,j) + 4.0*y(i+1,j) - y(i+2,j))/grid_ptr->r_step;
    xs(i,j) = 0.5*(-3.0*x(i,j) + 4.0*x(i,j+1) - x(i,j+2))/grid_ptr->s_step;
    ys(i,j) = 0.5*(-3.0*y(i,j) + 4.0*y(i,j+1) - y(i,j+2))/grid_ptr->s_step;

    for (j=2; j<grid_ptr->s_dim; j++){
      xr(i,j) = 0.5*(-3.0*x(i,j) + 4.0*x(i+1,j) - x(i+2,j))/grid_ptr->r_step;
      yr(i,j) = 0.5*(-3.0*y(i,j) + 4.0*y(i+1,j) - y(i+2,j))/grid_ptr->r_step;
      xs(i,j) = 0.5*(x(i,j+1) - x(i,j-1))/grid_ptr->s_step;
      ys(i,j) = 0.5*(y(i,j+1) - y(i,j-1))/grid_ptr->s_step;
    }

    j=grid_ptr->s_dim;
    xr(i,j) = 0.5*(-3.0*x(i,j) + 4.0*x(i+1,j) - x(i+2,j) )/grid_ptr->r_step;
    yr(i,j) = 0.5*(-3.0*y(i,j) + 4.0*y(i+1,j) - y(i+2,j) )/grid_ptr->r_step;
    xs(i,j) = 0.5*( 3.0*x(i,j) - 4.0*x(i,j-1) + x(i,j-2) )/grid_ptr->s_step;
    ys(i,j) = 0.5*( 3.0*y(i,j) - 4.0*y(i,j-1) + y(i,j-2) )/grid_ptr->s_step;

/* interior in the i-direction */
    for (i=2; i<grid_ptr->r_dim; i++){

      j = 1;
      xr(i,j) = 0.5*(x(i+1,j) - x(i-1,j))/grid_ptr->r_step;
      yr(i,j) = 0.5*(y(i+1,j) - y(i-1,j))/grid_ptr->r_step;
      xs(i,j) = 0.5*(-3.0*x(i,j) + 4.0*x(i,j+1) - x(i,j+2))/grid_ptr->s_step;
      ys(i,j) = 0.5*(-3.0*y(i,j) + 4.0*y(i,j+1) - y(i,j+2))/grid_ptr->s_step;

      for (j=2; j<grid_ptr->s_dim; j++){
	xr(i,j) = 0.5*(x(i+1,j) - x(i-1,j))/grid_ptr->r_step;
	xs(i,j) = 0.5*(x(i,j+1) - x(i,j-1))/grid_ptr->s_step;
	yr(i,j) = 0.5*(y(i+1,j) - y(i-1,j))/grid_ptr->r_step;
	ys(i,j) = 0.5*(y(i,j+1) - y(i,j-1))/grid_ptr->s_step;
      }

      j=grid_ptr->s_dim;
      xr(i,j) = 0.5*(x(i+1,j) - x(i-1,j))/grid_ptr->r_step;
      yr(i,j) = 0.5*(y(i+1,j) - y(i-1,j))/grid_ptr->r_step;
      xs(i,j) = 0.5*( 3.0*x(i,j) - 4.0*x(i,j-1) + x(i,j-2) )/grid_ptr->s_step;
      ys(i,j) = 0.5*( 3.0*y(i,j) - 4.0*y(i,j-1) + y(i,j-2) )/grid_ptr->s_step;
    }

/* right side */
    i=grid_ptr->r_dim;
    j=1;
    xr(i,j) = 0.5*( 3.0*x(i,j) - 4.0*x(i-1,j) + x(i-2,j))/grid_ptr->r_step;
    yr(i,j) = 0.5*( 3.0*y(i,j) - 4.0*y(i-1,j) + y(i-2,j))/grid_ptr->r_step;
    xs(i,j) = 0.5*(-3.0*x(i,j) + 4.0*x(i,j+1) - x(i,j+2))/grid_ptr->s_step;
    ys(i,j) = 0.5*(-3.0*y(i,j) + 4.0*y(i,j+1) - y(i,j+2))/grid_ptr->s_step;

    for (j=2; j<grid_ptr->s_dim; j++){
      xr(i,j) = 0.5*( 3.0*x(i,j) - 4.0*x(i-1,j) + x(i-2,j))/grid_ptr->r_step;
      yr(i,j) = 0.5*( 3.0*y(i,j) - 4.0*y(i-1,j) + y(i-2,j))/grid_ptr->r_step;
      xs(i,j) = 0.5*(x(i,j+1) - x(i,j-1))/grid_ptr->s_step;
      ys(i,j) = 0.5*(y(i,j+1) - y(i,j-1))/grid_ptr->s_step;
    }

    j=grid_ptr->s_dim;
    xr(i,j) = 0.5*( 3.0*x(i,j) - 4.0*x(i-1,j) + x(i-2,j))/grid_ptr->r_step;
    yr(i,j) = 0.5*( 3.0*y(i,j) - 4.0*y(i-1,j) + y(i-2,j))/grid_ptr->r_step;
    xs(i,j) = 0.5*( 3.0*x(i,j) - 4.0*x(i,j-1) + x(i,j-2) )/grid_ptr->s_step;
    ys(i,j) = 0.5*( 3.0*y(i,j) - 4.0*y(i,j-1) + y(i,j-2) )/grid_ptr->s_step;

/*  }*/
  return;
}

#undef xr
#undef xs
#undef yr
#undef ys

void 
save_grid_plot3d(generic_mapping *map_ptr, FILE *fp, int i_min, int i_max, 
		 int j_min, int j_max ){
  int i, j, tmp;
  real dr, ds;
  grid_point gp;

  dr = 1.0/ (map_ptr->r_points - 1);
  ds = 1.0/ (map_ptr->s_points - 1);

/* periodic in r */
  if (i_min > i_max && map_ptr->r_period){
    printf("The grid seems to be periodic in r\n");
    if (j_min > j_max){
      tmp = j_min;
      j_min = j_max;
      j_max = tmp;
    }
/* formatted plot3d file with one grid in the file */
    fprintf( fp, "%i %i %i\n", map_ptr->r_points - i_min + i_max, j_max - j_min + 1, 1 );

/* save the x coordinates */
    for (j = j_min; j <= j_max; j++){
      for (i = i_min; i < map_ptr->r_points; i++){
	gp.r = (i-1)*dr; gp.s = (j-1)*ds;
	forward_mapping( &gp, map_ptr );
	fprintf( fp, "%.18G\n", gp.x);
      }
      for (i = 1; i <= i_max; i++){
	gp.r = (i-1)*dr; gp.s = (j-1)*ds;
	forward_mapping( &gp, map_ptr );
	fprintf( fp, "%.18G\n", gp.x);
      }
    }

/* save the y coordinates */
    for (j = j_min; j <= j_max; j++){
      for (i = i_min; i < map_ptr->r_points; i++){
	gp.r = (i-1)*dr; gp.s = (j-1)*ds;
	forward_mapping( &gp, map_ptr );
	fprintf( fp, "%.18G\n", gp.y);
      }
      for (i = 1; i <= i_max; i++){
	gp.r = (i-1)*dr; gp.s = (j-1)*ds;
	forward_mapping( &gp, map_ptr );
	fprintf( fp, "%.18G\n", gp.y);
      }
    }
    
  } /* end periodic in r */

/* periodic in s */
  else if (j_min > j_max && map_ptr->s_period){
    printf("The grid seems to be periodic in s. Sorry, not implemented yet\n");
  }

/* non-periodic */
  else{
/* formatted plot3d file with one grid in the file */
    fprintf( fp, "%i %i %i\n", map_ptr->r_points, map_ptr->s_points, 1 );

/* save the x coordinates */
    for (j = 1; j <= map_ptr->s_points; j++){
      for (i = 1; i <= map_ptr->r_points; i++){
	gp.r = (i-1)*dr; gp.s = (j-1)*ds;
	forward_mapping( &gp, map_ptr );
	fprintf( fp, "%.18G\n", gp.x);
      }
    }

/* save the y coordinates */
    for (j = 1; j <= map_ptr->s_points; j++){
      for (i = 1; i <= map_ptr->r_points; i++){
	gp.r = (i-1)*dr; gp.s = (j-1)*ds;
	forward_mapping( &gp, map_ptr );
	fprintf( fp, "%.18G\n", gp.y);
      }
    }
  } /* end non-periodic */

/* done */
}

