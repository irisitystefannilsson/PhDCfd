#ifdef THE_GAR 

#include "mappings.h"
#define MAX_POINTS 1000

/* Theodorsen-Garrick mapping */

/* private member functions */
static void 
compute_tg_mapping( tg_grid_info *info, generic_mapping *grid_ptr );
static void 
set_tg_mapping(input_output *io_ptr, generic_mapping *grid_ptr, generic_curve *first_c);
static int 
tg_mapping( grid_point *gp_ptr, void *info);
static void *
cleanup_tg_mapping( void *data_ptr );
static void 
write_tg_mapping( int32 dir, void *tg_data_ptr );
static void *
copy_tg_mapping( void *mapping_data_ptr );
static real_array_2d *
create_naca0012(void);
/* end private member function declaration */

#define off(i,j) compute_index_2d(info->coordinate_list,i,j)

static void 
set_tg_mapping(input_output *io_ptr, generic_mapping *grid_ptr, generic_curve *first_c){
  char prompt[80];
  int icom, quit, replot;
  tg_grid_info *info;
  int no_indent=0, n;
  real xte,yte,x1,y1,slope;

#include "tg_grid_com.h"

  sprintf(prompt, "%s: theodorsen-garrick mapping>", grid_ptr->grid_name);

/* cast the data pointer to the right type */
  info = (tg_grid_info *) grid_ptr->mapping_data_ptr;

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
      replot = FALSE;
    }

    icom = get_command(io_ptr, prompt, COMMAND, LAST_COM, 
		       LEVEL, SAVE_ON_COPY, ARGUMENT);

    switch (icom) {

    case TE_RADIUS:
      info->tgcoeff.rte = real_max( 0.0, 
				   get_real( io_ptr, 
					    "trailing edge radius >= 0 (but small!): ", 
					    info->tgcoeff.rte, no_indent));

      n = info->coordinate_list->n2;
#define coord_list(i,j) compute_index_2d(info->coordinate_list,i,j)
      xte = coord_list(1,1);
      yte = coord_list(2,1);
      x1 = (coord_list(1,2)+coord_list(1,n-1))/2.0;
      y1 = (coord_list(2,2)+coord_list(2,n-1))/2.0;
#undef coord_list
      slope = (y1-yte)/(x1-xte);

      info->tgcoeff.dtau = -0.52;
      info->tgcoeff.dxte =-50.0*info->tgcoeff.rte;
      info->tgcoeff.dyte = slope*info->tgcoeff.dxte;
      compute_tg_mapping( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      replot = TRUE;
      break;

    case S_SCALING:
      info->s_scaling = real_max(1.e-7, 
				 get_real( io_ptr, "Scaling of s-direction > 0: ", 
					  info->s_scaling, no_indent));
      compute_tg_mapping( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      replot = TRUE;
      break;

    case QUALITY:
/* plot the mapping with the directional arrows, all grid lines and quality */
/* information which is computed by a call to grid_quality from plot_generic_mapping */
      PL_window(1, grid_ptr->xyab);
      PL_erase(1);
      PL_start_plot(1);
      plot_generic_mapping( grid_ptr, 2+8+256 );
      PL_end_plot();
      replot = FALSE;
      break;

    case R_STRETCH:
/* r-stretching */
      if (info->r_stretch == NULL){
	info->r_stretch = choose_stretching( io_ptr, NULL, &grid_ptr->r_points );
/* the inverse is no longer valid */
	grid_ptr->inverse_known = 0; grid_ptr->inverse_grid_mapping = NULL;
/* update the plot mode */
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	replot = TRUE;
      }
      else{
	if ( set_stretching( io_ptr, info->r_stretch, NULL, &grid_ptr->r_points) ){
	  grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	  replot = TRUE;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->r_stretch ) );
	}
      }
      break;

    case NO_R_STRETCH:
/* no-r-stretching */
      if (info->r_stretch != NULL){
	info->r_stretch = delete_generic_stretching( info->r_stretch );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	replot = TRUE;
      }
      else{
	printf("There is no stretching in the s-direction.\n");
      }
      break;

    case S_STRETCH:
/* s-stretching */
      if (info->s_stretch == NULL){
	info->s_stretch = choose_stretching( io_ptr, NULL, &grid_ptr->s_points );
/* the inverse is no longer valid */
	grid_ptr->inverse_known = 0; grid_ptr->inverse_grid_mapping = NULL;
/* update the plot mode */
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	replot = TRUE;
      }
      else{
	if ( set_stretching( io_ptr, info->s_stretch, NULL, &grid_ptr->s_points) ){
	  grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	  replot = TRUE;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->s_stretch ) );
	}
      }
      break;

    case NO_S_STRETCH:
/* no-s-stretching */
      if (info->s_stretch != NULL){
	info->s_stretch = delete_generic_stretching( info->s_stretch );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	replot = TRUE;
      }
      else{
	printf("There is no stretching in the s-direction.\n");
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
      replot = TRUE;
      break;

    case S_LINES:
/* s-lines */
      grid_ptr->s_points = 
	int_max(2, get_int(io_ptr, 
			   "Enter number of gridlines in s>= 2: ", grid_ptr->s_points,
			   no_indent));

      grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      replot = TRUE;
      break;

    case PLOT_MODE:
/* plot-mode */
      set_grid_plot_mode( io_ptr, grid_ptr );
      break;

    case SHOW:
/* show */
      show_mapping_parameters( grid_ptr );
      printf("Trailing edge radius parameter = %e.\n", info->tgcoeff.rte);
      printf("%s in the r-direction and\n%s in the s-direction.\n", 
	     stretching_name( info->r_stretch ), 
	     stretching_name( info->s_stretch ) );
      replot = 0;
      break;

    case HELP:
/* help */
      while ( (icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				  COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1 );
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
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

}

static int 
tg_mapping( grid_point *gp_ptr, void *data_ptr){
  tg_grid_info *info;
  real r_stretch, r_ratio, s_stretch, s_ratio, dummy;
  real xp, xm, yp, ym, rp, rm, sp, sm;
  const real dr=0.01, ds=0.01;

/* cast the data to the right type */
  info = (tg_grid_info *) data_ptr;

/* r-stretching */
  uniform_to_stretch( gp_ptr->r, info->r_stretch, &r_stretch, 
		     &r_ratio, &dummy );

/* s-stretching */
  uniform_to_stretch( gp_ptr->s, info->s_stretch, &s_stretch, 
		     &s_ratio, &dummy );

/* scale in s */
  s_stretch *= info->s_scaling;

  mktggrid(&(gp_ptr->x),&(gp_ptr->y),r_stretch,s_stretch,info->tgcoeff);

/* approximation of the Jacobian... */
  rp = real_min(1.0, r_stretch+dr);
  rm = real_max(0.0, r_stretch-dr);
  mktggrid(&xp,&yp,rp,s_stretch,info->tgcoeff);
  mktggrid(&xm,&ym,rm,s_stretch,info->tgcoeff);
  gp_ptr->xr= r_ratio * (xp-xm)/(rp-rm);
  gp_ptr->yr= r_ratio * (yp-ym)/(rp-rm);

  sp = real_min(info->s_scaling, s_stretch+ds);
  sm = real_max(0.0, s_stretch-ds);
  mktggrid(&xp,&yp,r_stretch,sp,info->tgcoeff);
  mktggrid(&xm,&ym,r_stretch,sm,info->tgcoeff);
  gp_ptr->xs= s_ratio * info->s_scaling * (xp-xm)/(sp-sm);
  gp_ptr->ys= s_ratio * info->s_scaling * (yp-ym)/(sp-sm);

  return 0;
}

static void 
compute_tg_mapping( tg_grid_info *info, generic_mapping *grid_ptr ){
  grid_point gp;
  const real dr=0.05;
/* compute the coefficients of the mapping */
  mktgcoeff(info->coordinate_list->n2,
	    info->coordinate_list->arrayptr,
	    &(info->tgcoeff));

/* compute the bounding box */
  grid_ptr->x_min_0 = 1.e10;
  grid_ptr->x_max_0 = -1.e10;
  grid_ptr->y_min_0 = 1.e10;
  grid_ptr->y_max_0 = -1.e10;

/* the foil is at s=0 and the outer boundary at s=1.0 */
  for (gp.s = 0.0; gp.s <= 1.05; gp.s +=1.0)
    for (gp.r = 0.0; gp.r <= 1.0+0.5*dr; gp.r+=dr){
      tg_mapping( &gp, info );
      if (gp.x < grid_ptr->x_min_0) grid_ptr->x_min_0 = gp.x;
      if (gp.x > grid_ptr->x_max_0) grid_ptr->x_max_0 = gp.x;
      if (gp.y < grid_ptr->y_min_0) grid_ptr->y_min_0 = gp.y;
      if (gp.y > grid_ptr->y_max_0) grid_ptr->y_max_0 = gp.y;
    }

  reset_view_port( grid_ptr );
}


static void 
write_tg_mapping( int32 dir, void *tg_data_ptr ){
  tg_grid_info *info;
  int32 r_stretch_dir, s_stretch_dir;

  info = (tg_grid_info *) tg_data_ptr;

  hput_real( info->xmin, "xmin", dir );
  hput_real( info->xmax, "xmax", dir );
  hput_real( info->ymin, "ymin", dir );
  hput_real( info->ymax, "ymax", dir );
  hput_real( info->s_scaling, "s_scaling", dir );
  hput_real_array_2d(info->coordinate_list, "coordinate_list", dir);
  hput_real( info->tgcoeff.te, "te", dir );
  hput_real( info->tgcoeff.beta, "beta", dir );
  hput_real( info->tgcoeff.z0.re, "z0.re", dir );
  hput_real( info->tgcoeff.z0.im, "z0.im", dir );
  hput_real( info->tgcoeff.z1.re, "z1.re", dir );
  hput_real( info->tgcoeff.z1.im, "z1.im", dir );
  hput_real( info->tgcoeff.c.re, "c.re", dir );
  hput_real( info->tgcoeff.c.im, "c.im", dir );
  hput_real( info->tgcoeff.rte, "rte", dir );
  hput_real( info->tgcoeff.dtau, "dtau", dir );
  hput_real( info->tgcoeff.dxte, "dxte", dir );
  hput_real( info->tgcoeff.dyte, "dyte", dir );


/* r-stretching */
  if ((r_stretch_dir = create_dir("r-stretching", "r-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( r_stretch_dir, info->r_stretch );
/* release the sub-directory */
    Vdetach(r_stretch_dir);
  }
  else
    printf("Error: write_tg_mapping: unable to create a directory for "
	   "`r-stretching'.\n");

/* s-stretching */
  if ((s_stretch_dir = create_dir("s-stretching", "s-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( s_stretch_dir, info->s_stretch );
/* release the sub-directory */
    Vdetach(s_stretch_dir);
  }
  else
    printf("Error: save_curves_mappings: unable to create a directory for "
	   "`s-stretching'.\n");
}

static void *
cleanup_tg_mapping( void *data_ptr ){
  tg_grid_info *info;

/* cast the data to the right type */
  info = (tg_grid_info *) data_ptr;

  if (info != NULL) {
    info->r_stretch = delete_generic_stretching( info->r_stretch );

    info->s_stretch = delete_generic_stretching( info->s_stretch );

    delete_real_array_2d(info->coordinate_list);
    free(info->tgcoeff.name);
    free(info->tgcoeff.a);
    free(info->tgcoeff.b);
    free(info);
  }
  return NULL;
}

tg_grid_info *
init_tg_mapping(input_output *io_ptr,generic_mapping *grid_ptr){

  tg_grid_info *info;
  FILE *fp;
  char *file_name;
  int i,istat, npoints, n;
  real xte, yte, x1, y1, slope;
  const int save_command=1;
  
/* this must be deallocated when the mapping is deleted */
  if ((info = (tg_grid_info *) 
    malloc (sizeof (tg_grid_info))) == NULL)
    printf("memory error in init_tg_mapping\n");

/* get an offset file and read the coordinates. */
/* Use the same format as for cubic splines. */
  if ((fp = open_ascii_file( io_ptr, "Enter wing section ascii file: ", 
			    "naca0012.spl", &file_name, 'r', 4, 
			    save_command)) != NULL){
/* Read number of offset points in the file */
    fscanf( fp, "%i", &npoints);

/* check npoints */
    if (npoints < 1 || npoints > 10000){
      printf("Error in init_tg_mapping: unreasonable number of nodepoints = %i\n", 
	     npoints);
      return NULL;
    }

    if ((info->coordinate_list = create_real_array_2d(2,npoints)) == NULL){
      printf("memory error in init_tg_mapping\n");
      return NULL;
    }

/* Read the offset points */
    for(i=1;(istat = fscanf( fp, "%lf%lf", &(off(1,i)),&(off(2,i))) != EOF) && 
	i<=npoints;i++);

/* check the data */
    if (i != npoints+1 || istat == EOF || off(1,1) != off(1,npoints)) {
      printf("Error in init_tg_mapping: illegal format in the input file\n");
      delete_real_array_2d(info->coordinate_list);
      return NULL;
    }
    fclose(fp);
  }
/* if the file could not be opened, give a naca0012 profile instead */
  else{
/* set default wing section to NACA 0012 */
    printf("Could not read your file, so instead I give you a NACA 0012\n");
    info->coordinate_list = create_naca0012();
  }

/* noname */
  info->tgcoeff.name = NULL;

/* set default radius for the T-G mapping */
  info->tgcoeff.rte = 5.e-4;

/* default scaling in s */
  info->s_scaling = 1.0;

  n = info->coordinate_list->n2;
#define coord_list(i,j) compute_index_2d(info->coordinate_list,i,j)
  xte = coord_list(1,1);
  yte = coord_list(2,1);
  x1 = (coord_list(1,2)+coord_list(1,n-1))/2.0;
  y1 = (coord_list(2,2)+coord_list(2,n-1))/2.0;
#undef coord_list
  slope = (y1-yte)/(x1-xte);
  
  info->tgcoeff.dtau = -0.52;
  info->tgcoeff.dxte =-50.0*info->tgcoeff.rte;
  info->tgcoeff.dyte = slope*info->tgcoeff.dxte;

/* bounding box (cheating!)*/
  info->xmin = -0.5;
  info->xmax = 1.5;
  info->ymin =-0.5;
  info->ymax = 0.5;

/* no stretching */
  info->r_stretch = NULL;
  info->s_stretch = NULL;

/* set the grid type and the periodicity flag */
  grid_ptr->grid_type= 7; /* Theodorsen-Garrick grid */
  grid_ptr->r_period = 1;
  grid_ptr->s_period = 0;

/* set the number of grid lines */
  grid_ptr->r_points = 50;
  grid_ptr->s_points = 10;

/* save the pointers to the data and the functions for
   evaluating the mapping */
  grid_ptr->grid_mapping = tg_mapping;
  grid_ptr->inverse_known = 0;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->mapping_data_ptr = info;
  grid_ptr->cleanup_mapping_data = cleanup_tg_mapping;
  grid_ptr->set_mapping = set_tg_mapping;
  grid_ptr->write_specific_mapping = write_tg_mapping;
  grid_ptr->copy_mapping_data = copy_tg_mapping;

/* Compute the coefficients of the Theodorsen-Garrick mapping.		*/
  compute_tg_mapping( info, grid_ptr );

  return info;
}

static real_array_2d *
create_naca0012(void){
  real_array_2d *naca_;

#define naca(i,j) compute_index_2d(naca_, i, j)

  naca_ = create_real_array_2d(2, 35);
  naca(1, 1) = 1.00;   naca(2, 1) = 0.00;
  naca(1, 2) = 0.95;   naca(2, 2) = 0.00807;
  naca(1, 3) = 0.90;   naca(2, 3) = 0.01448;
  naca(1, 4) = 0.80;   naca(2, 4) = 0.02623;
  naca(1, 5) = 0.70;   naca(2, 5) = 0.03664;
  naca(1, 6) = 0.60;   naca(2, 6) = 0.04563;
  naca(1, 7) = 0.50;   naca(2, 7) = 0.05294;
  naca(1, 8) = 0.40;   naca(2, 8) = 0.05803;
  naca(1, 9) = 0.30;   naca(2, 9) = 0.06002;
  naca(1,10) = 0.25;   naca(2,10) = 0.05941;
  naca(1,11) = 0.20;   naca(2,11) = 0.05737;
  naca(1,12) = 0.15;   naca(2,12) = 0.05345;
  naca(1,13) = 0.10;   naca(2,13) = 0.04683;
  naca(1,14) = 0.075;  naca(2,14) = 0.04200;
  naca(1,15) = 0.05;   naca(2,15) = 0.03555;
  naca(1,16) = 0.025;  naca(2,16) = 0.02615;
  naca(1,17) = 0.0125; naca(2,17) = 0.01894;
  naca(1,18) = 0.00;   naca(2,18) = 0.0000;
  naca(1,19) = 0.0125; naca(2,19) = -0.01894;
  naca(1,20) = 0.025;  naca(2,20) = -0.02615;
  naca(1,21) = 0.05;   naca(2,21) = -0.03555;
  naca(1,22) = 0.075;  naca(2,22) = -0.04200;
  naca(1,23) = 0.10;   naca(2,23) = -0.04683;
  naca(1,24) = 0.15;   naca(2,24) = -0.05345;
  naca(1,25) = 0.20;   naca(2,25) = -0.05737;
  naca(1,26) = 0.25;   naca(2,26) = -0.05941;
  naca(1,27) = 0.30;   naca(2,27) = -0.06002;
  naca(1,28) = 0.40;   naca(2,28) = -0.05803;
  naca(1,29) = 0.50;   naca(2,29) = -0.05294;
  naca(1,30) = 0.60;   naca(2,30) = -0.04563;
  naca(1,31) = 0.70;   naca(2,31) = -0.03664;
  naca(1,32) = 0.80;   naca(2,32) = -0.02623;
  naca(1,33) = 0.90;   naca(2,33) = -0.01448;
  naca(1,34) = 0.95;   naca(2,34) = -0.00807;
  naca(1,35) = 1.00;   naca(2,35) = -0.00;

  return naca_;
#undef naca
}

int
read_tg_mapping(int32 dir, generic_mapping *grid_ptr, generic_curve *first_c ){
  tg_grid_info *info;
  int32 stretch_dir;

  info = (tg_grid_info *) malloc( sizeof(tg_grid_info) );

  hget_real( &(info->xmin), "xmin", dir );
  hget_real( &(info->xmax), "xmax", dir );
  hget_real( &(info->ymin), "ymin", dir );
  hget_real( &(info->ymax), "ymax", dir );

  if (!hget_real( &(info->s_scaling), "s_scaling", dir )){
    printf("Warning: Assigning default value s_scaling=1.0.\n");
    info->s_scaling = 1.0;
  }

  hget_real( &(info->tgcoeff.te),    "te",    dir );
  hget_real( &(info->tgcoeff.beta),  "beta",  dir );
  hget_real( &(info->tgcoeff.z0.re), "z0.re", dir );
  hget_real( &(info->tgcoeff.z0.im), "z0.im", dir );
  hget_real( &(info->tgcoeff.z1.re), "z1.re", dir );
  hget_real( &(info->tgcoeff.z1.im), "z1.im", dir );
  hget_real( &(info->tgcoeff.c.re),  "c.re",  dir );
  hget_real( &(info->tgcoeff.c.im),  "c.im",  dir );
  hget_real( &(info->tgcoeff.rte),   "rte",   dir );
  hget_real( &(info->tgcoeff.dtau),  "dtau",  dir );
  hget_real( &(info->tgcoeff.dxte),  "dxte",  dir );
  hget_real( &(info->tgcoeff.dyte),  "dyte",  dir );

  info->coordinate_list = hget_real_array_2d( "coordinate_list", dir);

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

/* Compute the coefficients of the Theodorsen-Garrick mapping.		*/
  compute_tg_mapping( info, grid_ptr );
/* Assign pointers */
  grid_ptr->mapping_data_ptr = (void *) info;
  grid_ptr->grid_mapping = tg_mapping;
  grid_ptr->inverse_known = 0;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->cleanup_mapping_data = cleanup_tg_mapping;
  grid_ptr->set_mapping = set_tg_mapping;
  grid_ptr->write_specific_mapping = write_tg_mapping;
  grid_ptr->copy_mapping_data = copy_tg_mapping;

  return OK;
}


static void *
copy_tg_mapping( void *mapping_data_ptr ){
  tg_grid_info *info, *old_info;
  int           nc,ne;
  
  info = (tg_grid_info *) malloc( sizeof(tg_grid_info) );

  old_info = (tg_grid_info *)mapping_data_ptr;

  info->xmin = old_info->xmin;
  info->xmax = old_info->xmax;
  info->ymin = old_info->ymin;
  info->ymax = old_info->ymax;

  info->s_scaling = old_info->s_scaling;

/* copy the coordinates */
  info->coordinate_list = create_real_array_2d(
    old_info->coordinate_list->n1,
    old_info->coordinate_list->n2);
  ne = old_info->coordinate_list->n1*old_info->coordinate_list->n2;
  memcpy(info->coordinate_list->arrayptr,
         old_info->coordinate_list->arrayptr,ne*sizeof(real));
/* copy the parameters */
  info->tgcoeff.name = (char *)malloc(strlen(old_info->tgcoeff.name));
  strcpy(info->tgcoeff.name,old_info->tgcoeff.name);
  info->tgcoeff.npoints = old_info->tgcoeff.npoints;
  nc = 2*info->tgcoeff.npoints+1;
  info->tgcoeff.a = (real *)malloc(nc*sizeof(real));
  info->tgcoeff.b = (real *)malloc(nc*sizeof(real));
  memcpy(info->tgcoeff.a,old_info->tgcoeff.a,nc*sizeof(real));
  memcpy(info->tgcoeff.b,old_info->tgcoeff.b,nc*sizeof(real));
  info->tgcoeff.te = old_info->tgcoeff.te;
  info->tgcoeff.beta = old_info->tgcoeff.beta;
  info->tgcoeff.z0.re = old_info->tgcoeff.z0.re;
  info->tgcoeff.z0.im = old_info->tgcoeff.z0.im;
  info->tgcoeff.z1.re = old_info->tgcoeff.z1.re;
  info->tgcoeff.z1.im = old_info->tgcoeff.z1.im;
  info->tgcoeff.c.re = old_info->tgcoeff.c.re;
  info->tgcoeff.c.im = old_info->tgcoeff.c.im;
  info->tgcoeff.rte = old_info->tgcoeff.rte;
  info->tgcoeff.dtau = old_info->tgcoeff.dtau;
  info->tgcoeff.dxte = old_info->tgcoeff.dxte;
  info->tgcoeff.dyte = old_info->tgcoeff.dyte;
/* make a copy of the stretching functions */
  info->r_stretch = copy_stretching( old_info->r_stretch );
  info->s_stretch = copy_stretching( old_info->s_stretch );

  return (void *)info;
}

#else /* add a dummy function to keep the compiler busy */

static void 
happy_compiler(void){
}

#endif
