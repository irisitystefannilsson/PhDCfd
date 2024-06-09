#include "mappings.h"
#include "log_cosh.h"

/* private member functions */
static void alloc_width_arrays( normal_curve_data *info );
static void free_width_arrays( normal_curve_data *info );
static void compute_normal_curve(normal_curve_data *info,
				 generic_mapping *grid_ptr);
static void 
eval_width(real w1, real *slope, real *width_jmp, real *r_node, 
		       real *v0, real *u0, int n, real width_trans, 
		       real r, real *w, real *w_r);
static void set_normal_curve_mapping(input_output *io_ptr, 
				     generic_mapping *grid_ptr, 
				  generic_curve *first_c);
static void 
write_normal_curve( int32 dir, void *normal_curve_data_ptr );
static int normal_curve_mapping( grid_point *gp_ptr, void *data_ptr);
static void *cleanup_normal_curve(void *data_ptr );
static void *copy_normal_curve( void *mapping_data_ptr );
/* mappings for smoothed polygon grids */

static void 
set_normal_curve_mapping(input_output *io_ptr, generic_mapping *grid_ptr, 
				  generic_curve *first_c){
  char prompt[80], question[80];
  const real eps=1.e-4;
  int icom, quit, i, replot;
  int level=2, no_indent=0;
  real default_width, *r_nodes;
  normal_curve_data *info;

#include "normal_curve_com.h"

  sprintf(prompt, "%s: normal-curve mapping>", grid_ptr->grid_name);

/* cast the data pointer to the right type */
  info = (normal_curve_data *) grid_ptr->mapping_data_ptr;

  quit = 0;
  replot = 1;
/* recompute the view port because it is possible that the curve has been changed */
  compute_normal_curve( info, grid_ptr );
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
/* change the curve */
      set_curve( io_ptr, info->curve_ptr );
      compute_normal_curve( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 1:
/* flip curve parametrization */
      info->flip = !info->flip;
      compute_normal_curve( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 2:
/* constant-width */
      if (!(info->is_constant)){
/* delete the data structure for the variable width */
	free_width_arrays( info );
      }
      info->width1 = 
	real_max(eps, get_real(io_ptr, "Constant width > 0: ", info->width1, 
			       no_indent));
      info->is_constant = 1;
      compute_normal_curve( info, grid_ptr );
      PL_window(1, grid_ptr->xyab);
      break;

    case 3:
/* variable width */
      if (!(info->is_constant)){
/* remove old arrays */
	free_width_arrays( info );
      }

/* call the curve module to get a copy of the parameter values at the node */
/* points of the curve and the number of node points. */
      r_nodes = copy_nodes( info->curve_ptr, &(info->n) );
      alloc_width_arrays( info );
/* unfortunately, we cannot used r_nodes directly, because we need one padding */
/* point outside each end point */
      for (i=1; i<= info->n; i++)
	info->r[i] = r_nodes[i];
/* free the r_nodes vector */
      free_vector( r_nodes, 1, info->n );

/* read the width */
      default_width = info->width1;
      if (info->flip){
	for (i=1; i<info->n; i++){
	  sprintf(question, "interval %i: start width > 0: ", i);
	  info->width_r[info->n-i] = get_real( io_ptr, question, default_width, 
					      level+1);
	  sprintf(question, "interval %i: end width > 0: ", i);
	  info->width_l[info->n-i] = get_real( io_ptr, question, default_width, 
					      level+1);
	}
      }
      else{
	for (i=1; i<info->n; i++){
	  sprintf(question, "interval %i: start width > 0: ", i);
	  info->width_l[i] = get_real( io_ptr, question, default_width, level+1);
	  sprintf(question, "interval %i: end width > 0: ", i);
	  info->width_r[i] = get_real( io_ptr, question, default_width, level+1);
	}
      }
      info->width_trans = 5.0;
/* flag that variable width is used */
      info->is_constant = 0;

/* compute the mapping */
      compute_normal_curve( info, grid_ptr ); 
      PL_window(1, grid_ptr->xyab);
      break;

    case 4:
/* modify variable width */
      if (info->is_constant){
	printf("You cant modify the varaible width because the width is constant.\n");
	replot = 0;
      }
      else{
	i = get_int(io_ptr , "Which interval: ", 1, no_indent);
	if (i < 1 || i >= info->n){
	  printf("That interval does not exist.\n");
	  replot = 0;
	}
	else{
	  if (info->flip){
	    sprintf(question, "interval %i: start width > 0: ", i);
	    info->width_r[info->n-i] = real_max(eps, get_real( io_ptr, question, 
						info->width_r[info->n-i], level));
	    sprintf(question, "interval %i: end width > 0: ", i);
	    info->width_l[info->n-i] = real_max(eps, get_real( io_ptr, question, 
						info->width_l[info->n-i], level));
	  }
	  else{
	    sprintf(question, "interval %i: start width > 0: ", i);
	    info->width_l[i] = real_max(eps, get_real(io_ptr, question, 
						      info->width_l[i], level));
	    sprintf(question, "interval %i: end width > 0: ", i);
	    info->width_r[i] = real_max(eps, get_real(io_ptr, question, 
						      info->width_r[i], level));
	  }
/* compute the mapping */
	  compute_normal_curve( info, grid_ptr ); 
	  PL_window(1, grid_ptr->xyab);
	}
      }
      break;

    case 5:
/* width-change sharpness */
      if (info->is_constant){
	printf("The width change parameter has no meaning because the width\n");
	printf("is constant.\n");
	replot = 0;
      }
      else{
	info->width_trans = 
	  real_max(eps,
		   get_real(io_ptr, 
			    "Sharpness of width changes > 0: ", info->width_trans,
			    no_indent));
	compute_normal_curve( info, grid_ptr );
	PL_window(1, grid_ptr->xyab);
      }
      break;

    case 6:
/* normal (s-) stretching */
      if (info->normal_stretch == NULL){
	info->normal_stretch = choose_stretching( io_ptr, NULL, 
						 &grid_ptr->s_points );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	if ( set_stretching( io_ptr, info->normal_stretch, NULL, 
			     &grid_ptr->s_points ) ){
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

    case 7:
/* no normal (s-) stretching */
      if (info->normal_stretch != NULL){
	info->normal_stretch = delete_generic_stretching( info->normal_stretch );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	printf("There is no stretching in the s-direction.\n");
	replot = 0;
      }
      break;

    case 8:
/* tangential (r-) stretching */
      if (info->tangential_stretch == NULL){
	info->tangential_stretch = choose_stretching( io_ptr, info->curve_ptr, 
						     &grid_ptr->r_points );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	if ( set_stretching( io_ptr, info->tangential_stretch, info->curve_ptr, 
			     &grid_ptr->r_points) ){
	  grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->tangential_stretch ) );
	  replot = 0;
	}
      }
      break;

    case 9:
/* no r-stretching */
      if (info->tangential_stretch != NULL){
	info->tangential_stretch = 
	  delete_generic_stretching( info->tangential_stretch );
	grid_ptr->grid_plot_mode = grid_ptr->grid_plot_mode | 7;
      }
      else{
	printf("There is no stretching in the r-direction.\n");
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
      set_grid_plot_mode( io_ptr, grid_ptr );
      replot = 0;
      break;

    case 13:
/* show parameters */
      show_curve_parameters( info->curve_ptr );
      if (info->is_constant)
	printf("Constant width: %e\n", info->width1);
      else{
	printf("Sharpness of width changes: %e\n", info->width_trans);
	printf("\n");
	printf("Number of corners: %i\n", info->n);
	printf("\n");
	printf("Widths:\n");
	if (info->flip){
	  for (i=1; i< info->n; i++)
	    printf("Interval %i, start width: %e, end width: %e\n", 
		   i, info->width_r[info->n-i], info->width_l[info->n-i]);
	}
	else{
	  for (i=1; i< info->n; i++)
	    printf("Interval %i, start width: %e, end width: %e\n", 
		   i, info->width_l[i], info->width_r[i]);
	}
	printf("\n");
      }
      show_mapping_parameters( grid_ptr );
      replot = 0;
      break;

    case 14:
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


static int 
normal_curve_mapping( grid_point *gp_ptr, void *data_ptr){
  normal_curve_data *info;
  curve_point cp;
  real sigma, sigma_r, w, w_r, r_stretch, dstr_duni, d2str_duni2, s_stretch, s_ratio
    , dummy;
/* cast the data to the right type */
  info = (normal_curve_data *) data_ptr;

/* flip the parameter? */
  if (info->flip)
    cp.r = 1.0 - gp_ptr->r;
  else
    cp.r = gp_ptr->r;

/* stretching in the tangential direction ? */
  uniform_to_stretch( cp.r, info->tangential_stretch, &r_stretch, &dstr_duni, 
		     &d2str_duni2 );

/* evaluate the width */
  if (info->is_constant){
    w = info->width1;
    w_r = 0;
  }
  else{
/* evaluate the width at the stretched coordinate */
    eval_width(info->width1, info->dndr, info->width_jmp, info->r, 
	       info->vn0, info->un0, info->n, info->width_trans, 
	       r_stretch, &w, &w_r);
/* scale the derivative to take the stretching into account */
    w_r = dstr_duni * w_r;
/* if the parametrization is reversed, change sign of the derivative */
    if (info->flip) w_r = - w_r;
  }

/* stretching in the normal direction? */
  uniform_to_stretch( gp_ptr->s, info->normal_stretch, &s_stretch, &s_ratio, 
		     &dummy );

  cp.r = r_stretch;
/* call the curve routine */
  curve_function( &cp, info->curve_ptr );
  if (info->flip){
    cp.xr = - cp.xr;
    cp.yr = - cp.yr;
  }

/* normalize the normal */
  sigma = sqrt( cp.xr*cp.xr + cp.yr*cp.yr);
/* the mapping will crash here if sigma is zero (should probbaly check for */
/* a small, but non-zero value instead */
  if (sigma == 0.0) return 1;

  sigma_r = dstr_duni * ( cp.xr * cp.xrr + cp.yr * cp.yrr )/ sigma;

/* Cartesian coordinates */
  gp_ptr->x = cp.x - cp.yr/sigma * s_stretch * w;
  gp_ptr->y = cp.y + cp.xr/sigma * s_stretch * w;

/* something is wrong in the jacobian */

/* Jacobian */
  gp_ptr->xr = cp.xr * dstr_duni - cp.yrr * dstr_duni / sigma * s_stretch * w 
      + cp.yr * sigma_r/(sigma*sigma) * s_stretch * w 
	- cp.yr/sigma * s_stretch * w_r;
  gp_ptr->xs = -cp.yr/sigma * w * s_ratio;
  gp_ptr->yr = cp.yr * dstr_duni + cp.xrr * dstr_duni / sigma * s_stretch * w
      - cp.xr * sigma_r/(sigma*sigma) * s_stretch * w 
	+ cp.xr/sigma * s_stretch * w_r;
  gp_ptr->ys = cp.xr/sigma * w * s_ratio;

  return 0;
}




static void 
compute_normal_curve(normal_curve_data *info, 
			  generic_mapping *grid_ptr ){
  int i,n;
  grid_point gp;

/* the grid is periodic if the curve is */
  grid_ptr->r_period = info->curve_ptr->periodic;
  grid_ptr->s_period = 0;

  if (!(info->is_constant)){
/* make sure info->r is a copy of that field in the smooth_polygon_curve structure */
    n = info->n;
/* assign the ghost cells */
    info->width_l[0] = info->width_r[0] = info->width_l[1];
    info->width_l[n] = info->width_r[n] = info->width_r[n-1];

/* compute slopes */
    for (i=1; i<n; i++){
      info->dndr[i] = (info->width_r[i] - info->width_l[i])/
	(info->r[i+1]-info->r[i]);
  }

/* compute jumps*/
    for (i=1; i<=n; i++){
      info->width_jmp[i] = info->width_l[i] - info->width_r[i-1];
    }

/* add padding nodes */
    info->r[0]    = info->r[1] - (info->r[2] - info->r[1]);
    info->r[n+1]  = info->r[n] + (info->r[n] - info->r[n-1]);
    info->dndr[0] = info->dndr[1];
    info->dndr[n] = info->dndr[n-1];

/* precompute vn_j(0) */
    for (i=0; i<=n; i++){
      info->vn0[i] = log_cosh(info->dndr[i], info->r[i], 
			      info->r[i+1], info->width_trans, 0.0 );
    }

/* precompute un_j(0) */
    for (i=1; i<=n; i++){
      info->un0[i] = 0.5 * info->width_jmp[i] * 
	tanh(-(info->width_trans)*info->r[i]);
    }

/* save width_l[1] in the structure */
    info->width1 = info->width_l[1];
  }

/* compute bounding box */
/* Important to use the local coordinates, i.e. to call normal_curve_mapping */
/* instead of forward_mapping, which adds rotation, translation and scaling */
  grid_ptr->x_min_0 = 1.e10;
  grid_ptr->x_max_0 = -1.e10;
  grid_ptr->y_min_0 = 1.e10;
  grid_ptr->y_max_0 = -1.e10;
/* s=0 */
  gp.s = 0.0;
  for (gp.r=0.0; gp.r<=1.05; gp.r+=0.1){
    normal_curve_mapping( &gp, grid_ptr->mapping_data_ptr );
    grid_ptr->x_min_0 = real_min( gp.x, grid_ptr->x_min_0 );
    grid_ptr->x_max_0 = real_max( gp.x, grid_ptr->x_max_0 );
    grid_ptr->y_min_0 = real_min( gp.y, grid_ptr->y_min_0 );
    grid_ptr->y_max_0 = real_max( gp.y, grid_ptr->y_max_0 );
  }
/* s=1 */
  gp.s = 1.0;
  for (gp.r=0.0; gp.r<=1.05; gp.r+=0.1){
    normal_curve_mapping( &gp, grid_ptr->mapping_data_ptr );
    grid_ptr->x_min_0 = real_min( gp.x, grid_ptr->x_min_0 );
    grid_ptr->x_max_0 = real_max( gp.x, grid_ptr->x_max_0 );
    grid_ptr->y_min_0 = real_min( gp.y, grid_ptr->y_min_0 );
    grid_ptr->y_max_0 = real_max( gp.y, grid_ptr->y_max_0 );
  }

  reset_view_port( grid_ptr );
}


static void *
cleanup_normal_curve( void *data_ptr ){
  normal_curve_data *info;

/* cast the data to the right type */
  info = (normal_curve_data *) data_ptr;

  if (info != NULL){
/* flag that the curve is no longer used */
    info->curve_ptr->used_by_grid--;
/* cleanup the width information */
    if (!(info->is_constant))
      free_width_arrays( info );

/* stretchings */
    delete_generic_stretching( info->tangential_stretch );
    delete_generic_stretching( info->normal_stretch );

/* cleanup the struct itself */
    free(info);
  }

  return NULL;
}

static void 
alloc_width_arrays( normal_curve_data *info ){
/* node points */
  info->r         = vector(0,info->n+1);
/* normal width */
  info->width_l   = vector(0,info->n);
  info->width_r   = vector(0,info->n);
  info->vn0       = vector(0,info->n);
  info->dndr      = vector(0,info->n);
  info->width_jmp = vector(1,info->n);
  info->un0       = vector(1,info->n);
}


static void 
free_width_arrays( normal_curve_data *info ){
  free_vector( info->r,         0, info->n+1 );
  free_vector( info->width_l,   0, info->n );
  free_vector( info->width_r,   0, info->n );
  free_vector( info->vn0,       0, info->n );
  free_vector( info->dndr,      0, info->n );
  free_vector( info->width_jmp, 1, info->n );
  free_vector( info->un0,       1, info->n );
}


static void 
eval_width(real w1, real *slope, real *width_jmp, real *r_node, 
		       real *v0, real *u0, int n, real width_trans, 
		       real r, real *w, real *w_r){
  int i;

/* function */
  *w = w1;
/* assume normalized input in the range 0 <= r<= 1.0 */
  r = r*r_node[n];
/* ramps */
  for (i=0; i<=n; i++)
    *w += log_cosh(slope[i], r_node[i], r_node[i+1], width_trans, r) - v0[i];
/* jumps */
  for (i=1; i<=n; i++)
    *w += 0.5 * width_jmp[i] * tanh( width_trans*(r - r_node[i]) ) - u0[i];

/* derivative */
  *w_r = 0.0;
/* ramps */
  for (i=0; i<=n; i++)
    *w_r +=  0.5 * slope[i] * (tanh(width_trans*(r-r_node[i])) - 
			    tanh(width_trans*(r-r_node[i+1])) );
/* jumps */
  for (i=1; i<=n; i++)
    *w_r += 0.5 * width_jmp[i] * width_trans * 
      1.0/cosh(width_trans*(r - r_node[i]))/cosh(width_trans*(r - r_node[i]));
  *w_r = *w_r * r_node[n];
}


static void 
write_normal_curve( int32 dir, void *normal_curve_data_ptr ){
  normal_curve_data *info;
  int32 stretch_dir;

  info = (normal_curve_data *) normal_curve_data_ptr;

/* save the name of the curve */
  hput_string( info->curve_ptr->curve_name, "curve_name", dir );

  hput_int( info->is_constant, "is_constant", dir );
  hput_int( info->flip, "flip", dir );
  hput_int( info->n, "n", dir );
  hput_real(info->width_trans, "width_trans", dir );
  hput_real(info->width1, "width1", dir );

/* variable width? */
  if (!info->is_constant){
    hput_vector( info->r,         0, info->n+1, "r",         dir );
    hput_vector( info->width_l,   0, info->n,   "width_l",   dir );
    hput_vector( info->width_r,   0, info->n,   "width_r",   dir );
    hput_vector( info->vn0,       0, info->n,   "vn0",       dir );
    hput_vector( info->dndr,      0, info->n,   "dndr",      dir );
    hput_vector( info->width_jmp, 1, info->n,   "width_jmp", dir );
    hput_vector( info->un0,       1, info->n,   "un0",       dir );
  }

/* tangential stretching */
  if ((stretch_dir = create_dir("tangential-stretching", 
				"tangential-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( stretch_dir, info->tangential_stretch );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }
  else
    printf("Error: save_curves_mappings: unable to create a directory for "
	   "`tangential-stretching'.\n");

/* normal stretching */
  if ((stretch_dir = create_dir("normal-stretching", 
				"normal-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( stretch_dir, info->normal_stretch );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }
  else
    printf("Error: save_curves_mappings: unable to create a directory for "
	   "`normal-stretching'.\n");

}


normal_curve_data *
init_normal_curve(input_output *io_ptr, generic_mapping *grid_ptr,
		  generic_curve *new_curve_ptr, char *status){
  normal_curve_data *info;
  
/* allocate memory */
  if ((info = (normal_curve_data *) malloc( sizeof(normal_curve_data) )) == NULL)
    printf("memory error in init_normal_curve\n");

/* point to the new curve */
  info->curve_ptr = new_curve_ptr;
/* flag that the new curve is used */
  info->curve_ptr->used_by_grid++; 

/* constant width 0.5 */
  info->is_constant = 1;
  info->width_trans = 5.0;
  info->width1 = 0.5;
  info->flip = 0;

/* no variable with arrays when the width is constant */
  info->n = 0;

/* initially no curve stretching */
  info->tangential_stretch = NULL;
  info->normal_stretch = NULL;

/* set the grid type and the periodicity flag */
  grid_ptr->grid_type= 4; /* normal-curve grid */

/* set the number of grid lines */
  grid_ptr->r_points = 30;
  grid_ptr->s_points = 7;

/* save the pointers to the data and the functions for
   evaluating the mapping */
  grid_ptr->grid_mapping = normal_curve_mapping;
  grid_ptr->inverse_known = 0;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->mapping_data_ptr = info;
  grid_ptr->cleanup_mapping_data = cleanup_normal_curve;
  grid_ptr->set_mapping = set_normal_curve_mapping;
  grid_ptr->write_specific_mapping = write_normal_curve;
  grid_ptr->copy_mapping_data = copy_normal_curve;

  compute_normal_curve( info, grid_ptr );

  return info;
}


int 
read_normal_curve( int32 dir, generic_mapping *grid_ptr, generic_curve *first_c ){
  normal_curve_data *info;
  char *curve_name;
  int low, high;
  int32 stretch_dir;

  info = (normal_curve_data *) malloc( sizeof(normal_curve_data) );

/* read the name */
  curve_name = hget_string( "curve_name", dir );

/* look for this curve in the curve list */
  for (info->curve_ptr = first_c; 
       info->curve_ptr != NULL && 
       strcmp(info->curve_ptr->curve_name, curve_name) != 0;
       info->curve_ptr = info->curve_ptr->next);

/* was the curve found? */
  if (info->curve_ptr == NULL){
    printf("ERROR: read_normal_curve. Could not find curve with name: %s\n",
	   curve_name);
    free( curve_name );
    free( info );
    return FALSE;
  }
  free( curve_name );

  hget_int(  &(info->is_constant), "is_constant", dir );
  hget_int(  &(info->flip),        "flip",        dir );
  hget_int(  &(info->n),           "n",           dir );
  hget_real( &(info->width_trans), "width_trans", dir );
  hget_real( &(info->width1),      "width1",      dir );

/* variable width ? */
  if (!info->is_constant){
/* read the vectors */
    info->r         = hget_vector( "r",         &low, &high, dir );
    info->width_l   = hget_vector( "width_l",   &low, &high, dir );
    info->width_r   = hget_vector( "width_r",   &low, &high, dir );
    info->vn0       = hget_vector( "vn0",       &low, &high, dir );
    info->dndr      = hget_vector( "dndr",      &low, &high, dir );
    info->width_jmp = hget_vector( "width_jmp", &low, &high, dir );
    info->un0       = hget_vector( "un0",       &low, &high, dir );
  }
  else{
    info->r         = NULL;
    info->width_l   = NULL;
    info->width_r   = NULL;
    info->vn0       = NULL;
    info->dndr      = NULL;
    info->width_jmp = NULL;
    info->un0       = NULL;
  }

/* read the stretchings */
  if ((stretch_dir = locate_dir("tangential-stretching", dir)) != -1){
/* save the generic and specific curve data */
    info->tangential_stretch = read_stretching( stretch_dir );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }

  if ((stretch_dir = locate_dir("normal-stretching", dir)) != -1){
/* save the generic and specific curve data */
    info->normal_stretch     = read_stretching( stretch_dir );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }

/* setup the pointers to all specific functions */
  grid_ptr->mapping_data_ptr = (void *) info;
  grid_ptr->grid_mapping = normal_curve_mapping;
  grid_ptr->inverse_known = 0;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->cleanup_mapping_data = cleanup_normal_curve;
  grid_ptr->set_mapping = set_normal_curve_mapping;
  grid_ptr->write_specific_mapping = write_normal_curve;
  grid_ptr->copy_mapping_data = copy_normal_curve;

  return OK;
}

static void *
copy_normal_curve( void *mapping_data_ptr ){
  normal_curve_data *info, *old_info;
  int i;
  
  info = (normal_curve_data *) malloc( sizeof(normal_curve_data) );

  old_info = (normal_curve_data *)mapping_data_ptr;

/* only make a shallow copy of the curve. */
  info->curve_ptr = old_info->curve_ptr;

/* remember that this curve is now used by both mappings */
  info->curve_ptr->used_by_grid++;

  info->is_constant = old_info->is_constant;
  info->flip = old_info->flip;
  info->n = old_info->n;
  info->width_trans = old_info->width_trans;
  info->width1 = old_info->width1;

/* make a copy of the stretching functions */
  info->tangential_stretch = copy_stretching( old_info->tangential_stretch );
  info->normal_stretch = copy_stretching( old_info->normal_stretch );

/* does the width vary? */
  if (!(info->is_constant)){
/* allocate the width arrays */
    alloc_width_arrays( info );
/* copy the values */
    for (i=1; i<= info->n; i++){
      info->r[i]       = old_info->r[i];
      info->width_l[i] = old_info->width_l[i];
      info->width_r[i] = old_info->width_r[i];
      info->vn0[i]     = old_info->vn0[i];
      info->dndr[i]    = old_info->dndr[i];
      info->width_jmp[i] = old_info->width_jmp[i];
      info->un0[i]     = old_info->un0[i];
    }
/* first and last elements */
    info->r[0] = old_info->r[0];  info->r[info->n+1] = old_info->r[info->n+1];
    info->width_l[0] = old_info->width_l[0];
    info->width_r[0] = old_info->width_r[0];
    info->vn0[0]     = old_info->vn0[0];
    info->dndr[0]    = old_info->dndr[0];
  }

  return (void *)info;
}
