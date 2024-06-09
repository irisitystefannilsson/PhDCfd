#include "curves.h"

/* private member functions */
static void 
alloc_spline_arrays( cubic_spline_data *info );
static void 
free_spline_arrays( cubic_spline_data *info );
static void 
compute_cubic_spline(cubic_spline_data *info, generic_curve *curve_ptr );
static int 
cubic_spline_curve( curve_point *cp_ptr, void *data_ptr );
static void *
cleanup_cubic_spline( void *data_ptr );
static void 
write_cubic_spline( int32 dir, void *cubic_spline_data_ptr );
static void 
plot_cubic_spline( generic_curve *curve_ptr );
static void 
set_cubic_spline( input_output *io_ptr, generic_curve *curve_ptr);
static real *
copy_spline_nodes( generic_curve *curve_ptr, int *n );
static void *
copy_cubic_spline( void *curve_data_ptr );
/* end private member functions */

static void 
set_cubic_spline( input_output *io_ptr, generic_curve *curve_ptr){
  const real eps=1.e-7;
  real new_x, new_y, xi, yi;
  char prompt[80], question[80];
  int icom, quit, i, replot;
  int level=2, no_indent=0;
  cubic_spline_data *info;
  curve_point cp;

#include "cubic_spline_com.h"

  sprintf(prompt, "%s: cubic-spline curve>", curve_ptr->curve_name);

  info = (cubic_spline_data *) curve_ptr->curve_data_ptr;

  quit = 0;
  PL_window(2, curve_ptr->xyab);
  replot = 1;
  do{

    if (replot){
/* plot the curve */
      PL_erase(2);
      PL_start_plot(2);
      plot_curve( curve_ptr, curve_ptr->plot_mode );
      PL_end_plot();
    }

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {
    case 0:
/* enter nodes */ 
      init_cubic_spline( io_ptr, curve_ptr, "new" );
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      PL_window(2, curve_ptr->xyab);
      break;

    case 1:
/* read nodes from file */
      if ((init_cubic_spline( io_ptr, curve_ptr, "file" )) != NULL){
	curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
	PL_window(2, curve_ptr->xyab);
      }
      else{
	printf("No new spline was created because of the above error(s).\n");
	replot = 0;
      }
      break;

    case 2:
/* change node */
      i = get_int(io_ptr , "Which node: ", 1, no_indent);
      if (i < 1 || i > info->n)
	printf("That node does not exist.\n");
      else{
	sprintf(question, "node %i: x: ", i);
	rot_trans_scale( info->x[i], info->y[i], &xi, &yi, curve_ptr );
	new_x = get_real( io_ptr, question, xi, level+1);
	sprintf(question, "node %i: y: ", i);
	new_y = get_real( io_ptr, question, yi, level+1);
	inv_rot_trans_scale( new_x, new_y, &new_x, &new_y, curve_ptr );
	while ((i > 1 && real_max(fabs(info->x[i-1] - new_x), 
				  fabs(info->y[i-1] - new_y) ) < eps ) ||
	       (i < info->n && real_max(fabs(info->x[i+1] - new_x), 
					fabs(info->y[i+1] - new_y) ) < eps )){
	  printf("Input Error: consecetive node points may not be identical!\n");
	  sprintf(question, "node %i: x: ", i);
	  new_x = get_real( io_ptr, question, xi, level+1);
	  sprintf(question, "node %i: y: ", i);
	  new_y = get_real( io_ptr, question, yi, level+1);
	  inv_rot_trans_scale( new_x, new_y, &new_x, &new_y, curve_ptr );
	}
	info->x[i] = new_x;
	info->y[i] = new_y;
	compute_cubic_spline( info, curve_ptr );
	PL_window(2, curve_ptr->xyab);
	curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      }
      break;

    case 3:
      info->bc_x1 = get_real(io_ptr, "Enter dx/ds(0) (>=1e30 means free): ",
                             info->bc_x1, no_indent);
      compute_cubic_spline( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 4:
      info->bc_y1 = get_real(io_ptr, "Enter dy/ds(0) (>=1e30 means free): ",
                             info->bc_y1, no_indent);
      compute_cubic_spline( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 5:
      info->bc_xn = get_real(io_ptr, "Enter dx/ds(smax) (>=1e30 means free): ",
                             info->bc_xn, no_indent);
      compute_cubic_spline( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 6:
      info->bc_yn = get_real(io_ptr, "Enter dy/ds(smax) (>=1e30 means free): ",
                             info->bc_yn, no_indent);
      compute_cubic_spline( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 7:
/* plot mode */
      set_curve_plot_mode( io_ptr, curve_ptr );
      replot = 0;
      break;

    case 8:
/* show parameters */
      printf("Number of nodes: %i\n", info->n);
      printf("\n");
      printf("Node coordinates:\n");
      for (i=1; i<= info->n; i++){
	rot_trans_scale( info->x[i], info->y[i], &xi, &yi, curve_ptr );
	printf("Node %i, x: %e, y: %e\n", i, xi, yi);
      }
      printf("\n");
/* evaluate the derivatives */
      printf("Boundary conditions:\n");

/* r=0 */
      cp.r=0.0;
      cubic_spline_curve( &cp, info );

/* x */
      if (info->bc_x1 > 0.99e30)
	printf("dx/dt(start): %e (free)\n", cp.xr/info->rmax);
      else
	printf("dx/dt(start): %e (clamped)\n", info->bc_x1);
/* y */
      if (info->bc_y1 > 0.99e30)
	printf("dy/dt(start): %e (free)\n", cp.yr/info->rmax);
      else
	printf("dy/dt(start): %e (clamped)\n", info->bc_y1);

/* r=1 */
      cp.r=1.0;
      cubic_spline_curve( &cp, info );

/* x */
      if (info->bc_xn > 0.99e30)
	printf("dx/dt(end):  %e (free)\n", cp.xr/info->rmax);
      else
	printf("dx/dt(end):  %e (clamped)\n", info->bc_xn);
/* y */
      if (info->bc_yn > 0.99e30)
	printf("dy/dt(end):  %e (free)\n", cp.yr/info->rmax);
      else
	printf("dy/dt(end):  %e (clamped)\n", info->bc_yn);

      replot = 0;
      break;

    case 9:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case 10:
/* quit */
      quit = 1;
      break;
    default:
      replot = 0;
    }
  }
  while(!quit);
}

void *init_cubic_spline( input_output *io_ptr, generic_curve *curve_ptr,
			     char *status ){
  cubic_spline_data *info=NULL;
  
  int i, n, level=3, no_indent=0, new_n, read_error;
  const int save_on_copy=1;
  real prev_x, prev_y;
  const real eps=1.e-7;
  char question[80], *file_name;
  FILE *fp;
  real *new_x, *new_y;
  
/* initialize a dummy grid */

  if (strcmp("init", status)==0){
    if ((info = (cubic_spline_data *) malloc( sizeof(cubic_spline_data ) )) == NULL)
      printf("memory error in init_cubic_spline\n");
    info->n = n = 2;
/* allocate space */
    alloc_spline_arrays( info );
/* a trivial but valid curve */
    info->x[1]=0.0; info->y[1]=0.0;
    info->x[2]=1.0; info->y[2]=0.0;
    info->bc_x1 = FREE_BC;
    info->bc_xn = FREE_BC;
    info->bc_y1 = FREE_BC;
    info->bc_yn = FREE_BC;
  }

/* assign all nodes */

  else if (strcmp("new", status)==0){
/* check that at least some things are ok */
    if (curve_ptr->curve_mapping != cubic_spline_curve ||
	curve_ptr->curve_data_ptr==NULL){
      printf("init_cubic_spline: status is `new' but the pointers are incorrect\n");
      exit(1);
    }
/* look up the existing structure */
    info = (cubic_spline_data *) curve_ptr->curve_data_ptr;
/* remove old nodes */
    free_spline_arrays( info );
/* read number of node points */
    info->n = n = 
      int_max(2, get_int(io_ptr , "Enter number of nodes >= 2: ", 2,
			 no_indent));
/* allocate space */
    alloc_spline_arrays( info );
/* read the nodes */
    prev_x = -1.0; prev_y = -1.0;
    for (i=1; i<=n; i++){
      sprintf(question, "node %i: x: ", i);
      info->x[i] = get_real( io_ptr, question, prev_x + 1.0, level);
      sprintf(question, "node %i: y: ", i);
      info->y[i] = get_real( io_ptr, question, prev_y + 1.0, level);
/* take the global transformation into account */
      inv_rot_trans_scale( info->x[i], info->y[i], 
			  &(info->x[i]), &(info->y[i]), curve_ptr );
/* check for uniqueness */
      while (i > 1 && real_max(fabs(info->x[i] - prev_x), 
			       fabs(info->y[i] - prev_y) ) < eps ){
	printf("Input Error: consecetive node points may not be identical!\n");
	sprintf(question, "node %i: x: ", i);
	info->x[i] = get_real( io_ptr, question, prev_x + 1.0, level);
	sprintf(question, "node %i: y: ", i);
	info->y[i] = get_real( io_ptr, question, prev_y + 1.0, level);
/* take the global transformation into account */
	inv_rot_trans_scale( info->x[i], info->y[i], 
			    &(info->x[i]), &(info->y[i]), curve_ptr );
      }
      prev_x = info->x[i];
      prev_y = info->y[i];
/* take the global transformation into account */
      rot_trans_scale( prev_x, prev_y, &prev_x, &prev_y, curve_ptr );
    }
  }
  else if (strcmp("file", status)==0){
/* check that at least some things are ok */
    if (curve_ptr->curve_mapping != cubic_spline_curve ||
	curve_ptr->curve_data_ptr==NULL){
      printf("init_cubic_spline: status is `file' but the pointers are incorrect\n");
      exit(1);
    }
/* open the file */
    if ((fp = open_ascii_file( io_ptr, "Enter node-point file: ", "test.spl", 
			      &file_name, 'r', 4, save_on_copy)) == NULL) return NULL;
/* read number of node points */
    if (fscanf( fp, "%i", &new_n ) == EOF){
      printf("Error: End-of-file occured while reading spline node-points\n");
      fclose(fp);
      return NULL;
    }
    if ((n = new_n) < 2){
      printf("Error: The number of node-points %i < 2\n", n);
      return NULL;
    }

/* allocate space */
    new_x = vector(1, new_n); 
    new_y = vector(1, new_n);

/* read the nodes */
    read_error = 0;
    for (i=1; i<=n && !read_error; i++){
#ifdef SINGLE
      if (fscanf( fp, "%g%g", &(new_x[i]), &(new_y[i]) ) == EOF){
#else
      if (fscanf( fp, "%lg%lg", &(new_x[i]), &(new_y[i]) ) == EOF){
#endif
	printf("Error: End-of-file occured while reading node-point %i\n",i);
	read_error = 1;
      }
    }
/* close the file */
    fclose( fp );
/* check if all consequtive spline points are unique */
    if ( !read_error ){
      for (i=2; i<=n && !read_error; i++)
	if (real_max(fabs(new_x[i] - new_x[i-1]), 
		     fabs(new_y[i] - new_y[i-1]) ) < eps){
	  printf("Error: Nodes %i and %i are identical\n", i-1, i);
	  read_error = 1;
	}
    }
/* everything seems ok! */
    if ( !read_error ){
/* look up the existing structure */
      info = (cubic_spline_data *) curve_ptr->curve_data_ptr;
/* it is not necessary to remove the space if the number of nodes is unchanged */
      if (new_n != info->n){
/* remove old nodes */
	free_spline_arrays( info );
/* allocate new arrays */
	info->n = new_n;
	alloc_spline_arrays( info );
      }
/* copy the new arrays and take the global transformation into account */
      for (i=1; i<=n; i++){
/* take the global transformation into account */
	inv_rot_trans_scale( new_x[i], new_y[i], 
			    &(info->x[i]), &(info->y[i]), curve_ptr );
      }
    }
/* free the temporary storage */
    free_vector( new_x, 1, new_n );
    free_vector( new_y, 1, new_n );

/* tell the calling routine if there was an error */
    if (read_error) return NULL;
  }
  else{
    printf("init_cubic_spline: Unknown status value.\n");
    exit(1);
  }

  curve_ptr->curve_type = 1;
  curve_ptr->curve_data_ptr = (void *) info;
  curve_ptr->curve_mapping = cubic_spline_curve;
  curve_ptr->cleanup_curve_data = cleanup_cubic_spline;
  curve_ptr->set_curve = set_cubic_spline;
  curve_ptr->write_specific_curve = write_cubic_spline;
  curve_ptr->plot_control_points = plot_cubic_spline;
  curve_ptr->copy_nodes = copy_spline_nodes;
  curve_ptr->copy_curve_data = copy_cubic_spline;

  compute_cubic_spline( info, curve_ptr );

  return info;

}

static void 
compute_cubic_spline(cubic_spline_data *info, generic_curve *curve_ptr ){
  int i,n;
  real box_size;
/*  arcl_stretch_info *stretch_ptr;*/

  n = info->n;
/* update the number of segments */
  curve_ptr->r_step = 0.1/n;

/* compute the pseudo arclength */
  info->r[1] = 0.0;
  for (i=2; i<=n; i++){
    info->r[i] = info->r[i-1] + 
      sqrt((info->x[i]-info->x[i-1])*(info->x[i]-info->x[i-1]) +
	   (info->y[i]-info->y[i-1])*(info->y[i]-info->y[i-1]));
  }
  info->rmax = info->r[n];

/* compute spline coefficients */
  spline(info->r, info->x, info->n, info->bc_x1, info->bc_xn, info->xcoeff);
  spline(info->r, info->y, info->n, info->bc_y1, info->bc_yn, info->ycoeff);

/* compute bounding box */
  curve_ptr->x_min_0 = 1.e10;
  curve_ptr->x_max_0 = -1.e10;
  curve_ptr->y_min_0 = 1.e10;
  curve_ptr->y_max_0 = -1.e10;
  for (i=1; i<=n; i++){
    curve_ptr->x_min_0 = real_min( info->x[i], curve_ptr->x_min_0 );
    curve_ptr->x_max_0 = real_max( info->x[i], curve_ptr->x_max_0 );
    curve_ptr->y_min_0 = real_min( info->y[i], curve_ptr->y_min_0 );
    curve_ptr->y_max_0 = real_max( info->y[i], curve_ptr->y_max_0 );
  }
/* estimate lengthscale of spline */
  box_size = real_max(curve_ptr->y_max_0 - curve_ptr->y_min_0, 
		      curve_ptr->x_max_0 - curve_ptr->x_min_0);


/* check the periodicity */
  if (fabs(info->x[n] - info->x[1]) + fabs(info->y[n] - info->y[1]) <=
      1.e-7 * box_size)
    curve_ptr->periodic = 1;
  else
    curve_ptr->periodic = 0;
  info->periodic = curve_ptr->periodic;

/* set view port */
  curve_view_port( curve_ptr );
}

static int 
cubic_spline_curve( curve_point *cp_ptr, void *data_ptr ){
  cubic_spline_data *info;
  real r, rmax, xr, xrr, yr, yrr, r_stretch, dstr_duni, d2str_duni2;

/* cast the data to the right type */
  info = (cubic_spline_data *) data_ptr;

  r_stretch = cp_ptr->r;

  if (info->periodic){
    if (r_stretch < 0.0)
      r_stretch += 1.0;
    else if (r_stretch > 1.0)
      r_stretch -= 1.0;
  }

  dstr_duni = 1.0;
  d2str_duni2 = 0.0;

/* scale the parameter */
  rmax = info->r[info->n]; 
  r = rmax * r_stretch;

  cseval( info->r, info->x, info->xcoeff, info->n, r, 
	 &(cp_ptr->x), &xr, &xrr );

  cseval( info->r, info->y, info->ycoeff, info->n, r,
	 &(cp_ptr->y), &yr, &yrr );

/* scale the derivatives with respect to the scaling by rmax */
  cp_ptr->xr  =  xr * rmax; 
  cp_ptr->yr  =  yr * rmax;
  cp_ptr->xrr = xrr * rmax * rmax; 
  cp_ptr->yrr = yrr * rmax * rmax;

  return 0;
}

static void *
cleanup_cubic_spline( void *data_ptr ){
  cubic_spline_data *info;

/* cast the data to the right type */
  info = (cubic_spline_data *) data_ptr;

  if (info != NULL){
    free_spline_arrays( info );
    free(info);
  }

  return NULL;
}

static void 
alloc_spline_arrays( cubic_spline_data *info ){
  info->r      = vector(1,info->n);
  info->x      = vector(1,info->n); 
  info->y      = vector(1,info->n);
  info->xcoeff = vector(1,info->n); 
  info->ycoeff = vector(1,info->n);
}


static void 
free_spline_arrays( cubic_spline_data *info ){
  if (info != NULL){
    free_vector( info->r     , 1, info->n);
    free_vector( info->x     , 1, info->n); 
    free_vector( info->y     , 1, info->n);
    free_vector( info->xcoeff, 1, info->n); 
    free_vector( info->ycoeff, 1, info->n);
  }
}

void 
read_cubic_spline( int32 dir, generic_curve *curve_ptr ){
  cubic_spline_data *info;
  int low, high;

  info = (cubic_spline_data *) malloc( sizeof(cubic_spline_data) );
  curve_ptr->curve_data_ptr = (void *) info;
  curve_ptr->curve_mapping = cubic_spline_curve;
  curve_ptr->cleanup_curve_data = cleanup_cubic_spline;
  curve_ptr->set_curve = set_cubic_spline;
  curve_ptr->write_specific_curve = write_cubic_spline;
  curve_ptr->plot_control_points = plot_cubic_spline;
  curve_ptr->copy_nodes = copy_spline_nodes;
  curve_ptr->copy_curve_data = copy_cubic_spline;

/* read the vector size and the periodicity */
  hget_int(&(info->n), "n", dir );
  hget_int(&(info->periodic), "periodic", dir );

/* read the vectors */
/* low should be 1, and high should be info->n */
  info->r = hget_vector( "r", &low, &high, dir ); 
  info->x = hget_vector( "x", &low, &high, dir );
  info->y = hget_vector( "y", &low, &high, dir );
  info->xcoeff = hget_vector( "xcoeff", &low, &high, dir );
  info->ycoeff = hget_vector( "ycoeff", &low, &high, dir );

/* boundary conditions */
  hget_real( &(info->bc_x1), "bc_x1", dir );
  hget_real( &(info->bc_xn), "bc_xn", dir );
  hget_real( &(info->bc_y1), "bc_y1", dir );
  hget_real( &(info->bc_yn), "bc_yn", dir );
  hget_real( &(info->rmax),  "rmax", dir );

}

static void 
write_cubic_spline( int32 dir, void *cubic_spline_data_ptr ){
  cubic_spline_data *info;

  info = (cubic_spline_data *) cubic_spline_data_ptr;

/* save the size and the periodicity */
  hput_int( info->n, "n", dir );
  hput_int( info->periodic, "periodic", dir );

/* write the vectors */
  hput_vector( info->r, 1, info->n, "r", dir );
  hput_vector( info->x, 1, info->n, "x", dir );
  hput_vector( info->y, 1, info->n, "y", dir );
  hput_vector( info->xcoeff, 1, info->n, "xcoeff", dir );
  hput_vector( info->ycoeff, 1, info->n, "ycoeff", dir );

  hput_real( info->bc_x1, "bc_x1", dir);
  hput_real( info->bc_xn, "bc_xn", dir);
  hput_real( info->bc_y1, "bc_y1", dir);
  hput_real( info->bc_yn, "bc_yn", dir);
  hput_real( info->rmax, "rmax", dir);

}

static void 
plot_cubic_spline( generic_curve *curve_ptr ){
  cubic_spline_data *info;
  int i;
  char cflag[10];
  real x_new, y_new;

  info = (cubic_spline_data *) curve_ptr->curve_data_ptr;

  for (i=1; i<= info->n; i++){
    rot_trans_scale( info->x[i], info->y[i], &x_new, &y_new, curve_ptr );
    PL_marker( x_new, y_new, ASTERISK);
    sprintf(cflag, "%i", i);
    PL_plot_string(cflag, x_new, y_new, 1, 1, 3);
  }
}


static real *
copy_spline_nodes( generic_curve *curve_ptr, int *n ){
  cubic_spline_data *info;
  real *nodes;
  int i;

  info = (cubic_spline_data *) curve_ptr->curve_data_ptr;

/* allocate a node vector */
  nodes = vector( 1, info->n );

/* copy node-points from the geometry description */
  for (i=1; i<=info->n; i++){ 
    nodes[i] = info->r[i]/info->r[info->n]; 
  } 

  *n = info->n;
  return nodes;
}


static void *
copy_cubic_spline( void *curve_data_ptr ){
  cubic_spline_data *info, *old_info;
  int i;

  info = (cubic_spline_data *) malloc( sizeof(cubic_spline_data ) );

  old_info = (cubic_spline_data *) curve_data_ptr;

  info->n     = old_info->n;
  info->bc_x1 = old_info->bc_x1;
  info->bc_xn = old_info->bc_xn;
  info->bc_y1 = old_info->bc_y1;
  info->bc_yn = old_info->bc_yn;
  info->rmax  = old_info->rmax;

/* allocate the arrays */
  alloc_spline_arrays( info );

/* copy the array elements */
  for (i=1; i<=info->n; i++){
    info->r[i]      = old_info->r[i];
    info->x[i]      = old_info->x[i];
    info->y[i]      = old_info->y[i];
    info->xcoeff[i] = old_info->xcoeff[i];
    info->ycoeff[i] = old_info->ycoeff[i];
  }

  return info;
}
