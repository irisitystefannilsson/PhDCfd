#include "curves.h"
#include "log_cosh.h"

/* private member functions */
static void alloc_poly_curve_arrays( smooth_poly_data *info );
static void free_poly_curve_arrays( smooth_poly_data *info );
static void eval_polygon(real x1, real *slope, real *r_node, real *v0, 
			 int n, real sharpness, real r, real *x, real *x_r, 
			 real *x_rr);
static int smooth_poly_curve( curve_point *cp_ptr, void *smooth_poly_data_ptr );
static void *cleanup_smooth_poly_curve( void *smooth_poly_data_ptr );
static void compute_smooth_poly_curve(smooth_poly_data *info, 
				      generic_curve *curve_ptr );
static void 
write_smooth_poly( int32 dir, void *smooth_poly_data_ptr );
static void plot_smooth_poly( generic_curve *curve_ptr );
static void set_smooth_polygon(input_output *io_ptr, generic_curve *curve_ptr);
static real *copy_polygon_nodes( generic_curve *curve_ptr, int *n );
static void *copy_smooth_poly( void *curve_data_ptr );
/* end private member functions */

static void set_smooth_polygon( input_output *io_ptr, generic_curve *curve_ptr){
  char prompt[80], question[120];
  const real eps=1.e-7;
  real new_x, new_y, xi, yi;
  int icom, quit, i, replot;
  int level=2, no_indent=0;
  smooth_poly_data *info;

#include "smooth_poly_curve_com.h"

  sprintf(prompt, "%s: smooth-polygon curve>", curve_ptr->curve_name);

  info = (smooth_poly_data *) curve_ptr->curve_data_ptr;

  quit = 0;
  replot = 1;
  PL_window(2, curve_ptr->xyab);
  do{

    if (replot){
/* plot the mapping */
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
/* enter geometry (corners and widths) */ 
      info = init_smooth_poly_curve( io_ptr, curve_ptr, "new" );
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      PL_window(2, curve_ptr->xyab);
      break;

    case 1:
/* change corner */
      i = get_int(io_ptr , "Which corner: ", 1, no_indent);
      if (i < 1 || i > info->n)
	printf("That corner does not exist.\n");
      else{
	sprintf(question, "corner %i: x: ", i);
	rot_trans_scale( info->x[i], info->y[i], &xi, &yi, curve_ptr );
	new_x = get_real( io_ptr, question, xi, level+1);
	sprintf(question, "corner %i: y: ", i);
	new_y = get_real( io_ptr, question, yi, level+1);
	inv_rot_trans_scale( new_x, new_y, &new_x, &new_y, curve_ptr );
	while ((i > 1 && real_max(fabs(info->x[i-1] - new_x), 
				  fabs(info->y[i-1] - new_y) ) < eps ) ||
	       (i < info->n && real_max(fabs(info->x[i+1] - new_x), 
					fabs(info->y[i+1] - new_y) ) < eps )){
	  printf("Input Error: consecutive corner points may not be identical!\n");
	  sprintf(question, "corner %i: x: ", i);
	  new_x = get_real( io_ptr, question, xi, level+1);
	  sprintf(question, "corner %i: y: ", i);
	  new_y = get_real( io_ptr, question, yi, level+1);
	  inv_rot_trans_scale( new_x, new_y, &new_x, &new_y, curve_ptr );
	}
	info->x[i] = new_x;
	info->y[i] = new_y;
	compute_smooth_poly_curve( info, curve_ptr );
	PL_window(2, curve_ptr->xyab);
	curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      }
      break;

    case 2:
/* corner sharpness */
      info->sharpness = 
	real_max(eps,
		 get_real(io_ptr, "Sharpness of corners > 0: ", info->sharpness,
			  no_indent));
      compute_smooth_poly_curve( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 3:
/* plot mode */
      set_curve_plot_mode( io_ptr, curve_ptr );
      replot = 0;
      break;

    case 4:
/* show parameters */
      printf("Sharpness of corners: %e\n", info->sharpness);
      printf("\n");
      printf("Number of corners: %i\n", info->n);
      printf("\n");
      printf("Corner coordinates:\n");
      for (i=1; i<= info->n; i++){
	rot_trans_scale( info->x[i], info->y[i], &xi, &yi, curve_ptr );
	printf("Corner %i, x: %e, y: %e\n", i, xi, yi);
      }
      printf("\n");
      replot = 0;
      break;

    case 5:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case 6:
/* quit */
      quit = 1;
      break;
    default:
      replot = 0;
    }
  }
  while(!quit);
}

static void compute_smooth_poly_curve(smooth_poly_data *info, 
				      generic_curve *curve_ptr ){
  int i,n;
  curve_point cp;

  n = info->n;

/* compute the pseudo arclength */
  info->r[1] = 0.0;
  for (i=2; i<=n; i++){
    info->r[i] = info->r[i-1] + 
      sqrt((info->x[i]-info->x[i-1])*(info->x[i]-info->x[i-1]) +
	   (info->y[i]-info->y[i-1])*(info->y[i]-info->y[i-1]));
  }

/* compute slopes */
  for (i=1; i<n; i++){
    info->dxdr[i] = (info->x[i+1]-info->x[i])/(info->r[i+1]-info->r[i]);
    info->dydr[i] = (info->y[i+1]-info->y[i])/(info->r[i+1]-info->r[i]);
  }

/* add padding nodes */
  info->r[0] = info->r[1] - (info->r[2] - info->r[1]);
  info->r[n+1] = info->r[n] + (info->r[n] - info->r[n-1]);
  info->dxdr[0] = info->dxdr[1];
  info->dxdr[n] = info->dxdr[n-1];
  info->dydr[0] = info->dydr[1];
  info->dydr[n] = info->dydr[n-1];

/* precompute vx_j(0), vy_j(0), vn_j(0) */
  for (i=0; i<=n; i++){
    info->vx0[i] = log_cosh(info->dxdr[i], info->r[i], 
			    info->r[i+1], info->sharpness, 0.0 );
    info->vy0[i] = log_cosh(info->dydr[i], info->r[i], 
			    info->r[i+1], info->sharpness, 0.0 );
  }

/* save x[1] and y[1] in the structure */
  info->x1 = info->x[1];
  info->y1 = info->y[1];

/* compute bounding box */
  curve_ptr->x_min_0 = 1.e10;
  curve_ptr->x_max_0 = -1.e10;
  curve_ptr->y_min_0 = 1.e10;
  curve_ptr->y_max_0 = -1.e10;

  for (cp.r=0.0; cp.r<=1.0 + 0.5*curve_ptr->r_step; cp.r+=curve_ptr->r_step){
    smooth_poly_curve( &cp, curve_ptr->curve_data_ptr );
    curve_ptr->x_min_0 = real_min( cp.x, curve_ptr->x_min_0 );
    curve_ptr->x_max_0 = real_max( cp.x, curve_ptr->x_max_0 );
    curve_ptr->y_min_0 = real_min( cp.y, curve_ptr->y_min_0 );
    curve_ptr->y_max_0 = real_max( cp.y, curve_ptr->y_max_0 );
  }

/* set view port */
  curve_view_port( curve_ptr );
}

static int smooth_poly_curve( curve_point *cp_ptr, void *data_ptr ){
  smooth_poly_data *info;

/* cast the data to the right type */
  info = (smooth_poly_data *) data_ptr;

/* evaluate the shape */
  eval_polygon(info->x1, info->dxdr, info->r, info->vx0, info->n, 
		 info->sharpness, cp_ptr->r, &(cp_ptr->x), &(cp_ptr->xr), &(cp_ptr->xrr));
  eval_polygon(info->y1, info->dydr, info->r, info->vy0, info->n, 
		 info->sharpness, cp_ptr->r, &(cp_ptr->y), &(cp_ptr->yr), &(cp_ptr->yrr));

  return 0;
}


static void *cleanup_smooth_poly_curve( void *data_ptr ){
  smooth_poly_data *info;

/* cast the data to the right type */
  info = (smooth_poly_data *) data_ptr;

  if (info != NULL){
    free_poly_curve_arrays( info );
    free(info);
  }

  return NULL;
}

static void alloc_poly_curve_arrays( smooth_poly_data *info ){
/* node points */
  info->r = vector(0,info->n+1);
/* x-coordinate */
  info->x    = vector(1,info->n);
  info->vx0  = vector(0,info->n);
  info->dxdr = vector(0,info->n);
/* y-coordinate */
  info->y    = vector(1,info->n);
  info->vy0  = vector(0,info->n);
  info->dydr = vector(0,info->n);
}

static void free_poly_curve_arrays( smooth_poly_data *info ){
/* node points */
  free_vector( info->r, 0, info->n+1 );
/* x-coordinate */
  free_vector( info->x,    1, info->n );
  free_vector( info->vx0,  0, info->n );
  free_vector( info->dxdr, 0, info->n );
/* y-coordinate */
  free_vector( info->y,    1, info->n );
  free_vector( info->vy0,  0, info->n );
  free_vector( info->dydr, 0, info->n );
}


static void eval_polygon(real x1, real *slope, real *r_node, real *v0, 
			 int n, real sharpness, real r, real *x, real *x_r, 
			 real *x_rr){
  int i;

/* function */
  *x = x1;
/* assume normalized input in the range 0 <= r<= 1.0 */
  r = r*r_node[n];
  for (i=0; i<=n; i++)
    *x += log_cosh(slope[i], r_node[i], r_node[i+1], sharpness, r) - v0[i];

/* derivative */
  *x_r = 0.0;
/* assume normalized input in the range 0 <= r<= 1.0 */
  for (i=0; i<=n; i++)
    *x_r += 0.5 * slope[i] * (tanh(sharpness*(r-r_node[i])) - 
			    tanh(sharpness*(r-r_node[i+1])) );
  *x_r = r_node[n]* *x_r;

/* second derivative */
  *x_rr = 0.0;
  for (i=0; i<=n; i++)
    *x_rr += 0.5 * slope[i] * sharpness * 
      (1.0/cosh(sharpness*(r-r_node[i]))/cosh(sharpness*(r-r_node[i])) - 
       1.0/cosh(sharpness*(r-r_node[i+1]))/cosh(sharpness*(r-r_node[i+1])) );
  *x_rr = r_node[n]*r_node[n]* *x_rr;
}


static void 
write_smooth_poly( int32 dir, void *smooth_poly_data_ptr ){
  smooth_poly_data *info;

  info = (smooth_poly_data *) smooth_poly_data_ptr;

/* write the size */
  hput_int( info->n, "n", dir );

/* write the vectors */
  hput_vector( info->r,    0, info->n+1, "r",    dir );
  hput_vector( info->x,    1, info->n,   "x",    dir );
  hput_vector( info->vx0,  0, info->n,   "vx0",  dir );
  hput_vector( info->dxdr, 0, info->n,   "dxdr", dir );
  hput_vector( info->y,    1, info->n,   "y",    dir );
  hput_vector( info->vy0,  0, info->n,   "vy0",  dir );
  hput_vector( info->dydr, 0, info->n,   "dydr", dir );

  hput_real( info->x1, "x1", dir );
  hput_real( info->y1, "y1", dir );
  hput_real( info->sharpness, "sharpness", dir );

}


static void plot_smooth_poly( generic_curve *curve_ptr ){
  smooth_poly_data *info;
  int i;
  char cflag[10];
  real x_new, y_new;

  info = (smooth_poly_data *) curve_ptr->curve_data_ptr;

  for (i=1; i<= info->n; i++){
    rot_trans_scale( info->x[i], info->y[i], &x_new, &y_new, curve_ptr );
    PL_marker( x_new, y_new, CIRCLE);
    sprintf(cflag, "%i", i);
    PL_plot_string(cflag, x_new, y_new, 1, 1, 3);
  }
}

void *init_smooth_poly_curve( input_output *io_ptr, generic_curve *curve_ptr,
			     char *status ){
  smooth_poly_data *info;
  
  int i, n, level=3, no_indent=0;
  real prev_x, prev_y;
  const real eps = 1.e-7;
  char question[80];
  
/* initialize a dummy grid */

  if (strcmp("init", status)==0){
    if ((info = (smooth_poly_data *) malloc( sizeof(smooth_poly_data ) )) == NULL)
      printf("memory error in init_smooth_poly\n");
    info->n = n = 2;
/* allocate space */
    alloc_poly_curve_arrays( info );
/* sharpness of corners and width change*/
    info->sharpness = 5.0;
/* trivial but valid curve */
    info->x[1]=0.0; info->y[1]=0.0;
    info->x[2]=1.0; info->y[2]=0.0;
  }

/* assign all corners and widths */

  else if (strcmp("new", status)==0){
/* check that at least some things are ok */
    if (curve_ptr->curve_mapping!=smooth_poly_curve ||
	curve_ptr->curve_data_ptr==NULL){
      printf("init_smooth_poly_curve: status is `new' but the pointers are screwed\n");
      exit(1);
    }
/* look up the existing structure */
    info = (smooth_poly_data *) curve_ptr->curve_data_ptr;
/* remove old corners and widths */
    free_poly_curve_arrays( info );
/* read number of node points */
    info->n = n = 
      int_max(2, get_int(io_ptr , "Enter number of corners >= 2: ", 2,
			 no_indent ));
/* allocate space */
    alloc_poly_curve_arrays( info );

/* read the corners */
    prev_x = -1.0; prev_y = -1.0;
    for (i=1; i<=n; i++){
      sprintf(question, "corner %i: x: ", i);
      info->x[i] = get_real( io_ptr, question, prev_x + 1.0, level);
      sprintf(question, "corner %i: y: ", i);
      info->y[i] = get_real( io_ptr, question, prev_y + 1.0, level);
/* take the global transformation into account */
      inv_rot_trans_scale( info->x[i], info->y[i], 
			  &(info->x[i]), &(info->y[i]), curve_ptr );
      while (i > 1 && real_max(fabs(info->x[i] - prev_x), 
			       fabs(info->y[i] - prev_y) ) < eps ){
	printf("Input Error: consecetive corner points may not be identical!\n");
	sprintf(question, "corner %i: x: ", i);
	info->x[i] = get_real( io_ptr, question, prev_x + 1.0, level);
	sprintf(question, "corner %i: y: ", i);
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
/* corner sharpness */
    info->sharpness = 5.0;
  }
  else{
    printf("init_smooth_poly_curve: Unknown status value.\n");
    exit(1);
  }

  curve_ptr->curve_type = 2;
  curve_ptr->curve_data_ptr = (void *) info;
  curve_ptr->curve_mapping = smooth_poly_curve;
  curve_ptr->cleanup_curve_data = cleanup_smooth_poly_curve;
  curve_ptr->set_curve = set_smooth_polygon;
  curve_ptr->write_specific_curve = write_smooth_poly;
  curve_ptr->plot_control_points = plot_smooth_poly;
  curve_ptr->copy_nodes = copy_polygon_nodes;
  curve_ptr->copy_curve_data = copy_smooth_poly;

  curve_ptr->r_step = 0.01;

  compute_smooth_poly_curve( info, curve_ptr );

  return info;

}

void 
read_smooth_poly( int32 dir, generic_curve *curve_ptr ){
  smooth_poly_data *info;
  int low, high;

  info = (smooth_poly_data *) malloc( sizeof(smooth_poly_data) );
  curve_ptr->curve_data_ptr = (void *) info;
  curve_ptr->curve_mapping = smooth_poly_curve;
  curve_ptr->cleanup_curve_data = cleanup_smooth_poly_curve;
  curve_ptr->set_curve = set_smooth_polygon;
  curve_ptr->write_specific_curve = write_smooth_poly;
  curve_ptr->plot_control_points = plot_smooth_poly;
  curve_ptr->copy_nodes = copy_polygon_nodes;
  curve_ptr->copy_curve_data = copy_smooth_poly;

  hget_int( &(info->n), "n", dir );

/* read the vectors */
  info->r    = hget_vector( "r",    &low, &high, dir );
  info->x    = hget_vector( "x",    &low, &high, dir );
  info->vx0  = hget_vector( "vx0",  &low, &high, dir );
  info->dxdr = hget_vector( "dxdr", &low, &high, dir );
  info->y    = hget_vector( "y",    &low, &high, dir );
  info->vy0  = hget_vector( "vy0",  &low, &high, dir );
  info->dydr = hget_vector( "dydr", &low, &high, dir );

  hget_real( &(info->x1), "x1", dir );
  hget_real( &(info->y1), "y1", dir );
  hget_real( &(info->sharpness), "sharpness", dir );
}


static real *copy_polygon_nodes( generic_curve *curve_ptr, int *n ){
  smooth_poly_data *info;
  real *nodes;
  int i;

  info = (smooth_poly_data *) curve_ptr->curve_data_ptr;

/* allocate a node vector */
  nodes = vector( 1, info->n );

/* copy node-points from the geometry description */
  for (i=1; i<=info->n; i++){ 
    nodes[i] = info->r[i]/info->r[info->n]; 
  } 

  *n = info->n;
  return nodes;
}


static void *copy_smooth_poly( void *curve_data_ptr ){
  smooth_poly_data *info, *old_info;
  int i;

  info = (smooth_poly_data *) malloc( sizeof(smooth_poly_data ) );

  old_info = (smooth_poly_data *) curve_data_ptr;

  info->n         = old_info->n;
  info->x1        = old_info->x1;
  info->y1        = old_info->y1;
  info->sharpness = old_info->sharpness;

/* allocate the arrays */
  alloc_poly_curve_arrays( info );

/* copy the array elements */
  for (i=1; i<=info->n; i++){
    info->r[i]    = old_info->r[i];

    info->x[i]    = old_info->x[i];
    info->vx0[i]  = old_info->vx0[i];
    info->dxdr[i] = old_info->dxdr[i];

    info->y[i]    = old_info->y[i];
    info->vy0[i]  = old_info->vy0[i];
    info->dydr[i] = old_info->dydr[i];
  }

/* end points */
  info->r[0] = old_info->r[0]; info->r[info->n+1] = old_info->r[info->n+1];

  info->vx0[0]  = old_info->vx0[0];
  info->dxdr[0] = old_info->dxdr[0];
  
  info->vy0[0]  = old_info->vy0[0];
  info->dydr[0] = old_info->dydr[0];

  return info;
}
