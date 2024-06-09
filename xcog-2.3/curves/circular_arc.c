#include "curves.h"

/* private member function declarations */
static void compute_circular_arc(circular_arc_data *info, 
				 generic_curve *curve_ptr );
static int circular_arc( curve_point *cp_ptr, void *circular_arc_data_ptr );
static void *cleanup_circular_arc( void *circular_arc_data_ptr );
static void 
write_circular_arc( int32 dir, void *circular_arc_data_ptr );
static void plot_circular_arc( generic_curve *curve_ptr );
static void set_circular_arc(input_output *io_ptr, generic_curve *curve_ptr);
static real *copy_circular_nodes( generic_curve *curve_ptr, int *n );
static void *copy_circular_arc( void *curve_data_ptr );
/* end private member function declarations */


static void set_circular_arc(input_output *io_ptr, generic_curve *curve_ptr){
  char prompt[80];
  const real eps=1.e-7;
  real pi, new_radius, xc, yc;
  int icom, quit, replot;
  circular_arc_data *info;
  int level=2, no_indent=0;

#include "circular_arc_com.h"

  pi = 4.0 * atan( 1.0 );

  sprintf(prompt, "%s: circular arc>", curve_ptr->curve_name);

  info = (circular_arc_data *) curve_ptr->curve_data_ptr;

  quit = 0;
  replot = 1;
  PL_window(2, curve_ptr->xyab);
  do{

    if (replot){
/* plot the curve */
      PL_erase(2);
      PL_start_plot(2);
      plot_curve( curve_ptr, curve_ptr->plot_mode );
      PL_end_plot();
    }

    icom = get_command(io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {
    case 0:
/* starting angle */
      info->theta_0 = pi/180 * 
	( get_real( io_ptr, "start at angle (degrees): ", 
		   info->theta_0*180.0/pi + curve_ptr->theta_rot , no_indent)
	 - curve_ptr->theta_rot );
      compute_circular_arc( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 1:
/* stopping angle */
      info->theta_1 = pi/180 * 
	( get_real( io_ptr, "stop at angle (degrees): ", 
		   info->theta_1*180.0/pi + curve_ptr->theta_rot , no_indent)
	 - curve_ptr->theta_rot );
      compute_circular_arc( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 2:
/* radius */
      do{
	new_radius = get_real( io_ptr, "radius>0: ", 
			      info->radius * curve_ptr->scaling, no_indent) 
	  / curve_ptr->scaling;
	if (new_radius < eps )
	  printf("Input Error: the radius must be positive!\n");
      } while (new_radius < eps );
      info->radius = new_radius;
      compute_circular_arc( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 3:
/* origin */
      rot_trans_scale( info->xcenter, info->ycenter, &xc, &yc, curve_ptr );
      info->xcenter = get_real( io_ptr, "x-center: ", xc, no_indent);
      info->ycenter = get_real( io_ptr, "y-center: ", yc, no_indent);
      inv_rot_trans_scale( info->xcenter, info->ycenter, 
			  &(info->xcenter), &(info->ycenter), curve_ptr );
      compute_circular_arc( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 4:
/* plot-mode */
      set_curve_plot_mode( io_ptr, curve_ptr );
      replot = 0;
      break;

    case 5:
/* show */
      printf("start-angle: %e, stop-angle: %e, radius: %e\n",
	     180/pi*info->theta_0 + curve_ptr->theta_rot, 
	     180/pi*info->theta_1 + curve_ptr->theta_rot, 
	     info->radius * curve_ptr->scaling);
      rot_trans_scale( info->xcenter, info->ycenter, &xc, &yc, curve_ptr );
      printf("x-center: %e, y-center: %e\n", xc, yc);
      replot = 0;
      break;

    case 6:
/* help */
      while ( (icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				  command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case 7:
/* exit */
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);

}

static int circular_arc( curve_point *cp_ptr, void *circular_arc_data_ptr ){
  circular_arc_data *info;
  real theta;
  
  info = (circular_arc_data *) circular_arc_data_ptr;

  theta = (1.0 - cp_ptr->r) * info->theta_0 + cp_ptr->r * info->theta_1;
  cp_ptr->x = info->xcenter + info->radius * cos( theta );
  cp_ptr->y = info->ycenter + info->radius * sin( theta );

  cp_ptr->xr =-info->radius * sin( theta) * (info->theta_1 - info->theta_0);
  cp_ptr->yr = info->radius * cos( theta) * (info->theta_1 - info->theta_0);

  cp_ptr->xrr = -info->radius * cos( theta ) * (info->theta_1 - info->theta_0) * 
    (info->theta_1 - info->theta_0);
  cp_ptr->yrr = -info->radius * sin( theta ) * (info->theta_1 - info->theta_0) * 
    (info->theta_1 - info->theta_0);

  return 0;
}

static void compute_circular_arc(circular_arc_data *info, 
			   generic_curve *curve_ptr ){
  const real eps=1.0e-7;
  curve_point cp;

/* periodic or not? */
  curve_ptr->periodic = (fabs(info->theta_1 - info->theta_0) + eps >= 
			 2.0*info->pi)? 1 : 0;

/* update the number of segments */
  curve_ptr->r_step = 1.0/real_max( 5.0, 50.0 * 
				   fabs(info->theta_1 - info->theta_0)/ (2.0*info->pi) );

/* find x,y - min,max */
  curve_ptr->x_min_0 = info->xcenter;
  curve_ptr->x_max_0 = info->xcenter;
  curve_ptr->y_min_0 = info->ycenter;
  curve_ptr->y_max_0 = info->ycenter;
  for (cp.r=0.0; cp.r <=1.0 + 0.5*curve_ptr->r_step; cp.r += curve_ptr->r_step){
    circular_arc( &cp, info );
    curve_ptr->x_min_0 = real_min( curve_ptr->x_min_0, cp.x );
    curve_ptr->x_max_0 = real_max( curve_ptr->x_max_0, cp.x );
    curve_ptr->y_min_0 = real_min( curve_ptr->y_min_0, cp.y );
    curve_ptr->y_max_0 = real_max( curve_ptr->y_max_0, cp.y );
  }

  curve_view_port( curve_ptr );
}

static void *cleanup_circular_arc( void *circular_arc_data_ptr ){
  circular_arc_data *info;
  
  info = (circular_arc_data *) circular_arc_data_ptr;
  if (info != NULL)
    free( info );

  return NULL;
}

static void 
write_circular_arc( int32 dir, void *circular_arc_data_ptr ){
  circular_arc_data *info;

  info = (circular_arc_data *) circular_arc_data_ptr;
  
  hput_real(info->pi, "pi", dir);
  hput_real(info->theta_0, "theta_0", dir);
  hput_real(info->theta_1, "theta_1", dir);
  hput_real(info->xcenter, "xcenter", dir);
  hput_real(info->ycenter, "ycenter", dir);
  hput_real(info->radius, "radius", dir);
}

static void plot_circular_arc( generic_curve *curve_ptr ){
  circular_arc_data *info;
  real xp, yp, theta, sx, sy, px, py, mx, my, length, xm, ym;
  real x_new, y_new, xc, yc;
  int i;

  info = (circular_arc_data *) curve_ptr->curve_data_ptr;

  rot_trans_scale(info->xcenter, info->ycenter, &x_new, &y_new, curve_ptr);
  PL_marker(x_new, y_new, CIRCLE);

/* draw ray to the start point */
  PL_move( x_new, y_new );
  xp = info->xcenter + info->radius * cos( info->theta_0 );
  yp = info->ycenter + info->radius * sin( info->theta_0 );
  rot_trans_scale( xp, yp, &x_new, &y_new, curve_ptr );
  PL_draw( x_new, y_new );
  PL_marker(x_new, y_new, INV_TRIANGLE);
  PL_plot_string("start", x_new, y_new, 1, -1, 3);

/* draw ray to the stop point */
  rot_trans_scale(info->xcenter, info->ycenter, &x_new, &y_new, curve_ptr);
  PL_move( x_new, y_new );
  xp = info->xcenter + info->radius * cos( info->theta_1 );
  yp = info->ycenter + info->radius * sin( info->theta_1 );
  rot_trans_scale( xp, yp, &x_new, &y_new, curve_ptr );
  PL_draw( x_new, y_new );
  PL_marker(x_new, y_new, INV_TRIANGLE);
  PL_plot_string("end", x_new, y_new, 1, 1, 3);

/* draw x-axis ( take global transformations into account )*/
  rot_trans_scale(info->xcenter, info->ycenter, &x_new, &y_new, curve_ptr);
  PL_move( x_new, y_new );
  x_new += 0.3*info->radius * curve_ptr->scaling;
  PL_draw( x_new, y_new );

/* draw an arrow on the x-axis */
  sx = 1.0;
  sy = 0.0;
  theta = 20.0*info->pi/180.0;
  px = cos(theta)*sx - sin(theta)*sy;
  py = sin(theta)*sx + cos(theta)*sy;
  theta = -theta;
  mx = cos(theta)*sx - sin(theta)*sy;
  my = sin(theta)*sx + cos(theta)*sy;
  length = 0.05*info->radius * curve_ptr->scaling;
  xp = x_new - length*px;
  yp = y_new - length*py;
  xm = x_new - length*mx;
  ym = y_new - length*my;

  PL_move( xp, yp);
  PL_draw( x_new, y_new );
  PL_draw( xm, ym);

/* mark the axis with an "x" */
  PL_plot_string("x", x_new, y_new, 1, 1, 3);

/* angle theta_0 */
  rot_trans_scale( info->xcenter, info->ycenter, &xc, &yc, curve_ptr );
  xp = xc + 0.2*info->radius * curve_ptr->scaling;
  yp = yc;
  PL_move( xp, yp );
  for (i=0; i<=40; i++){
    theta = i * (info->theta_0 + curve_ptr->theta_rot * ((double)M_PI)/180.0)/ 40.0;
    xp = xc + 0.2*info->radius * curve_ptr->scaling * cos( theta );
    yp = yc + 0.2*info->radius * curve_ptr->scaling * sin( theta );
    PL_draw( xp, yp );
  }

/* angle theta_1 */
  xp = xc + 0.1*info->radius * curve_ptr->scaling;
  yp = yc;
  PL_move( xp, yp );
  for (i=0; i<=40; i++){
    theta = i * (info->theta_1 + curve_ptr->theta_rot * ((double)M_PI)/180.0)/ 40.0;
    xp = xc + 0.1*info->radius * curve_ptr->scaling * cos( theta );
    yp = yc + 0.1*info->radius * curve_ptr->scaling * sin( theta );
    PL_draw( xp, yp );
  }
}

circular_arc_data *init_circular_arc(generic_curve *curve_ptr,
				       real theta_0, real theta_1, real radius,
				       real xcenter, real ycenter ){
  circular_arc_data *info;
  
  if ((info = (circular_arc_data *) malloc( sizeof(circular_arc_data ) )) == NULL)
    printf("memory error in init_circular_arc\n");

  info->pi = 4.0 * atan( 1.0 );
  info->theta_0 = theta_0;
  info->theta_1 = theta_1;
  info->radius  = radius;
  info->xcenter = xcenter;
  info->ycenter = ycenter;

  curve_ptr->curve_type = 3;
  curve_ptr->curve_data_ptr = (void *) info;
  curve_ptr->curve_mapping = circular_arc;
  curve_ptr->cleanup_curve_data = cleanup_circular_arc;
  curve_ptr->set_curve = set_circular_arc;
  curve_ptr->write_specific_curve = write_circular_arc;
  curve_ptr->plot_control_points = plot_circular_arc;
  curve_ptr->copy_nodes = copy_circular_nodes;
  curve_ptr->copy_curve_data = copy_circular_arc;

  compute_circular_arc( info, curve_ptr );

  return info;
}

void 
read_circular_arc( int32 dir, generic_curve *curve_ptr ){
  circular_arc_data *info;

  info = (circular_arc_data *) malloc( sizeof(circular_arc_data) );
  curve_ptr->curve_data_ptr = (void *) info;
  curve_ptr->curve_mapping = circular_arc;
  curve_ptr->cleanup_curve_data = cleanup_circular_arc;
  curve_ptr->set_curve = set_circular_arc;
  curve_ptr->write_specific_curve = write_circular_arc;
  curve_ptr->plot_control_points = plot_circular_arc;
  curve_ptr->copy_nodes = copy_circular_nodes;
  curve_ptr->copy_curve_data = copy_circular_arc;
  
  hget_real( &(info->pi), "pi", dir );
  hget_real( &(info->theta_0), "theta_0", dir );
  hget_real( &(info->theta_1), "theta_1", dir );
  hget_real( &(info->xcenter), "xcenter", dir );
  hget_real( &(info->ycenter), "ycenter", dir );
  hget_real( &(info->radius), "radius", dir );
}


static real *
copy_circular_nodes( generic_curve *curve_ptr, int *n ){
  real *nodes;

/* allocate a node vector */
  nodes = vector( 1, 2 );

/* there are only 2 node points for a straight line */
  nodes[1] = 0.0;
  nodes[2] = 1.0;

  *n = 2;
  return nodes;
}


static void *copy_circular_arc( void *curve_data_ptr ){
  circular_arc_data *info, *old_info;

  info = (circular_arc_data *) malloc( sizeof(circular_arc_data ) );

  old_info = (circular_arc_data *) curve_data_ptr;

  info->pi      = old_info->pi;
  info->theta_0 = old_info->theta_0;
  info->theta_1 = old_info->theta_1;
  info->xcenter = old_info->xcenter;
  info->ycenter = old_info->ycenter;
  info->radius  = old_info->radius;

  return info;

}
