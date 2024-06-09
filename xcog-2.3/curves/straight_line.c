#include "curves.h"

/* private prototypes */
static int straight_line( curve_point *cp_ptr, void *straight_line_data_ptr );
static void compute_straight_line(straight_line_data *info, 
			   generic_curve *curve_ptr );
static void *cleanup_straight_line( void *straight_line_data_ptr );
static void 
write_straight_line( int32 dir, void *straight_line_data_ptr );
static void plot_straight_line( generic_curve *curve_ptr );
static void set_straight_line(input_output *io_ptr, generic_curve *curve_ptr );
static real *copy_line_nodes( generic_curve *curve_ptr, int *n );
static void *copy_straight_line( void *curve_data_ptr );
/* end prototypes */

static void set_straight_line(input_output *io_ptr, generic_curve *curve_ptr){
  char prompt[80];
  const real eps=1.e-7;
  real new_x, new_y, xb, yb, xe, ye;
  int icom, quit, replot;
  straight_line_data *info;
  int level=2, no_indent=0;

#include "straight_line_com.h"

  sprintf(prompt, "%s: straight line>", curve_ptr->curve_name);

/* cast the data pointer to the right type */
  info = (straight_line_data *) curve_ptr->curve_data_ptr;

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
/* start point */
      rot_trans_scale( info->x_begin, info->y_begin, &xb, &yb, curve_ptr );
      new_x = get_real( io_ptr, "x-start: ", xb, no_indent);
      new_y = get_real( io_ptr, "y-start: ", yb, no_indent);
      inv_rot_trans_scale( new_x, new_y, &new_x, &new_y, curve_ptr );
      while (real_max(fabs(info->x_end - new_x), fabs(info->y_end - new_y) ) 
	     < eps ){
	printf("Input Error: the starting and ending points may not be");
	printf(" identical!\n");
	new_x = get_real( io_ptr, "x-start: ", xb, no_indent);
	new_y = get_real( io_ptr, "y-start: ", yb, no_indent);
	inv_rot_trans_scale( new_x, new_y, &new_x, &new_y, curve_ptr );
      }
      info->x_begin = new_x;
      info->y_begin = new_y;
      compute_straight_line( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 1:
/* end point */
      rot_trans_scale( info->x_end,   info->y_end,   &xe, &ye, curve_ptr );
      new_x = get_real( io_ptr, "x-end: ", xe, no_indent);
      new_y = get_real( io_ptr, "y-end: ", ye, no_indent);
      inv_rot_trans_scale( new_x, new_y, &new_x, &new_y, curve_ptr );
      while (real_max(fabs(info->x_begin - new_x),fabs(info->y_begin - new_y) )
	     < eps ){
	printf("Input Error: the starting and ending points may not be");
	printf(" identical!\n");
	new_x = get_real( io_ptr, "x-end: ", xe, no_indent);
	new_y = get_real( io_ptr, "y-end: ", ye, no_indent);
	inv_rot_trans_scale( new_x, new_y, &new_x, &new_y, curve_ptr );
      }
      info->x_end = new_x;
      info->y_end = new_y;
      compute_straight_line( info, curve_ptr );
      PL_window(2, curve_ptr->xyab);
      curve_ptr->plot_mode = curve_ptr->plot_mode | 3;
      break;

    case 2:
/* plot-mode */
      set_curve_plot_mode( io_ptr, curve_ptr );
      replot = 0;
      break;

    case 3:
/* show */
      rot_trans_scale( info->x_begin, info->y_begin, &xb, &yb, curve_ptr );
      rot_trans_scale( info->x_end,   info->y_end,   &xe, &ye, curve_ptr );
      printf("x-begin: %e, y-begin: %e,\n"
	     "x-end: %e, y-end: %e.\n", xb, yb, xe, ye);
      replot = 0;
      break;

    case 4:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case 5:
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);

}

static int straight_line( curve_point *cp_ptr, void *straight_line_data_ptr ){
  straight_line_data *info;
  
  info = (straight_line_data *) straight_line_data_ptr;

  cp_ptr->x = (1.0 - cp_ptr->r)*info->x_begin + cp_ptr->r*info->x_end;
  cp_ptr->y = (1.0 - cp_ptr->r)*info->y_begin + cp_ptr->r*info->y_end;

  cp_ptr->xr = info->x_end - info->x_begin;
  cp_ptr->yr = info->y_end - info->y_begin;

  cp_ptr->xrr = 0.0;
  cp_ptr->yrr = 0.0;

  return 0;
}


static void compute_straight_line(straight_line_data *info, 
			   generic_curve *curve_ptr ){

  curve_ptr->x_min_0 = real_min( info->x_begin, info->x_end );
  curve_ptr->x_max_0 = real_max( info->x_begin, info->x_end );
  curve_ptr->y_min_0 = real_min( info->y_begin, info->y_end );
  curve_ptr->y_max_0 = real_max( info->y_begin, info->y_end );

  curve_view_port( curve_ptr );
}

static void *cleanup_straight_line( void *straight_line_data_ptr ){
  straight_line_data *info;
  
  info = (straight_line_data *) straight_line_data_ptr;
  if (info != NULL)
    free( info );

  return NULL;
}

static void 
write_straight_line( int32 dir, void *straight_line_data_ptr ){
  straight_line_data *info;

  info = (straight_line_data *) straight_line_data_ptr;
  hput_real( info->x_begin, "x_begin", dir);
  hput_real( info->x_end, "x_end", dir);
  hput_real( info->y_begin, "y_begin", dir);
  hput_real( info->y_end, "y_end", dir);
}

static void plot_straight_line( generic_curve *curve_ptr ){
  straight_line_data *info;
  curve_point cp;

/* we have to change the calling structure here to make curve_ptr accessible */

  info = (straight_line_data *) curve_ptr->curve_data_ptr;

  cp.r = 0.0;
  curve_function( &cp, curve_ptr );
  PL_marker(cp.x, cp.y, TRIANGLE);
  PL_plot_string("start", cp.x, cp.y, 0, 1, 3);

  cp.r = 1.0;
  curve_function( &cp, curve_ptr );
  PL_marker(cp.x, cp.y, TRIANGLE);
  PL_plot_string("end", cp.x, cp.y, 0, 1, 3);
}

straight_line_data *init_straight_line(generic_curve *curve_ptr,
				       real x_begin, real y_begin, 
				       real x_end, real y_end ){
  straight_line_data *info;
  
  info = (straight_line_data *) malloc( sizeof(straight_line_data ) );

  info->x_begin = x_begin;
  info->y_begin = y_begin;
  info->x_end   = x_end;
  info->y_end   = y_end;

  curve_ptr->curve_type = 0;
  curve_ptr->curve_data_ptr = (void *) info;
  curve_ptr->curve_mapping = straight_line;
  curve_ptr->cleanup_curve_data = cleanup_straight_line;
  curve_ptr->set_curve = set_straight_line;
  curve_ptr->write_specific_curve = write_straight_line;
  curve_ptr->plot_control_points = plot_straight_line;
  curve_ptr->copy_nodes = copy_line_nodes;
  curve_ptr->copy_curve_data = copy_straight_line;

  curve_ptr->r_step = 0.05;

/* not periodic */
  curve_ptr->periodic = 0;

  compute_straight_line( info, curve_ptr );

  return info;
}

void 
read_straight_line( int32 dir, generic_curve *curve_ptr ){
  straight_line_data *info;

  curve_ptr->curve_mapping = straight_line;
  curve_ptr->cleanup_curve_data = cleanup_straight_line;
  curve_ptr->set_curve = set_straight_line;
  curve_ptr->write_specific_curve = write_straight_line;
  curve_ptr->plot_control_points = plot_straight_line;
  curve_ptr->copy_nodes = copy_line_nodes;
  curve_ptr->copy_curve_data = copy_straight_line;

  info = (straight_line_data *) malloc( sizeof(straight_line_data) );
  curve_ptr->curve_data_ptr = (void *) info;

  hget_real( &(info->x_begin), "x_begin", dir );
  hget_real( &(info->x_end),   "x_end", dir );
  hget_real( &(info->y_begin), "y_begin", dir );
  hget_real( &(info->y_end),   "y_end", dir );
}

static real *
copy_line_nodes( generic_curve *curve_ptr, int *n ){
  real *nodes;

/* allocate a node vector */
  nodes = vector( 1, 2 );

/* there are only 2 node points for a straight line */
  nodes[1] = 0.0;
  nodes[2] = 1.0;

  *n = 2;
  return nodes;
}


static void *
copy_straight_line( void *curve_data_ptr ){
  straight_line_data *info, *old_info;

  info = (straight_line_data *) malloc( sizeof(straight_line_data ) );

  old_info = (straight_line_data *) curve_data_ptr;

  info->x_begin = old_info->x_begin;
  info->x_end   = old_info->x_end;
  info->y_begin = old_info->y_begin;
  info->y_end   = old_info->y_end;

  return info;

}
