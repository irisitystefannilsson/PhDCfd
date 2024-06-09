#include "curves.h"
#include "hdf_stuff.h"

generic_curve *
read_curve( int32 curve_dir ){

  generic_curve *curve_ptr;
  int32 specific_dir;

/* create a curve structure */
  curve_ptr = (generic_curve *) malloc( sizeof(generic_curve) );

/* it is important to initialize the pointers */
  curve_ptr->next = curve_ptr->prev = NULL;

/* read the name */
  curve_ptr->curve_name = hget_string( "curve_name", curve_dir );

  hget_int( &(curve_ptr->curve_type), "curve_type", curve_dir );
  hget_int( &(curve_ptr->periodic), "periodic", curve_dir );

  hget_real( &(curve_ptr->x_min_0), "x_min_0", curve_dir );
  hget_real( &(curve_ptr->x_max_0), "x_max_0", curve_dir );
  hget_real( &(curve_ptr->y_min_0), "y_min_0", curve_dir );
  hget_real( &(curve_ptr->y_max_0), "y_max_0", curve_dir );

  hget_real( &(curve_ptr->x_min), "x_min", curve_dir );
  hget_real( &(curve_ptr->x_max), "x_max", curve_dir );
  hget_real( &(curve_ptr->y_min), "y_min", curve_dir );
  hget_real( &(curve_ptr->y_max), "y_max", curve_dir );

  hget_int( &(curve_ptr->used_by_grid), "used_by_grid", curve_dir );
  hget_int( &(curve_ptr->plot_mode), "plot_mode", curve_dir );

  hget_real( &(curve_ptr->xyab[0][0]), "xyab00", curve_dir);
  hget_real( &(curve_ptr->xyab[1][0]), "xyab10", curve_dir);
  hget_real( &(curve_ptr->xyab[0][1]), "xyab01", curve_dir);
  hget_real( &(curve_ptr->xyab[1][1]), "xyab11", curve_dir);

/* color info */
  hget_double( &(curve_ptr->curve_color.r), "curve_color.r", curve_dir );
  hget_double( &(curve_ptr->curve_color.g), "curve_color.g", curve_dir );
  hget_double( &(curve_ptr->curve_color.b), "curve_color.b", curve_dir );

  hget_real( &(curve_ptr->r_step), "r_step", curve_dir );
  hget_int( &(curve_ptr->n_ticks), "n_ticks", curve_dir );

/* rotation, translation, scaling */
  hget_real( &(curve_ptr->theta_rot), "theta_rot", curve_dir );
  hget_real( &(curve_ptr->cos_theta), "cos_theta", curve_dir );
  hget_real( &(curve_ptr->sin_theta), "sin_theta", curve_dir );
  hget_real( &(curve_ptr->x_trans), "x_trans", curve_dir );
  hget_real( &(curve_ptr->y_trans), "y_trans", curve_dir );
  hget_real( &(curve_ptr->scaling), "scaling", curve_dir );

/* get the directory for the curve-specific parameters */
  if ((specific_dir = locate_dir("specific_params", curve_dir)) != -1){

/* assign the appropriate curve-specific pointers */
    switch( curve_ptr->curve_type ){

    case 0:
/* straight line */
      read_straight_line( specific_dir, curve_ptr );
      break;

    case 1:
/* cubic spline */
      read_cubic_spline( specific_dir, curve_ptr );
      break;

    case 2:
/* smoothed polygon */
      read_smooth_poly( specific_dir, curve_ptr );
      break;

    case 3:
/* circular arc */
      read_circular_arc( specific_dir, curve_ptr );
      break;

    default:
      printf("Unknown curve type: %i in read_curve\n", curve_ptr->curve_type);
/* free the memory */
      curve_ptr->cleanup_curve_data = NULL;
      delete_generic_curve( curve_ptr );
      curve_ptr = NULL;

    } /* end switch */

/* terminate access to the specific dir */
    Vdetach( specific_dir );
  } /* end if locate dir */

  return curve_ptr;
} /* end read_curve */

void 
write_curve( int32 curve_dir, generic_curve *curve_ptr ){

  int32 specific_dir;

/* save the name */
  hput_string(curve_ptr->curve_name, "curve_name", curve_dir);

  hput_int(curve_ptr->curve_type, "curve_type", curve_dir );

  hput_int(curve_ptr->periodic, "periodic", curve_dir );

  hput_real(curve_ptr->x_min_0, "x_min_0", curve_dir );
  hput_real(curve_ptr->x_max_0, "x_max_0", curve_dir );
  hput_real(curve_ptr->y_min_0, "y_min_0", curve_dir );
  hput_real(curve_ptr->y_max_0, "y_max_0", curve_dir );

  hput_real(curve_ptr->x_min, "x_min", curve_dir );
  hput_real(curve_ptr->x_max, "x_max", curve_dir );
  hput_real(curve_ptr->y_min, "y_min", curve_dir );
  hput_real(curve_ptr->y_max, "y_max", curve_dir );

  hput_int(curve_ptr->used_by_grid, "used_by_grid", curve_dir );
  hput_int(curve_ptr->plot_mode, "plot_mode", curve_dir );

  hput_real(curve_ptr->xyab[0][0], "xyab00", curve_dir);
  hput_real(curve_ptr->xyab[1][0], "xyab10", curve_dir);
  hput_real(curve_ptr->xyab[0][1], "xyab01", curve_dir);
  hput_real(curve_ptr->xyab[1][1], "xyab11", curve_dir);

/* color info */
  hput_double(curve_ptr->curve_color.r, "curve_color.r", curve_dir );
  hput_double(curve_ptr->curve_color.g, "curve_color.g", curve_dir );
  hput_double(curve_ptr->curve_color.b, "curve_color.b", curve_dir );

  hput_real(curve_ptr->r_step, "r_step", curve_dir );
  hput_int(curve_ptr->n_ticks, "n_ticks", curve_dir );

/* rotation, translation, scaling */
  hput_real(curve_ptr->theta_rot, "theta_rot", curve_dir );
  hput_real(curve_ptr->cos_theta, "cos_theta", curve_dir );
  hput_real(curve_ptr->sin_theta, "sin_theta", curve_dir );
  hput_real(curve_ptr->x_trans, "x_trans", curve_dir );
  hput_real(curve_ptr->y_trans, "y_trans", curve_dir );
  hput_real(curve_ptr->scaling, "scaling", curve_dir );

/* make a directory for the curve-specific parameters */
  if ((specific_dir = create_dir("specific_params", "specific_params", 
				 curve_dir)) != -1){
/* save the curve-type dependent data */
    (*curve_ptr->write_specific_curve)( specific_dir, curve_ptr->curve_data_ptr );
    Vdetach(specific_dir);
  }
  else
    printf("Error: save_curve: unable to create a directory for the type-specific "
	   "stuff\n");

}



generic_curve *choose_curve( input_output *io_ptr, char *name ){
  char prompt[80];
  int level=1;
  int icom, quit;
  generic_curve *curve_ptr=NULL;

#include "choose_curve_com.h"

  sprintf(prompt,"set curve type %s>", name);

  quit = 0;
  do{

    icom = get_command(io_ptr, prompt, command, ncom, level,
		       save_on_copy, argument);

    switch (icom) {

/* curves */
    case 0: 
/* make curve */
      curve_ptr = new_generic_curve(name);
/* straight line */
      init_straight_line( curve_ptr, 0.0, 0.0, 1.0, 0.0 );
      quit = 1;
      break;

    case 1: 
/* make curve */
      curve_ptr = new_generic_curve(name);
/* circular arc */
      init_circular_arc( curve_ptr, 0.0, M_PI/3.0, 1.0, 0.0, 0.0 );
      quit = 1;
      break;

    case 2: 
/* make curve */
      curve_ptr = new_generic_curve(name);
/* cubic spline */
      init_cubic_spline( io_ptr, curve_ptr, "init" );
      quit = 1;
      break;

    case 3: 
/* make curve */
      curve_ptr = new_generic_curve(name);
/* smooth polygon */
      init_smooth_poly_curve( io_ptr, curve_ptr, "init" );
      quit = 1;
      break;

    case 4:
      while ( (icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				  command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );

      break;

    case 5:
/* cancel, unsucessful curve generation */
      return NULL;
      break;
    default:
      ;
    }
  }
  while(!quit);

/* call the specific setup routine to enter all the curve specific data */
  set_curve( io_ptr, curve_ptr ); 


/* sucessful curve generation */
  return curve_ptr;
#undef ncom
}



generic_curve *new_generic_curve(char *name){
  generic_curve *curve_ptr;

  if ((curve_ptr = (generic_curve *)malloc(sizeof(generic_curve))) == NULL)
    printf("memory error in new_generic_curve\n");
  
/* initialize all fields */
  curve_ptr->prev = NULL;
  curve_ptr->next = NULL;

  if ((curve_ptr->curve_name = (char *) malloc( (strlen(name)+1) * sizeof(char) )) == NULL)
    printf("memory error in new_generic_curve, curve_name\n");
  curve_ptr->curve_name = strcpy( curve_ptr->curve_name, name );

  curve_ptr->curve_type = 0;
  curve_ptr->periodic = 0;

  curve_ptr->x_min_0 = 0.0;
  curve_ptr->x_max_0 = 1.0;
  curve_ptr->y_min_0 = 0.0;
  curve_ptr->y_max_0 = 1.0;

  curve_ptr->x_min = 0.0;
  curve_ptr->x_max = 1.0;
  curve_ptr->y_min = 0.0;
  curve_ptr->y_max = 1.0;

  curve_ptr->curve_data_ptr = NULL;
  curve_ptr->curve_mapping = NULL;
  curve_ptr->cleanup_curve_data = NULL;
  curve_ptr->set_curve = NULL;
  curve_ptr->write_specific_curve = NULL;
  curve_ptr->plot_control_points = NULL;
  curve_ptr->copy_nodes = NULL;

  curve_ptr->used_by_grid = 0;

/* set the plot mode */
  curve_ptr->plot_mode = 1 + 2 + 16;

  curve_ptr->xyab[0][0] = 0.0;
  curve_ptr->xyab[0][1] = 1.0;
  curve_ptr->xyab[1][0] = 0.0;
  curve_ptr->xyab[1][1] = 1.0;

  curve_ptr->r_step = 0.05;

  curve_ptr->n_ticks = 20;

/* get a color */
  get_new_color( &(curve_ptr->curve_color), 'c' );

/* rotation, translation, scaling */
  curve_ptr->theta_rot = 0.0;
  curve_ptr->cos_theta = 1.0;
  curve_ptr->sin_theta = 0.0;
  curve_ptr->x_trans = 0.0;
  curve_ptr->y_trans = 0.0;
  curve_ptr->scaling = 1.0;

  return curve_ptr;

}


void 
delete_generic_curve( generic_curve *curve_ptr){

/* deallocate the name */
  free( curve_ptr->curve_name );

/* free the data associated with the specific curve */
  if (curve_ptr->cleanup_curve_data != NULL){
    curve_ptr->curve_data_ptr = 
      (*curve_ptr->cleanup_curve_data)(curve_ptr->curve_data_ptr);
  }
/* free the structure itself */
  free(curve_ptr);

}



void plot_curve( generic_curve *curve_ptr, int plot_mode ){
  curve_point cp;
  int i;
  real sigma, sx, sy, pi, theta, px, py, mx, my, length, xp, yp, xm, ym
    , r_step;

/* plot_mode controls what is plotted according to:
   plot_mode = 1   the curve
   plot_mode = 2   a directional arrow
   plot_mode = 4   tickmarks to indicate parametrization
   plot_mode = 8   Name of the grid on the plot
   plot_mode = 16  Plot control points, i.e. node points for cubic splines, corners 
                   for smoothed polygons, start and endpoints for straight lines, and
                   origin, radius and angles for circular arcs
   Different options are combined by adding, i.e. plot_mode=1+2 gives both
   option 1 and 2. */

/* set the color */
  PL_color( curve_ptr->curve_color );

/* draw the curve in solid */
  if ((plot_mode & 1) && curve_ptr->curve_mapping != NULL){
    cp.r = 0.0;
    curve_function( &cp, curve_ptr );
    PL_move(cp.x, cp.y);
    for (cp.r = curve_ptr->r_step; cp.r < 1.0;
	 cp.r += curve_ptr->r_step){
      curve_function( &cp, curve_ptr );
      PL_draw(cp.x , cp.y);
    }
/* draw the last segment */
    cp.r = 1.0;
    curve_function( &cp, curve_ptr );
    PL_draw(cp.x , cp.y);
  }

/* draw an arrow that indicates the direction of the curve */
  if ((plot_mode & 2) && curve_ptr->curve_mapping != NULL){
    cp.r = 0.5;
    curve_function( &cp, curve_ptr );
    sigma = 1.0/sqrt(cp.xr*cp.xr + cp.yr*cp.yr);
    sx = sigma*cp.xr;
    sy = sigma*cp.yr;
    pi = 4.0*atan(1.0);
    theta = 20.0*pi/180.0;
    px = cos(theta)*sx - sin(theta)*sy;
    py = sin(theta)*sx + cos(theta)*sy;
    theta = -theta;
    mx = cos(theta)*sx - sin(theta)*sy;
    my = sin(theta)*sx + cos(theta)*sy;
    length = 0.05*real_max(curve_ptr->x_max-curve_ptr->x_min, 
			   curve_ptr->y_max-curve_ptr->y_min);
    xp = cp.x - length*px;
    yp = cp.y - length*py;
    xm = cp.x - length*mx;
    ym = cp.y - length*my;
    PL_move(xp, yp);
    PL_draw(cp.x, cp.y);
    PL_draw(xm, ym);
  }

/* draw tickmarks to indicate the parametrization */
  if ((plot_mode & 4) && curve_ptr->curve_mapping != NULL){
    length = 0.01*real_max(curve_ptr->x_max-curve_ptr->x_min, 
			   curve_ptr->y_max-curve_ptr->y_min);
    r_step = 1.0/(curve_ptr->n_ticks-1);
    for (i=1; i<= curve_ptr->n_ticks; i++){
      cp.r = (i-1) * r_step;
      curve_function( &cp, curve_ptr );
      sigma = 1.0/sqrt(cp.xr*cp.xr + cp.yr*cp.yr);
/* tangent */
      sx = sigma*cp.xr;
      sy = sigma*cp.yr;
/* start and end points */
      xp = cp.x - length*sy;
      yp = cp.y + length*sx;
      xm = cp.x + length*sy;
      ym = cp.y - length*sx;
      PL_move(xp, yp);
      PL_draw(xm, ym);
    }
  }

  if (plot_mode & 8)
    PL_title( 2, curve_ptr->curve_name);

  if (plot_mode & 16)
    plot_control_points( curve_ptr );

} /* end curve_plot */



void set_curve_plot_mode( input_output *io_ptr, generic_curve *curve_ptr ){

  char *prompt, *file_name;
  int level=3;
  const int save_command=1;
  int icom, quit, replot, col_flag;
  static int coordinate_axis=1;
  FILE *fp;

#include "set_curve_plot_com.h"

  prompt = "curve plot-mode>";

  quit = 0;
  replot = 0;
  PL_window(2, curve_ptr->xyab);
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
      if (!(curve_ptr->plot_mode & 1))
	curve_ptr->plot_mode = curve_ptr->plot_mode | 1;
      else
	curve_ptr->plot_mode = curve_ptr->plot_mode & ~1;
      break;

    case 1:
      if (!(curve_ptr->plot_mode & 2))
	curve_ptr->plot_mode = curve_ptr->plot_mode | 2;
      else
	curve_ptr->plot_mode = curve_ptr->plot_mode & ~2;
      break;

    case 2:
      if (!(curve_ptr->plot_mode & 4))
	curve_ptr->plot_mode = curve_ptr->plot_mode | 4;
      else
	curve_ptr->plot_mode = curve_ptr->plot_mode & ~4;
      break;

    case 3:
      if (!(curve_ptr->plot_mode & 16))
	curve_ptr->plot_mode = curve_ptr->plot_mode | 16;
      else
	curve_ptr->plot_mode = curve_ptr->plot_mode & ~16;
      break;

    case 4:
/* toggle coordinate_axis */
      coordinate_axis = !coordinate_axis;
      PL_scale( 2, coordinate_axis );
      break;

    case 5:
      if (!(curve_ptr->plot_mode & 8))
	curve_ptr->plot_mode = curve_ptr->plot_mode | 8;
      else
	curve_ptr->plot_mode = curve_ptr->plot_mode & ~8;
      break;

    case 6:
/* save postscript */
      if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
				 "curve.ps", &file_name, 'w', 0, 
				 save_command) ) != NULL){
	col_flag = get_yes_no( io_ptr, "Do you want to save color "
				 "information>", level );
	PL_postscript( 2, fp, col_flag );
	fclose( fp );
      }
      replot = 0;
      break;

    case 7:
      printf("Curve Plotmode:\n");
      printf("Curve: %s\n", (curve_ptr->plot_mode & 1)? "on" : "off");
      printf("Directional arrows: %s\n", (curve_ptr->plot_mode & 2)? "on" : "off");
      printf("Parametrization tickmarks: %s\n", 
	     (curve_ptr->plot_mode & 4)? "on" : "off");
      printf("Coordinate axis: %s\n", (curve_ptr->plot_mode & 8)? "on" : "off");
      printf("Control points: %s\n", (curve_ptr->plot_mode & 16)? "on" : "off");
/* don't redraw the plot */
      replot = 0;
      break;

    case 8:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case 9:
      quit = 1;
/* don't redraw the plot */
      replot = 0;
      break;

    default:
/* don't redraw the plot */
      replot = 0;
    }


  }
  while(!quit);
#undef ncom
}



void show_curve_parameters( generic_curve *curve_ptr ){
/* curve name */
  printf("Curve name: `%s'\n",curve_ptr->curve_name);

/* mapping type */
  printf("\t");
  switch( curve_ptr->curve_type ){

  case 0:
    printf("Straight line\n");
    break;

  case 1:
    printf("Cubic spline\n");
    break;

  case 2:
    printf("Smoothed polygon\n");
    break;

  case 3:
    printf("Circular arc\n");
    break;

  default:
    printf("Unknown curve type\n");
  }

/* one separator line */
/*  printf("\n"); */
}


generic_curve *get_curve_ptr( input_output *io_ptr, generic_curve *first_c){
  generic_curve *curve_ptr;
  char **command;
  int i, icom, ncom, *save_on_copy=NULL, *argument=NULL, level=0, quit, n_curves;

  if (first_c == NULL) return NULL;

/* count the number of curves */
  n_curves = 0;
  for (curve_ptr = first_c; curve_ptr != NULL; curve_ptr = curve_ptr->next)
    n_curves++;

  ncom = n_curves + 1;
  if ((command = (char **) malloc( (ncom+1)*sizeof(char*) )) == NULL)
    printf("memory error in get_curve_ptr\n");

  i=0;
  for (curve_ptr=first_c; curve_ptr != NULL; 
       curve_ptr = curve_ptr->next){
    command[i] = curve_ptr->curve_name;
    i++;
  }
  command[ncom-1] = "help";
  command[ncom]   = "cancel";

  quit = 0;
  do{
    icom = get_command(io_ptr, "curve name: ", command, ncom, level, 
		       save_on_copy, argument);

/* default is to quit */
    quit = 1;

    if (icom == -1)
/* not a valid command */
      quit = 0;
    else if (icom == ncom-1){
/* help */
      printf(
" o  Enter one of the curve names if you wish to proceed.\n"
" o  Enter cancel if you do not wish to proceed. You will then return to the\n"
"    previous command level.\n");
      quit = 0;
    }

  }
  while (!quit);

/* return the pointer to the selected curve */
  i=0;
  for (curve_ptr=first_c; curve_ptr != NULL; 
       curve_ptr = curve_ptr->next){
    if (i==icom) break;
    i++;
  }

  free(command);

  return curve_ptr;

}


int get_2curve_ptr( input_output *io_ptr, generic_curve *first_c, 
		   generic_curve **curve_1_ptr, generic_curve **curve_2_ptr, 
		   int level){
  generic_curve *curve_ptr;
  char **command;
  int i, icom, ncom, *save_on_copy=NULL, *argument=NULL, quit, n_curves;

/* count the number of curves */
  n_curves = 0;
  for (curve_ptr = first_c; curve_ptr != NULL; curve_ptr = curve_ptr->next)
    n_curves++;

  ncom = n_curves + 1;
  if ((command = (char **) malloc( (ncom+1)*sizeof(char*) )) == NULL)
    printf("memory error in get_2curve_ptr\n");

  i=0;
  for (curve_ptr=first_c; curve_ptr != NULL; 
       curve_ptr = curve_ptr->next){
    command[i] = curve_ptr->curve_name;
    i++;
  }
  command[ncom-1] = "help";
  command[ncom]   = "cancel";

  quit = 0;
  do{
    icom = get_command(io_ptr, "s=0 curve name: ", command, ncom, level, 
		       save_on_copy, argument);

/* default is to quit */
    quit = 1;

    if (icom == -1)
/* not a valid command */
      quit = 0;
    else if (icom == ncom-1){
/* help */
      printf(
" o  If you wish to proceed, enter the name of the curve corresponding to s=0.\n"
" o  If you do not wish to proceed, enter cancel. You will then return to the\n"
"    previous command level.\n");
      quit = 0;
    }
    else if (icom == ncom){
/* cancel ? */
      free(command);
      return 0;
    }

  }
  while (!quit);

/* find the pointer to the selected curve */
  i=0;
  for (*curve_1_ptr=first_c; *curve_1_ptr != NULL; 
       *curve_1_ptr = (*curve_1_ptr)->next){
    if (i==icom) break;
    i++;
  }

/* the second curve must not be the same as the first */
  ncom += -1;
  i=0;
  for (curve_ptr=first_c; curve_ptr != NULL; 
       curve_ptr = curve_ptr->next){
    if (curve_ptr != *curve_1_ptr){
      command[i] = curve_ptr->curve_name;
      i++;
    }
  }
  command[ncom-1] = "help";
  command[ncom] = "cancel";

  quit = 0;
  do{
    icom = get_command(io_ptr, "s=1 curve name: ", command, ncom, level, 
		       save_on_copy, argument);

/* default is to quit */
    quit = 1;

    if (icom == -1)
/* not a valid command */
      quit = 0;
    else if (icom == ncom-1){
/* help */
      printf(
" o  If you wish to proceed, enter the name of the curve corresponding to s=1.\n"
" o  If you do not wish to proceed, enter cancel. You will then return to the\n"
"    previous command level.\n");
      quit = 0;
    }
    else if (icom == ncom){
/* cancel ? */
      free(command);
      return 0;
    }

  }
  while (!quit);

/* find the pointer to the selected curve */
  for (*curve_2_ptr=first_c; *curve_2_ptr != NULL; 
       *curve_2_ptr = (*curve_2_ptr)->next)
    if ( strcmp(command[icom], (*curve_2_ptr)->curve_name) == 0) break;

/* free the space taken by command */
  free(command);

/* sucessful completion */
  return 1;

}


int curve_function( curve_point *cp_ptr, generic_curve *curve_ptr ){
  curve_point cp_old;
  int msg;

  cp_old.r = cp_ptr->r;

/* take care of any errors in the specific curve routine and pass them on */
  if ( (msg=(*curve_ptr->curve_mapping)( &cp_old, curve_ptr->curve_data_ptr )) != 0)
    return msg;

/* the transformation is done in the order: rotation - translation - scaling */
  cp_ptr->x = curve_ptr->scaling* (curve_ptr->x_trans + 
				   cp_old.x * curve_ptr->cos_theta -
				   cp_old.y * curve_ptr->sin_theta);

  cp_ptr->y = curve_ptr->scaling* (curve_ptr->y_trans + 
				   cp_old.x * curve_ptr->sin_theta +
				   cp_old.y * curve_ptr->cos_theta);

/* first derivative */
  cp_ptr->xr = curve_ptr->scaling* (cp_old.xr * curve_ptr->cos_theta -
				    cp_old.yr * curve_ptr->sin_theta);

  cp_ptr->yr = curve_ptr->scaling* (cp_old.xr * curve_ptr->sin_theta +
				    cp_old.yr * curve_ptr->cos_theta);

/* second derivative */
  cp_ptr->xrr = curve_ptr->scaling* (cp_old.xrr * curve_ptr->cos_theta -
				     cp_old.yrr * curve_ptr->sin_theta);

  cp_ptr->yrr = curve_ptr->scaling* (cp_old.xrr * curve_ptr->sin_theta +
				     cp_old.yrr * curve_ptr->cos_theta);

  return 0;
}

void rot_trans_scale( real x_old, real y_old, real *x_new, real *y_new, 
		     generic_curve *curve_ptr ){

/* the transformation is done in the order: rotation - translation - scaling */
  *x_new = curve_ptr->scaling* (curve_ptr->x_trans + 
				x_old * curve_ptr->cos_theta -
				y_old * curve_ptr->sin_theta);

  *y_new = curve_ptr->scaling* (curve_ptr->y_trans + 
				x_old * curve_ptr->sin_theta +
				y_old * curve_ptr->cos_theta);

}


void inv_rot_trans_scale( real x_old, real y_old, real *x_new, real *y_new, 
			 generic_curve *curve_ptr ){

/* invert the rotation - translation - scaling */
  *x_new = (x_old/curve_ptr->scaling - curve_ptr->x_trans) * curve_ptr->cos_theta +
    (y_old/curve_ptr->scaling - curve_ptr->y_trans) * curve_ptr->sin_theta;

  *y_new =-(x_old/curve_ptr->scaling - curve_ptr->x_trans) * curve_ptr->sin_theta +
    (y_old/curve_ptr->scaling - curve_ptr->y_trans) * curve_ptr->cos_theta;

}


void set_curve_trans(input_output *io_ptr, generic_curve *curve_ptr){
  char prompt[80];
  int icom, quit, replot;
  const int level=1, no_indent=0, plot_mode=1+2+16;
  const real eps=1.e-7;

#include "set_curve_trans_com.h"

  sprintf(prompt, "%s: transformation>", curve_ptr->curve_name);
  quit = 0;
  replot = 1;
  do{

    if (replot){
/* update the view port */
      curve_view_port( curve_ptr );
      PL_window(2, curve_ptr->xyab);
/* plot the mapping */
      PL_erase(2);
      PL_start_plot(2);
      plot_curve( curve_ptr, plot_mode );
      PL_end_plot();
    }

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    replot = 1;
    switch (icom) {

    case 0: 
      curve_ptr->theta_rot = get_real( io_ptr, "Angle of rotation (degrees): ", 
				    curve_ptr->theta_rot, no_indent);
      curve_ptr->cos_theta = cos( M_PI/180.0 * curve_ptr->theta_rot );
      curve_ptr->sin_theta = sin( M_PI/180.0 * curve_ptr->theta_rot );
      break;

    case 1: 
      curve_ptr->x_trans = get_real( io_ptr, "Horizontal translation: ", 
				    curve_ptr->x_trans, no_indent);
      break;

    case 2: 
      curve_ptr->y_trans = get_real( io_ptr, "Vertical translation: ", 
				    curve_ptr->y_trans, no_indent);
      break;

    case 3:
      curve_ptr->scaling = real_max(eps, 
				  get_real( io_ptr, "Scaling (>0): ", 
				    curve_ptr->scaling, no_indent) );
      break;

    case 4:
      printf("Rotation angle = %f\n"
	     "Horizontal translation = %f\n"
	     "Vertical translation = %f\n"
	     "Scaling factor = %f\n", curve_ptr->theta_rot, curve_ptr->x_trans, 
	     curve_ptr->y_trans, curve_ptr->scaling);
      replot = 0;
      break;

    case 5:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      break;

    case 6:
      quit = 1;
      break;

    default:
      replot = 0;
    }
  }
  while(!quit);
#undef ncom
}


void curve_view_port( generic_curve *curve_ptr ){
  real x_pad, y_pad, x1, x2, x3, x4, y1, y2, y3, y4;

/* compute the transformed coordinates of the corners of the bounding box */
  x1 = curve_ptr->scaling* (curve_ptr->x_trans + 
			  curve_ptr->x_min_0 * curve_ptr->cos_theta -
			  curve_ptr->y_min_0 * curve_ptr->sin_theta);
  y1 = curve_ptr->scaling* (curve_ptr->y_trans + 
			  curve_ptr->x_min_0 * curve_ptr->sin_theta +
			  curve_ptr->y_min_0 * curve_ptr->cos_theta);

  x2 = curve_ptr->scaling* (curve_ptr->x_trans + 
			  curve_ptr->x_max_0 * curve_ptr->cos_theta -
			  curve_ptr->y_min_0 * curve_ptr->sin_theta);
  y2 = curve_ptr->scaling* (curve_ptr->y_trans + 
			  curve_ptr->x_max_0 * curve_ptr->sin_theta +
			  curve_ptr->y_min_0 * curve_ptr->cos_theta);

  x3 = curve_ptr->scaling* (curve_ptr->x_trans + 
			  curve_ptr->x_max_0 * curve_ptr->cos_theta -
			  curve_ptr->y_max_0 * curve_ptr->sin_theta);
  y3 = curve_ptr->scaling* (curve_ptr->y_trans + 
			  curve_ptr->x_max_0 * curve_ptr->sin_theta +
			  curve_ptr->y_max_0 * curve_ptr->cos_theta);

  x4 = curve_ptr->scaling* (curve_ptr->x_trans + 
			  curve_ptr->x_min_0 * curve_ptr->cos_theta -
			  curve_ptr->y_max_0 * curve_ptr->sin_theta);
  y4 = curve_ptr->scaling* (curve_ptr->y_trans + 
			  curve_ptr->x_min_0 * curve_ptr->sin_theta +
			  curve_ptr->y_max_0 * curve_ptr->cos_theta);

/* get the min and max from the transformed corners */
  curve_ptr->x_min = real_min( x1, real_min( x2, real_min( x3, x4 ) ) );
  curve_ptr->y_min = real_min( y1, real_min( y2, real_min( y3, y4 ) ) );
  
  curve_ptr->x_max = real_max( x1, real_max( x2, real_max( x3, x4 ) ) );
  curve_ptr->y_max = real_max( y1, real_max( y2, real_max( y3, y4 ) ) );
  
/* compute padding */  
  x_pad = 0.05*(curve_ptr->x_max - curve_ptr->x_min);
  y_pad = 0.05*(curve_ptr->y_max - curve_ptr->y_min);
  curve_ptr->xyab[0][0] = curve_ptr->x_min - x_pad;
  curve_ptr->xyab[0][1] = curve_ptr->x_max + x_pad;
  curve_ptr->xyab[1][0] = curve_ptr->y_min - y_pad;
  curve_ptr->xyab[1][1] = curve_ptr->y_max + y_pad;

}


generic_curve *copy_curve( char *name, generic_curve *old_curve_ptr){
  generic_curve *curve_ptr;

  curve_ptr = (generic_curve *)malloc(sizeof(generic_curve));
  
/* initialize all fields */
  curve_ptr->prev = NULL;
  curve_ptr->next = NULL;

/* copy the name */
  curve_ptr->curve_name = (char *) malloc( (strlen(name)+1)*sizeof(char) );
  curve_ptr->curve_name = strcpy( curve_ptr->curve_name, name );

/* set the grid type */
  curve_ptr->curve_type = old_curve_ptr->curve_type; 

  curve_ptr->periodic = old_curve_ptr->periodic; 

  curve_ptr->x_min_0 = old_curve_ptr->x_min_0;
  curve_ptr->x_max_0 = old_curve_ptr->x_max_0;
  curve_ptr->y_min_0 = old_curve_ptr->y_min_0;
  curve_ptr->y_max_0 = old_curve_ptr->y_max_0;

  curve_ptr->x_min = old_curve_ptr->x_min;
  curve_ptr->x_max = old_curve_ptr->x_max;
  curve_ptr->y_min = old_curve_ptr->y_min;
  curve_ptr->y_max = old_curve_ptr->y_max;

/* copy the specific data */
  curve_ptr->curve_data_ptr = 
    (*old_curve_ptr->copy_curve_data)( old_curve_ptr->curve_data_ptr );

/* save the pointers to the data and the functions for
   evaluating the curve */
  curve_ptr->curve_mapping = old_curve_ptr->curve_mapping;
  curve_ptr->cleanup_curve_data = old_curve_ptr->cleanup_curve_data;
  curve_ptr->set_curve = old_curve_ptr->set_curve;
  curve_ptr->write_specific_curve = old_curve_ptr->write_specific_curve;
  curve_ptr->plot_control_points = old_curve_ptr->plot_control_points;
  curve_ptr->copy_nodes = old_curve_ptr->copy_nodes;
  curve_ptr->copy_curve_data = old_curve_ptr->copy_curve_data;

/* the new curve can not yet be used by a mapping */
  curve_ptr->used_by_grid = 0;

  curve_ptr->plot_mode = old_curve_ptr->plot_mode;

  curve_ptr->xyab[0][0] = old_curve_ptr->xyab[0][0];
  curve_ptr->xyab[0][1] = old_curve_ptr->xyab[0][1];
  curve_ptr->xyab[1][0] = old_curve_ptr->xyab[1][0];
  curve_ptr->xyab[1][1] = old_curve_ptr->xyab[1][1];

/* get a new color */
  get_new_color( &(curve_ptr->curve_color), 'c' ); 

  curve_ptr->r_step = old_curve_ptr->r_step;
  curve_ptr->n_ticks = old_curve_ptr->n_ticks;

/* rotation, translation, scaling */
  curve_ptr->theta_rot = old_curve_ptr->theta_rot;
  curve_ptr->cos_theta = old_curve_ptr->cos_theta;
  curve_ptr->sin_theta = old_curve_ptr->sin_theta;
  curve_ptr->x_trans = old_curve_ptr->x_trans;
  curve_ptr->y_trans = old_curve_ptr->y_trans;
  curve_ptr->scaling = old_curve_ptr->scaling;

  return curve_ptr;

}

void set_curve( input_output *io_ptr, generic_curve *curve_ptr ){
  if (curve_ptr == NULL || curve_ptr->set_curve == NULL){
    printf("ERROR: set_curve was called with inconsistent pointers.\n");
  }
  else
    (*curve_ptr->set_curve)( io_ptr, curve_ptr );
}

real *copy_nodes( generic_curve *curve_ptr, int *n_nodes ){
  if (curve_ptr == NULL){
    *n_nodes = 0;
    return NULL;
  }
  else
    return (*curve_ptr->copy_nodes)( curve_ptr, n_nodes );
}

void plot_control_points( generic_curve *curve_ptr ){
  if (curve_ptr == NULL || curve_ptr->plot_control_points == NULL)
    return;
  else
    (*curve_ptr->plot_control_points)( curve_ptr );
}
