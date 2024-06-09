#include "stretchings.h"
#include "hdf_stuff.h"

generic_stretching *
choose_stretching(input_output *io_ptr, generic_curve *curve_ptr, int *grid_points ){
  char prompt[80];
  const int level=3;
  int icom, quit, i;
  generic_stretching *stretch_ptr=NULL;

#include "choose_stretching_com.h"

  if (curve_ptr == NULL)
    sprintf(prompt,"choose a stretching>");
  else
    sprintf(prompt,"choose a stretching %s %s>", "for", curve_ptr->curve_name );


  do{

    icom = get_command(io_ptr, prompt, command, ncom, level,
		       save_on_copy, argument);

    quit = 0;
    switch (icom) {

    case 0: 
/* arclength stretching */
      if (curve_ptr == NULL)
	printf("The arclength stretching needs a curve, but none is supplied!\n");
      else{
	stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );
	init_arclength_stretch( stretch_ptr, curve_ptr );
	quit = 1;
      }
      break;
      
    case 1:
/* exponential stretching */
      stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );
      init_exp_stretch( stretch_ptr, *grid_points );
      quit = 1;
      break;
      
    case 2:
/* layer stretching */
      stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );
      init_layer_stretch( stretch_ptr, curve_ptr );
      quit = 1;
      break;
      
    case 3:
/* tanh stretching */ 
      stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );
      init_tanh_stretch( stretch_ptr, *grid_points );
      quit = 1;
      break;

    case 4:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );

      break;

    case 5:
/* cancel */
      return NULL;
      break;

    default:
      ;
    }
  }
  while(!quit);

/* set the parameter value */
  set_stretching( io_ptr, stretch_ptr, curve_ptr, grid_points );

  return stretch_ptr;
}


char *
stretching_name( generic_stretching *stretch_ptr ){
  char *name;

  if (stretch_ptr == NULL){

    name = "no stretching";

  }
  else{
    
    switch( stretch_ptr->stretch_type ){

    case 1:
      name = "arclength stretching";
      break;

    case 2:
      name = "exponential stretching";
      break;

    case 3:
      name = "layer stretching";
      break;

    default:
      name = "unknown stretching";
    }

  }

  return name;
}


void 
plot_stretch(generic_stretching *stretch_ptr, generic_curve *curve_ptr, 
	     int grid_points, int plot_mode ){
  curve_point cp;
  int i;
  real sigma, sx, sy, pi, theta, px, py, mx, my, length, xp, yp, xm, ym
    , r_step, dstr, d2str, r_stretch;

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
  if (curve_ptr != NULL)
    PL_color( curve_ptr->curve_color );

/* draw the curve in solid */
  if (plot_mode & 1){
    if (curve_ptr == NULL){
      PL_move( 0.0, 0.5 );
      PL_draw( 1.0, 0.5 );
    }
    else if(curve_ptr != NULL && curve_ptr->curve_mapping != NULL){
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
  }

/* draw an arrow that indicates the direction of the stretching */
  if (curve_ptr == NULL){
    cp.x = 0.5; cp.y = 0.5;
    cp.xr = 1.0; cp.yr = 0.0; sigma = 1.0;
/* tangent */
    sx = sigma*cp.xr;
    sy = sigma*cp.yr;
    pi = 4.0*atan(1.0);
    theta = 20.0*pi/180.0;
    px = cos(theta)*sx - sin(theta)*sy;
    py = sin(theta)*sx + cos(theta)*sy;
    theta = -theta;
    mx = cos(theta)*sx - sin(theta)*sy;
    my = sin(theta)*sx + cos(theta)*sy;
    length = 0.05;
    xp = cp.x - length*px;
    yp = cp.y - length*py;
    xm = cp.x - length*mx;
    ym = cp.y - length*my;
    PL_move(xp, yp);
    PL_draw(cp.x, cp.y);
    PL_draw(xm, ym);
  }

/* draw tickmarks to indicate the parametrization */
  if (plot_mode & 4){
    if (curve_ptr == NULL){
      length = 0.01;
      r_step = 1.0/(grid_points - 1);
/* constant tangent */
      sx = 1.0;
      sy = 0.0;
      for (i=1; i<= grid_points; i++){
 	uniform_to_stretch( (i-1)*r_step, stretch_ptr, &r_stretch, &dstr, &d2str );
/* start and end points */
	xp = r_stretch - length*sy;
	yp = 0.5 + length*sx;
	xm = r_stretch + length*sy;
	ym = 0.5 - length*sx;
	PL_move(xp, yp);
	PL_draw(xm, ym);
      }
    }
    else if (curve_ptr != NULL && curve_ptr->curve_mapping != NULL){
      length = 0.01*real_max(curve_ptr->x_max-curve_ptr->x_min, 
			     curve_ptr->y_max-curve_ptr->y_min);
      r_step = 1.0/(grid_points-1);
      for (i=1; i<= grid_points; i++){
	uniform_to_stretch( (i-1)*r_step, stretch_ptr, &(cp.r), &dstr, &d2str );
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
  }

  if (plot_mode & 8){
    PL_title( 2, (curve_ptr == NULL)? "no_curve": curve_ptr->curve_name );
  }

  if (plot_mode & 16)
    plot_control_points( curve_ptr );

} /* end curve_plot */


generic_stretching *
read_stretching( int32 dir ){
  generic_stretching *stretch_ptr=NULL;
  int stretch_type;

  hget_int( &(stretch_type), "stretch_type", dir );

  if (stretch_type == 0)
    return NULL;
  else{
    stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );
    stretch_ptr->stretch_type = stretch_type;
  }

/* assign the appropriate curve-specific pointers */
  switch( stretch_ptr->stretch_type ){

  case 1:
/* arclength */
    read_arcl_stretch( dir, stretch_ptr );

    break;

  case 2:
/* exponential */
    read_exp_stretch( dir, stretch_ptr );
    break;

  case 3:
/* layer stretching */
    read_layer_stretch( dir, stretch_ptr );
    break;

  case 4:
/* hyperbolic tangent */
    read_tanh_stretch( dir, stretch_ptr );
    break;

  default:
    printf("Unknown stretch type: %i in read stretch\n", 
	   stretch_ptr->stretch_type);
    exit(1);
  }

  return stretch_ptr;
}

void 
write_stretch( int32 dir, generic_stretching *stretch_ptr ){

/* save the stretch type */
  hput_int( (stretch_ptr? stretch_ptr->stretch_type: 0), "stretch_type", dir );

/* only save stretching parameters if there is a stretching! */
  if (!stretch_ptr) return;

/* write the specific data (in the same directory to make it simple) */
  (*stretch_ptr->write_stretch_data)( dir, stretch_ptr );

}


generic_stretching *
delete_generic_stretching( generic_stretching *stretch_ptr ){
  if (stretch_ptr != NULL){
    (*stretch_ptr->cleanup_stretching)( stretch_ptr );
    free( stretch_ptr );
  }
  return NULL;
}

void 
uniform_to_stretch( real uniform, generic_stretching *stretch_ptr, 
			real *stretch, real *d_stretch, real *d2_stretch ){
/* *d_stretch will contain the derivataive of the stretching function with respect */
/* to the uniform parameter, and *d2_strech will contain the second derivative. */
  if ( stretch_ptr == NULL ){
    *stretch = uniform;
    *d_stretch = 1.0;
    *d2_stretch = 0.0;
  }
  else
    (*stretch_ptr->uniform_to_stretch)( uniform, stretch_ptr, stretch, 
				       d_stretch, d2_stretch );
}


generic_stretching *
copy_stretching( generic_stretching *old_stretch_ptr ){
  generic_stretching *stretch_ptr;

/* it is legal to make a copy of "no stretching" */
  if (old_stretch_ptr == NULL) return NULL;

  stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );

  stretch_ptr->stretch_type       = old_stretch_ptr->stretch_type;
  stretch_ptr->uniform_to_stretch = old_stretch_ptr->uniform_to_stretch;
  stretch_ptr->write_stretch_data = old_stretch_ptr->write_stretch_data;
  stretch_ptr->set_stretching     = old_stretch_ptr->set_stretching;
  stretch_ptr->cleanup_stretching = old_stretch_ptr->cleanup_stretching;
  stretch_ptr->copy_stretching    = old_stretch_ptr->copy_stretching;

/* copy the specific data */
  stretch_ptr->stretch_data_ptr = 
    (*old_stretch_ptr->copy_stretching)( old_stretch_ptr->stretch_data_ptr );

  return stretch_ptr;
}


int 
set_stretching(input_output *io_ptr, generic_stretching *stretch_ptr,
	       generic_curve *curve_ptr, int *grid_points ){
  if (stretch_ptr == NULL || stretch_ptr->set_stretching == NULL)
    return 0;
  else{
    (*stretch_ptr->set_stretching)( io_ptr, stretch_ptr, curve_ptr, grid_points );
    return 1;
  }
}
