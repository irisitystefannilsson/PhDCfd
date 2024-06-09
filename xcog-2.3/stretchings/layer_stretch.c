#include "stretchings.h"

/* functions for a layered stretching */

/* private member functions */
static void layer_stretch_to_uniform( real r_stretch, layer_stretch_info *info, 
				     real *r_uniform, real *duni_dstr, 
				     real *d2uni_dstr2 );
static void set_layer_stretch( input_output *io_ptr, 
			      generic_stretching *stretch_ptr, 
			      generic_curve *curve_ptr, 
			      int *grid_points );
static void layer_uniform_to_stretch( real r_uniform, 
				     generic_stretching *stretch_ptr, 
				     real *r_stretch, real *dstr_duni, 
				     real *d2str_duni2);
static void compute_layer_stretch_info( layer_stretch_info *info );
static void alloc_layer_stretch_arrays( layer_stretch_info *info );
static void free_layer_stretch_arrays( layer_stretch_info *info );
static void 
write_layer_stretch( int32 dir, generic_stretching *stretch_ptr );
static void cleanup_layer_stretching( generic_stretching *stretch_ptr );
static void default_strength_width( layer_stretch_info *info );
static void *copy_layer_stretch( void *stretch_data_ptr );
/* end */


void init_layer_stretch( generic_stretching *stretch_ptr, 
			generic_curve *curve_ptr ){
  layer_stretch_info *info;

/* allocate memory for the arcl_stretch_info structure */
  info = (layer_stretch_info *) malloc( sizeof(layer_stretch_info) );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 3;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = layer_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_layer_stretch;
  stretch_ptr->set_stretching = set_layer_stretch;
  stretch_ptr->cleanup_stretching = cleanup_layer_stretching;
  stretch_ptr->copy_stretching = copy_layer_stretch;

  if (curve_ptr != NULL){
/* get a copy of the parameter values at the node points of the curve */
    info->focus = copy_nodes( curve_ptr, &(info->n_tanh) );
  }
  else{
/* default node points are located at the ends */
    info->n_tanh = 2;
    info->focus = vector( 1, info->n_tanh );
    info->focus[1] = 0.0;
    info->focus[2] = 1.0;
  }

/* allocate the remaining vectors in the layer stretching */
  alloc_layer_stretch_arrays( info );

/* initialize the strength and width vectors */
  default_strength_width( info );

/* setup the stretching parameters */
  compute_layer_stretch_info( info );
}

void 
read_layer_stretch( int32 dir, generic_stretching *stretch_ptr ){
  layer_stretch_info *info;
  int low, high;

  info = (layer_stretch_info *) malloc( sizeof(layer_stretch_info) );

  hget_int(  &(info->n_tanh),  "n_tanh",  dir );
  hget_int(  &(info->n_nodes), "n_nodes", dir );
  hget_real( &(info->scale),   "scale",   dir );

/* read in the vectors */
  info->strength  = hget_vector( "strength",  &low, &high, dir );
  info->str_sharp = hget_vector( "str_sharp", &low, &high, dir );
  info->focus     = hget_vector( "focus",     &low, &high, dir );
  info->u0        = hget_vector( "u0",        &low, &high, dir );
  info->str_node  = hget_vector( "str_node",  &low, &high, dir );
  info->uni_node  = hget_vector( "uni_node",  &low, &high, dir );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 3;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = layer_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_layer_stretch;
  stretch_ptr->set_stretching = set_layer_stretch;
  stretch_ptr->cleanup_stretching = cleanup_layer_stretching;
  stretch_ptr->copy_stretching = copy_layer_stretch;

}

static void 
default_strength_width( layer_stretch_info *info ){
  int i;
  real quarter_dist;

  for (i=1; i<=info->n_tanh; i++){

/* default strength */
    if (info->n_tanh == 2)
      info->strength[i] = 0.0;
    else if (i==1 || i==info->n_tanh)
      info->strength[i] = 0.0;
    else
      info->strength[i] = 0.0; /*0.5/(info->n_tanh - 2);*/

/* default values for the width of the layer */
    if (i==1)
      quarter_dist = 0.25*( info->focus[2] - info->focus[1] );
    else if (i==info->n_tanh)
      quarter_dist = 0.25*( info->focus[info->n_tanh] - 
			   info->focus[info->n_tanh-1] );
    else
      quarter_dist = 0.25* real_min(info->focus[i+1] - info->focus[i],
				    info->focus[i] - info->focus[i-1]);

    info->str_sharp[i] = 2.0/quarter_dist;
  }
}

static void
set_layer_stretch( input_output *io_ptr, generic_stretching *stretch_ptr, 
		  generic_curve *curve_ptr, int *grid_points ){
  char question[120];
  char *prompt;
  const real eps=1.e-7;
  int level=3, icom, quit, i, no_indent=0, replot;
  real width, xyab[2][2], d_stretch, dummy;
  layer_stretch_info *info;

#include "set_layer_com.h"

  info = (layer_stretch_info *) stretch_ptr->stretch_data_ptr;

  prompt = "change layer stretching>";
  replot = 1;
  if (curve_ptr == NULL){
    xyab[0][0] = -0.05; xyab[1][0] = -0.05;
    xyab[0][1] =  1.05; xyab[1][1] =  1.05;
    PL_window(2, xyab);
  }
  else{
    PL_window(2, curve_ptr->xyab );
  }
  do{

    if (replot){
/* plot the mapping */
      PL_erase(2);
      PL_start_plot(2);
/* no arrow! It is only confusing when the parametrization is reversed by */
/* the mapping. */
      plot_stretch( stretch_ptr, curve_ptr, *grid_points, 29 );
      PL_end_plot();
    }

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    quit = 0;
    replot = 1;
    switch (icom) {

    case 0:
/* change strength */
      i = get_int(io_ptr, "Which corner: ", 2, no_indent );
      if (i<1 || i>info->n_tanh){
 	printf("That node does not exist!\n");
      }
      else{
 	sprintf( question, "Relative # of grid points clustered around "
		"node %i: ", i );
 	info->strength[i] = 
 	  real_max(0.0, get_real( io_ptr, question, 
 				 info->strength[i], level) );
/* setup the stretching parameters */
	compute_layer_stretch_info( info );
      }
      break;

    case 1:
/* change width */
      i = get_int(io_ptr, "Which corner: ", 2, no_indent );
      if (i<1 || i>info->n_tanh){
 	printf("That node does not exist!\n");
      }
      else{
 	sprintf( question, "Relative width of layer around corner %i: ", i);
 	width = get_real( io_ptr, question, 2.0/info->str_sharp[i], level);
 	if (width < eps)
 	  width = eps;
 	else if (width >0.99)
 	  width = 0.99;
 	info->str_sharp[i] = 2.0/width;
/* setup the stretching parameters */
	compute_layer_stretch_info( info );
      }
      break;

    case 2:
/* default strength and width */
      default_strength_width( info );
/* setup the stretching parameters */
      compute_layer_stretch_info( info );
      break;

    case 3:
/* change number of grid points */
      *grid_points = int_max(2, get_int(io_ptr, 
					"Enter number of gridlines >= 2: ", 
					*grid_points,
					no_indent));
      break;

    case 4:
/* show stretching parameters */
      printf("Stretching parameters:\n");
      printf("Number of grid points: %i\n", *grid_points );
      for (i=1; i<= info->n_tanh; i++){
	layer_uniform_to_stretch( info->focus[i], stretch_ptr, &dummy, 
				 &d_stretch, &dummy );
	printf("Node point %i: Grid size: %e, Strength: %e, Width: %e\n", 
	       i, d_stretch/(*grid_points - 1), info->strength[i],
	       2.0/info->str_sharp[i]);
      }
      printf("\n");
      replot = 0;
      break;

    case 5:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
			 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case 6:
/* exit */
      quit = 1;
      break;

    default:
      ;

    }
  }
  while( !quit );

}


static void 
layer_uniform_to_stretch( real r_uniform, generic_stretching *stretch_ptr, 
			 real *r_stretch, real *dstr_duni, real *d2str_duni2){
/*
   inverse stretching function
*/
  int i;
  real r_stretch0, r_str, r_uni, duni_dstr, d2uni_dstr2, str_low, str_high, uni_low
    , uni_high;
  layer_stretch_info *info;

/* cast the input pointer to the right type */
  info = (layer_stretch_info *) stretch_ptr->stretch_data_ptr;
/* find the right interval */
  for ( i=1; i <= info->n_nodes && r_uniform > info->uni_node[i]; i++);
/* now 1 <= i <= n+1 */

/* first interval */
  if ( i==1 ){
/* between 0.0 and node 1 */
    r_stretch0 = info->str_node[i] * r_uniform / info->uni_node[i];
  }
  else if ( i > info->n_nodes ){
/* between node n and 1.0 */
    r_stretch0 = (1.0 * (r_uniform - info->uni_node[i-1]) -
      info->str_node[i-1] * (r_uniform - 1.0))/
	(1.0 - info->uni_node[i-1]);
  }
  else
/* between nodes i and i-1 */
    r_stretch0 = (info->str_node[i] * (r_uniform - info->uni_node[i-1]) -
      info->str_node[i-1] * (r_uniform - info->uni_node[i]))/
	(info->uni_node[i] - info->uni_node[i-1]);

/* do a few newton iterations to get the exact solution */
  r_str = r_stretch0;
  for (i=1; i<=5; i++){
    layer_stretch_to_uniform( r_str, info, &r_uni, &duni_dstr, &d2uni_dstr2 );
    r_str = r_str - ( r_uni - r_uniform ) / duni_dstr;
    if (fabs( r_uni - r_uniform ) <= NEWTON_EPS) break;
  }
/* save result of newton iteration */
  *r_stretch = r_str;
  *dstr_duni = 1.0/ duni_dstr;
  *d2str_duni2 = - d2uni_dstr2 / duni_dstr / duni_dstr / duni_dstr;

/* check the convergence */
  if ( fabs( r_uni - r_uniform ) > NEWTON_EPS){
/*    printf("Newton failed at r_uniform = %e\n", r_uniform);*/
/* use the fail safe `halve the interval method' (intervall-halvering) */
/* we know that r_uni is a monotonously increasing function of r_str */
/* initial guess */
    r_str = r_stretch0;
    layer_stretch_to_uniform( r_str, info, &r_uni, &duni_dstr, &d2uni_dstr2 ); 
    if (r_uni < r_uniform){
      str_low = r_str;
      uni_low = r_uni;
      for (r_str = str_low; r_str <= 1.0 && r_uni < r_uniform; r_str += 0.01)
	layer_stretch_to_uniform( r_str, info, &r_uni, &duni_dstr, &d2uni_dstr2 );
      str_high = r_str;
      uni_high = r_uni;
    }
    else{
      str_high = r_str;
      uni_high = r_uni;
      for (r_str = str_high; r_str >= 0.0 && r_uni > r_uniform; r_str += -0.01)
	layer_stretch_to_uniform( r_str, info, &r_uni, &duni_dstr, &d2uni_dstr2 );
      str_low = r_str;
      uni_low = r_uni;
    }
/* subtract off r_uniform */
    uni_low += -r_uniform;
    uni_high += -r_uniform;

/* check the initial interval */
    if (uni_high * uni_low > 0.0){
      printf("Bad initial guess for the interval halving\n");
      return;
    }

/* the root finding algorithm */
    do{
      r_str = 0.5*( str_high + str_low );
      layer_stretch_to_uniform( r_str, info, &r_uni, &duni_dstr, &d2uni_dstr2 ); 
      r_uni += -r_uniform;
      if (fabs(r_uni) == 0.0)
	break;
      else if (uni_low * r_uni > 0.0){
	uni_low = r_uni;
	str_low = r_str;
      }
      else if (uni_high * r_uni > 0.0){
	uni_high = r_uni;
	str_high = r_str;
      }
    }
    while (fabs(uni_high - uni_low) > NEWTON_EPS);

/* save the result of the iteration */
    *r_stretch = 0.5 * ( str_high + str_low );
    *dstr_duni = 1.0 / duni_dstr;
    *d2str_duni2 = - d2uni_dstr2 / duni_dstr / duni_dstr / duni_dstr;
  }

}


static void 
layer_stretch_to_uniform( real r_stretch, layer_stretch_info *info, real *r_uniform, 
			real *duni_dstr, real *d2uni_dstr2 ){
  int i; real f, fp;
/* 
   stretching function 
*/
  *r_uniform = r_stretch;
  for (i=1; i<=info->n_tanh; i++)
    *r_uniform += 0.5 * info->strength[i] *
      tanh(info->str_sharp[i] * (r_stretch - info->focus[i]) ) - info->u0[i];
/* 
   first and second derivative 
*/
  *duni_dstr = 1.0;
  *d2uni_dstr2 = 0.0;
  for (i=1; i<=info->n_tanh; i++){
    f = 0.5 * info->strength[i] * info->str_sharp[i] /
      cosh(info->str_sharp[i] * (r_stretch - info->focus[i]) ) /
	cosh(info->str_sharp[i] * (r_stretch - info->focus[i]) );
    fp = -2.0 * f * info->str_sharp[i] *
      tanh( info->str_sharp[i] * (r_stretch - info->focus[i]) );
    *duni_dstr += f;
    *d2uni_dstr2 += fp;
  }

/* scaling factor */
  *r_uniform   = info->scale * *r_uniform;
  *duni_dstr   = info->scale * *duni_dstr;
  *d2uni_dstr2 = info->scale * *d2uni_dstr2;
}


static void 
compute_layer_stretch_info( layer_stretch_info *info ){
  real r_uniform, duni_dstr, d2uni_dstr2;
  int i;
/* input: n_tanh, focus, str_sharp, strength
   output: u0, scale */

/* evaluate u0 */
  for (i=1; i<=info->n_tanh; i++)
    info->u0[i] = 0.5 * info->strength[i] * 
      tanh( -(info->str_sharp[i] * info->focus[i]) );

/* compute the scale */
  info->scale = 1.0;
  layer_stretch_to_uniform( 1.0, info, &r_uniform, &duni_dstr, &d2uni_dstr2 );
  info->scale = 1.0/r_uniform;

/* compute the piecewise linear approximation to the stretching function */
  for (i=1; i<=info->n_tanh; i++){
    layer_stretch_to_uniform( info->focus[i], info, &r_uniform, &duni_dstr, 
			     &d2uni_dstr2 );
    info->str_node[2*i-1] = info->focus[i] - 1.0/info->str_sharp[i];
    info->uni_node[2*i-1] = r_uniform - duni_dstr/info->str_sharp[i];
    info->str_node[2*i]   = info->focus[i] + 1.0/info->str_sharp[i];
    info->uni_node[2*i]   = r_uniform + duni_dstr/info->str_sharp[i];
  }
}


static void 
alloc_layer_stretch_arrays( layer_stretch_info *info ){
/* allocate space for all arrays in the stretch_info structure 
   input: n_tanh
   output: space for all arrays */

  info->n_nodes = 2*info->n_tanh;

  info->strength  = vector( 1, info->n_tanh );
  info->str_sharp = vector( 1, info->n_tanh );
/* focus is allocated by `copy_nodes' */
/*  info->focus     = vector( 1, info->n_tanh );*/
  info->u0        = vector( 1, info->n_tanh );
  info->str_node  = vector( 1, 2*info->n_tanh);
  info->uni_node  = vector( 1, 2*info->n_tanh);
}

static void 
free_layer_stretch_arrays( layer_stretch_info *info ){
/* free the space for arrays in the stretch_info structure 
   input: n_tanh and all arrays
   output: none                */

  free_vector( info->strength,  1, info->n_tanh );
  free_vector( info->str_sharp, 1, info->n_tanh );
  free_vector( info->focus,     1, info->n_tanh );
  free_vector( info->u0,        1, info->n_tanh );
  free_vector( info->str_node,  1, 2*info->n_tanh );
  free_vector( info->uni_node,  1, 2*info->n_tanh );
}


static void 
write_layer_stretch( int32 dir, generic_stretching *stretch_ptr ){
  layer_stretch_info *info;

  if (stretch_ptr == NULL) return;

/* cast the generic pointer to the specific type */
  info = (layer_stretch_info *) stretch_ptr->stretch_data_ptr;

  hput_int(  info->n_tanh,  "n_tanh",  dir );
  hput_int(  info->n_nodes, "n_nodes", dir );
  hput_real( info->scale,   "scale",   dir );

/* write the vectors */
  hput_vector( info->strength,  1, info->n_tanh,   "strength",  dir  );
  hput_vector( info->str_sharp, 1, info->n_tanh,   "str_sharp", dir  );
  hput_vector( info->focus,     1, info->n_tanh,   "focus",     dir  );
  hput_vector( info->u0,        1, info->n_tanh,   "u0",        dir  );
  hput_vector( info->str_node,  1, 2*info->n_tanh, "str_node",  dir );
  hput_vector( info->uni_node,  1, 2*info->n_tanh, "uni_node",  dir );
}

static void
cleanup_layer_stretching( generic_stretching *stretch_ptr ){
  layer_stretch_info *info;

  if (stretch_ptr == NULL) return;

  info = (layer_stretch_info *) stretch_ptr->stretch_data_ptr;

  free_layer_stretch_arrays( info );
  free( info );
}

static void *copy_layer_stretch( void *stretch_data_ptr ){
  layer_stretch_info *info, *old_info;
  int i;
  
  old_info = (layer_stretch_info *)stretch_data_ptr;

  info = (layer_stretch_info *) malloc( sizeof(layer_stretch_info) );

  info->n_tanh = old_info->n_tanh;
  info->n_nodes = old_info->n_nodes;
  info->scale = old_info->scale;

/* allocate space for the vectors */
  alloc_layer_stretch_arrays( info );
  info->focus = vector( 1, info->n_tanh );

/* copy the array elements */
  for( i=1; i<=info->n_tanh; i++){
    info->strength[i]  = old_info->strength[i];
    info->str_sharp[i] = old_info->str_sharp[i];
    info->focus[i]     = old_info->focus[i];
    info->u0[i]        = old_info->u0[i];
  }
  for( i=1; i<=info->n_nodes; i++){
    info->str_node[i] = old_info->str_node[i];
    info->uni_node[i] = old_info->uni_node[i];
  }

  return (void *)info;
}
