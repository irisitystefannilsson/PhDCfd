#include "stretchings.h"

/* prototypes for private functions */
static void 
tanh_uniform_to_stretch(real uniform, 
			generic_stretching *stretch_ptr, 
			real *stretch, real *d_stretch, 
			real *d2_stretch );
static void 
cleanup_tanh_stretch( generic_stretching *stretch_ptr );
static void 
write_tanh_stretch( int32 dir, generic_stretching *stretch_ptr );
static void 
set_tanh_stretch(input_output *io_ptr, generic_stretching *stretch_ptr, 
		 generic_curve *curve_ptr, int *grid_points );
static void *
copy_tanh_stretch( void *stretch_data_ptr );
static void 
compute_tanh_stretch( tanh_stretch_info *info );
/* end private prototypes */


void 
init_tanh_stretch( generic_stretching *stretch_ptr, int grid_points ){
  tanh_stretch_info *info;

/* allocate memory for the arcl_stretch_info structure */
  info = (tanh_stretch_info *) malloc( sizeof(tanh_stretch_info) );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 4;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = tanh_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_tanh_stretch;
  stretch_ptr->set_stretching = set_tanh_stretch;
  stretch_ptr->cleanup_stretching = cleanup_tanh_stretch;
  stretch_ptr->copy_stretching = copy_tanh_stretch;

/* default stretching is very mild */
  info->ds_dx_0 = 0.5;
  info->ds_dx_1 = 0.5;

/* compute the coefficients */
  compute_tanh_stretch( info );
}

void 
read_tanh_stretch( int32 dir, generic_stretching *stretch_ptr ){
  tanh_stretch_info *info;

  info = (tanh_stretch_info *) malloc( sizeof(tanh_stretch_info) );

  hget_real( &(info->a),       "a",       dir );
  hget_real( &(info->delta),   "delta",   dir );
  hget_real( &(info->ds_dx_0), "ds_dx_0", dir );
  hget_real( &(info->ds_dx_1), "ds_dx_1", dir );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 4;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = tanh_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_tanh_stretch;
  stretch_ptr->set_stretching = set_tanh_stretch;
  stretch_ptr->cleanup_stretching = cleanup_tanh_stretch;
  stretch_ptr->copy_stretching = copy_tanh_stretch;

}

static void
set_tanh_stretch( input_output *io_ptr, generic_stretching *stretch_ptr, 
		  generic_curve *curve_ptr, int *grid_points ){
  char *prompt;
  const real eps=1.e-7;
  int icom, quit, replot;
  const int level=3, no_indent=0;
  real xyab[2][2], sigma0, sigma1;
  tanh_stretch_info *info;
  curve_point cp;

#include "set_tanh_com.h"

  info = (tanh_stretch_info *) stretch_ptr->stretch_data_ptr;

  prompt = "change tanh stretching>";
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
/* call the curve function to get the ratio between its parameter and the arclength */
      if (curve_ptr != NULL){
	cp.r = 0.0;
	curve_function( &cp, curve_ptr );
	sigma0 = sqrt( cp.xr*cp.xr + cp.yr*cp.yr );
      }
      else
	sigma0 = 1.0;

/* start grid size */
      info->ds_dx_0 = (*grid_points - 1) / sigma0 * 
	real_max( eps, get_real( io_ptr, "Starting grid size > 0:",
				info->ds_dx_0 * sigma0 / (*grid_points - 1), 
				no_indent ) );
      compute_tanh_stretch( info );
      break;

    case 1:
/* call the curve function to get the ratio between its parameter and the arclength */
      if (curve_ptr != NULL){
	cp.r = 1.0;
	curve_function( &cp, curve_ptr );
	sigma1 = sqrt( cp.xr*cp.xr + cp.yr*cp.yr );
      }
      else
	sigma1 = 1.0;

/* ending grid size */
      info->ds_dx_1 = (*grid_points - 1) / sigma1 * 
	real_max( eps, get_real( io_ptr, "Ending grid size > 0:",
				info->ds_dx_1 * sigma1 / (*grid_points - 1), 
				no_indent ) );
      compute_tanh_stretch( info );
      break;

    case 2:
/* default strength */
      info->ds_dx_0 = info->ds_dx_1 = 0.5;
      compute_tanh_stretch( info );
      break;

    case 3:
/* change number of grid points */
      *grid_points = int_max(2, get_int(io_ptr, 
					"Enter number of gridlines >= 2: ", 
					*grid_points,
					no_indent));
      break;

    case 4:
/* call the curve function to get the ratio between its parameter and the arclength */
      if (curve_ptr != NULL){
	cp.r = 0.0;
	curve_function( &cp, curve_ptr );
	sigma0 = sqrt( cp.xr*cp.xr + cp.yr*cp.yr );
	cp.r = 1.0;
	curve_function( &cp, curve_ptr );
	sigma1 = sqrt( cp.xr*cp.xr + cp.yr*cp.yr );
      }
      else
	sigma0 = sigma1 = 1.0;

/* show stretching parameters */
      printf("Stretching parameters:\n");
      printf("Number of grid points: %i\n", *grid_points );
      printf("Starting grid size: %f\n", info->ds_dx_0 * sigma0 / (*grid_points - 1));
      printf("Ending grid size: %f\n", info->ds_dx_1 * sigma1 / (*grid_points - 1));
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
write_tanh_stretch( int32 dir, generic_stretching *stretch_ptr ){
  tanh_stretch_info *info;

  if (stretch_ptr == NULL) return;

  info = (tanh_stretch_info *) stretch_ptr->stretch_data_ptr;

  hput_real( info->a,       "a",       dir );
  hput_real( info->delta,   "delta",   dir );
  hput_real( info->ds_dx_0, "ds_dx_0", dir );
  hput_real( info->ds_dx_1, "ds_dx_1", dir );
}


static void 
tanh_uniform_to_stretch(real uniform, generic_stretching *stretch_ptr, 
			real *stretch, real *d_stretch, real *d2_stretch ){
  tanh_stretch_info *info;
  real ds0, ds1;
  const real u_min = -0.2, u_max = 1.2;

  info = (tanh_stretch_info *) stretch_ptr->stretch_data_ptr;

/* help function w */
#define W(U) 0.5 * ( 1.0 + tanh( info->delta*( U - 0.5 ) ) / tanh( 0.5*info->delta ) )
#define D_W(U) 0.5 * info->delta / tanh( 0.5*info->delta ) / \
  cosh( info->delta*( U - 0.5 ) ) / cosh( info->delta*( U - 0.5 ) )
#define D2_W(U) -info->delta * info->delta * tanh( info->delta*( U - 0.5 ) )/ \
  tanh( 0.5*info->delta ) / cosh( info->delta*( U - 0.5 ) ) / \
  cosh( info->delta*( U - 0.5 ) )

/* hyperbolic tangent stretching function */
#define STR(U) W(U) / ( info->a + (1.0 - info->a) * W(U) )
#define D_STR(U) info->a * D_W(U) / (info->a + (1.0-info->a)*W(U)) / \
  (info->a + (1.0-info->a)*W(U))
#define D2_STR(U) info->a * D2_W(U) / (info->a + (1.0-info->a)*W(U)) / \
  (info->a + (1.0-info->a)*W(U)) - \
  2.0 * info->a * (1.0 - info->a) * D_W(U) * D_W(U) / (info->a + (1.0-info->a)*W(U)) / \
  (info->a + (1.0-info->a)*W(U)) / (info->a + (1.0-info->a)*W(U))

  if ( uniform < u_min ){
    ds0 = D_STR(u_min);
    *stretch    = STR(u_min) + (uniform-u_min)*ds0;
    *d_stretch  = ds0;
    *d2_stretch = 0.0;
  }
  else if ( uniform > u_max ){
    ds1 = D_STR(u_max);
    *stretch    = STR(u_max) + (uniform-u_max)*ds1;
    *d_stretch  = ds1;
    *d2_stretch = 0.0;
  }
  else{
    *stretch    = STR(uniform);
    *d_stretch  = D_STR(uniform);
    *d2_stretch = D2_STR(uniform);
  }
#undef W
#undef D_W
#undef D2_W

#undef STR
#undef D_STR
#undef D2_STR
}

static void 
compute_tanh_stretch( tanh_stretch_info *info ){
  real b, d, res;
  const real b_min = 1.0 + 1.e-4;
  int i;
  const int n_iter=10;

  info->a = sqrt( info->ds_dx_1 / info->ds_dx_0 );

  b = 1.0 / sqrt( info->ds_dx_1 * info->ds_dx_0 );

/* solve sinh( d ) = d * b */
/* this equation is only solvable if b >= 1 */

  if (b < b_min){
    printf("Error in compute_tanh_stretch: b < %f\n", b_min);
    b = b_min;
    info->ds_dx_1 = 1.0 / info->ds_dx_0 / b / b;
  }

/* use Taylor series expansion to approximate sinh d = d + d^3/6 + ... */
/* in order to get an initial guess */    
  d = sqrt( 6.0*( b - 1.0 ) );

/* newton iteration */
  for (i=1; i<=n_iter; i++){
    d = d - (sinh(d) - b*d)/(cosh(d) - b);
  }

/* check for reasonable convergence */
  if ((res=fabs( sinh(d) - b*d )) > NEWTON_EPS)
    printf("Warning: Newton iteration in compute_tanh_stretch converged "
	   "poorly.\nResidual after %i iterations: %e\n", n_iter, res);

  info->delta = d;
}

static void 
cleanup_tanh_stretch( generic_stretching *stretch_ptr ){
  tanh_stretch_info *info;

  if (stretch_ptr == NULL) return;

/* cast the pointer to the right type */
  info = (tanh_stretch_info *) stretch_ptr->stretch_data_ptr;

  free( info );
}

static void *
copy_tanh_stretch( void *stretch_data_ptr ){
  tanh_stretch_info *info, *old_info;
  
  old_info = (tanh_stretch_info *)stretch_data_ptr;

  info = (tanh_stretch_info *) malloc( sizeof(tanh_stretch_info) );

  info->a = old_info->a;
  info->delta = old_info->delta;
  info->ds_dx_0 = old_info->ds_dx_0;
  info->ds_dx_1 = old_info->ds_dx_1;

  return (void *)info;
}
