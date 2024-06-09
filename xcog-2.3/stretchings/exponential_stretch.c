#include "stretchings.h"

/* prototypes for private functions */
static void 
exp_uniform_to_stretch( real uniform, generic_stretching *stretch_ptr, 
		       real *stretch, real *d_stretch, 
		       real *d2_stretch );
static void 
cleanup_exp_stretch( generic_stretching *stretch_ptr );
static void 
write_exp_stretch( int32 dir, generic_stretching *stretch_ptr );
static void
set_exp_stretch( input_output *io_ptr, generic_stretching *stretch_ptr, 
		  generic_curve *curve_ptr, int *grid_points );
static void *
copy_exp_stretch( void *stretch_data_ptr );
static void 
compute_exp_stretch( exp_stretch_info *info );
/* end private prototypes */


void 
init_exp_stretch( generic_stretching *stretch_ptr, int grid_points ){
  exp_stretch_info *info;

/* allocate memory for the arcl_stretch_info structure */
  info = (exp_stretch_info *) malloc( sizeof(exp_stretch_info) );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 2;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = exp_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_exp_stretch;
  stretch_ptr->set_stretching = set_exp_stretch;
  stretch_ptr->cleanup_stretching = cleanup_exp_stretch;
  stretch_ptr->copy_stretching = copy_exp_stretch;

/* focus at r=0 */
  info->focus = 0;

/* give a default starting grid step */
  info->ds_dx_0 = 0.5;

/* compute alpha */
  compute_exp_stretch( info );
}

void 
read_exp_stretch( int32 dir, generic_stretching *stretch_ptr ){
  exp_stretch_info *info;

  info = (exp_stretch_info *) malloc( sizeof(exp_stretch_info) );

  hget_int( &(info->focus),    "focus",   dir );
  hget_real( &(info->alpha),   "alpha",   dir );
  hget_real( &(info->ds_dx_0), "ds_dx_0", dir );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 2;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = exp_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_exp_stretch;
  stretch_ptr->set_stretching = set_exp_stretch;
  stretch_ptr->cleanup_stretching = cleanup_exp_stretch;
  stretch_ptr->copy_stretching = copy_exp_stretch;

}

static void
set_exp_stretch( input_output *io_ptr, generic_stretching *stretch_ptr, 
		  generic_curve *curve_ptr, int *grid_points ){
  char *prompt;
  const real eps=1.e-7;
  int icom, quit, replot;
  const int level=3, no_indent=0;
  real xyab[2][2], sigma;
  exp_stretch_info *info;
  curve_point cp;

#include "set_exp_com.h"

  info = (exp_stretch_info *) stretch_ptr->stretch_data_ptr;

  prompt = "change exp stretching>";
  replot = 1;
  if (curve_ptr == NULL){
    xyab[0][0] = -0.05; xyab[1][0] = -0.05;
    xyab[0][1] =  1.05; xyab[1][1] =  1.05;
    PL_window(2, xyab);
  }
  else
    PL_window(2, curve_ptr->xyab );
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
	cp.r = (info->focus == 0)? 0.0 : 1.0;
	curve_function( &cp, curve_ptr );
	sigma = sqrt( cp.xr*cp.xr + cp.yr*cp.yr );
      }
      else
	sigma = 1.0;
/* change strength */
      info->ds_dx_0 = (*grid_points - 1) / sigma * 
	real_max( eps, get_real( io_ptr, "Starting grid size > 0:",
				info->ds_dx_0 * sigma/ (*grid_points - 1), 
				no_indent ) );
      compute_exp_stretch( info );
      break;

    case 1:
/* reverse parametrization */
      info->focus = !(info->focus);
      break;

    case 2:
/* change number of grid points */
      *grid_points = int_max(2, get_int(io_ptr, 
					"Enter number of gridlines >= 2: ", 
					*grid_points,
					no_indent));
      break;

    case 3:
/* default strength */
      info->ds_dx_0 = 0.5;
      compute_exp_stretch( info );
      break;

    case 4:
/* call the curve function to get the ratio between its parameter and the arclength */
      if (curve_ptr != NULL){
	cp.r = (info->focus == 0)? 0.0 : 1.0;
	curve_function( &cp, curve_ptr );
	sigma = sqrt( cp.xr*cp.xr + cp.yr*cp.yr );
      }
      else
	sigma = 1.0;
	
/* show stretching parameters */
      printf("Stretching parameters:\n");
      printf("The stretching is focused at %i\n", info->focus );
      printf("Number of grid points: %i\n", *grid_points );
      printf("Starting grid size: %f\n", info->ds_dx_0 * sigma / (*grid_points-1) );
      printf("Relative increase in grid size: %f.\n", 
	     info->alpha / (*grid_points - 1));
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
write_exp_stretch( int32 dir, generic_stretching *stretch_ptr ){
  exp_stretch_info *info;

  if (stretch_ptr == NULL) return;

  info = (exp_stretch_info *) stretch_ptr->stretch_data_ptr;

  hput_int( info->focus, "focus", dir );
  hput_real( info->alpha, "alpha", dir );
  hput_real( info->ds_dx_0, "ds_dx_0", dir );
}

static void 
exp_uniform_to_stretch(real uniform, generic_stretching *stretch_ptr, 
		       real *stretch, real *d_stretch, real *d2_stretch ){
  exp_stretch_info *info;
  real s0, s1, ds0, ds1;
/* the stretching function is extrapolated outside of u_min < uniform < u_max */
  const real u_min=-0.2, u_max=1.2; 

  info = (exp_stretch_info *) stretch_ptr->stretch_data_ptr;

/* flip the parameter? */
  if (info->focus == 1)
    uniform = 1.0 - uniform;

/* exponential stretching function */
  if (info->alpha > 0.0){
    if (uniform < u_min ){
      s0 = (exp(info->alpha * u_min) - 1.0) / (exp(info->alpha) - 1.0);
      ds0 = info->alpha * exp(info->alpha*u_min) / (exp(info->alpha) - 1.0);
      *stretch    = s0 + (uniform-u_min)*ds0;
      *d_stretch  = ds0;
      *d2_stretch = 0.0;
    }
    else if (uniform > u_max){
      s1 = (exp(info->alpha * u_max) - 1.0) / (exp(info->alpha) - 1.0);
      ds1 = info->alpha * exp(info->alpha*u_max) / (exp(info->alpha) - 1.0);
      *stretch    = s1 + (uniform-u_max)*ds1;
      *d_stretch  = ds1;
      *d2_stretch = 0.0;
    }
    else{
      *stretch = (exp( info->alpha * uniform ) - 1.0)/
	(exp( info->alpha ) - 1.0);
      *d_stretch = info->alpha * exp( info->alpha * uniform )/
	(exp( info->alpha ) - 1.0);
      *d2_stretch = info->alpha * (*d_stretch);
    }
  } /* end alpha > 0 */
  else{
    *stretch = uniform;
    *d_stretch = 1.0;
    *d2_stretch = 0.0;
  }

  if (info->focus == 1){
    *stretch = 1.0 - *stretch;
    *d2_stretch = - *d2_stretch;
  }
}

static void 
cleanup_exp_stretch( generic_stretching *stretch_ptr ){
  exp_stretch_info *info;

  if (stretch_ptr == NULL) return;

/* cast the pointer to the right type */
  info = (exp_stretch_info *) stretch_ptr->stretch_data_ptr;

  free( info );
}

static void *
copy_exp_stretch( void *stretch_data_ptr ){
  exp_stretch_info *info, *old_info;
  
  old_info = (exp_stretch_info *)stretch_data_ptr;

  info = (exp_stretch_info *) malloc( sizeof(exp_stretch_info) );

  info->focus = old_info->focus;
  info->alpha = old_info->alpha;

  return (void *)info;
}

static void 
compute_exp_stretch( exp_stretch_info *info ){
  real b, a, res;
  const real b_min = 1.0 + 1.e-4;
  int i;
  const int n_iter=10;

  b = 1.0 / info->ds_dx_0;

/* solve a: exp(a) - 1 = a * b */
/* this equation is only solvable if b > 1 */

  if (b < b_min){
    printf("Error in compute_exp_stretch: b < %f\n", b_min);
    b = b_min;
    info->ds_dx_0 = 1.0 / b;
  }

/* use Taylor series expansion to approximate exp(a)-1 = a + a^2/2 + a^3/6 + ... */
/* in order to get an initial guess */    
  a = -1.5 + sqrt( 1.5*1.5 + 6.0*( b - 1.0 ) );

/* newton iteration */
  for (i=1; i<=n_iter; i++){
    a = a - (exp(a) - 1.0 - b*a)/(exp(a) - b);
  }

/* check for reasonable convergence */
  if ((res=fabs( exp(a) - 1.0 - b*a )) > NEWTON_EPS)
    printf("Warning: Newton iteration in compute_exp_stretch converged "
	   "poorly.\nResidual after %i iterations: %e\n", n_iter, res);

  info->alpha = a;
}

