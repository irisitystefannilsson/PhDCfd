#include "stretchings.h"

/* prototypes for private functions */
static void arcl_uniform_to_stretch( real r_uniform, generic_stretching *stretch_ptr, 
				    real *r_stretch, real *dstr_duni, 
				    real *d2str_duni2);
static void compute_arcl_stretch_info( generic_curve *curve_ptr, 
				      arcl_stretch_info *info );
static void alloc_arcl_stretch_arrays( arcl_stretch_info *info );
static void free_arcl_stretch_arrays( arcl_stretch_info *info );
static void 
write_arcl_stretch( int32 dir, generic_stretching *stretch_ptr );
static void cleanup_arcl_stretching( generic_stretching *stretch_ptr );
static void *copy_arcl_stretch( void *stretch_data_ptr );
/* end private prototypes */

/* functions for arclength stretching */

void init_arclength_stretch( generic_stretching *stretch_ptr,
			    generic_curve *curve_ptr ){
  const int n_reparametrization = 50;
  arcl_stretch_info *info;

/* allocate memory for the arcl_stretch_info structure */
  info = (arcl_stretch_info *) malloc( sizeof(arcl_stretch_info) );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 1;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = arcl_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_arcl_stretch;
  stretch_ptr->set_stretching = NULL;
  stretch_ptr->cleanup_stretching = cleanup_arcl_stretching;
  stretch_ptr->copy_stretching = copy_arcl_stretch;

/* initialize the stretching */
  info->n = n_reparametrization;
  alloc_arcl_stretch_arrays( info );
  compute_arcl_stretch_info( curve_ptr, info );

}

void 
read_arcl_stretch( int32 dir, generic_stretching *stretch_ptr ){
  arcl_stretch_info *info;
  int low, high;

  info = (arcl_stretch_info *) malloc( sizeof(arcl_stretch_info) );

  hget_int(  &(info->n),          "n",          dir );
  hget_real( &(info->dr_darcl_0), "dr_darcl_0", dir );
  hget_real( &(info->dr_darcl_1), "dr_darcl_1", dir );

/* read in the vectors */
  info->arcl   = hget_vector( "arcl",   &low, &high, dir );
  info->r      = hget_vector( "r",      &low, &high, dir );
  info->rcoeff = hget_vector( "rcoeff", &low, &high, dir );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 1;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = arcl_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_arcl_stretch;
  stretch_ptr->set_stretching = NULL;
  stretch_ptr->cleanup_stretching = cleanup_arcl_stretching;
  stretch_ptr->copy_stretching = copy_arcl_stretch;

}

/* end public functions */

/* private functions start here */

static void 
arcl_uniform_to_stretch( real arcl, generic_stretching *stretch_ptr, 
			real *r_stretch, real *dstr_duni, real *d2str_duni2){
  arcl_stretch_info *info;

/* cast the input pointer to the right type */
  info = (arcl_stretch_info *) stretch_ptr->stretch_data_ptr;

/* evaluate the spline */
  cseval( info->arcl, info->r, info->rcoeff, info->n, arcl, r_stretch, 
	 dstr_duni, d2str_duni2 );

}


static void 
compute_arcl_stretch_info( generic_curve *curve_ptr, arcl_stretch_info *info ){
  int i;
  real dr, darcl_dr_0, darcl_dr_1, smax;
  curve_point cp;
/* assume the arrays already have been allocated */
  dr = 1.0/ (info->n - 1);
/* integrate the arclength - parameter relation */

/* start point */
  info->r[1] = 0.0;
  info->arcl[1] = 0.0;
  cp.r = 0.0;
  curve_function( &cp, curve_ptr );
  darcl_dr_1 = darcl_dr_0 = sqrt( cp.xr*cp.xr + cp.yr*cp.yr );
  info->dr_darcl_0 = 1.0/ darcl_dr_0;

  for (i=2; i<=info->n; i++){
    cp.r = (i-1)*dr;
    curve_function( &cp, curve_ptr );
    darcl_dr_1 = sqrt( cp.xr*cp.xr + cp.yr*cp.yr );
    info->r[i] = cp.r;
    info->arcl[i] = info->arcl[i-1] + 0.5 * dr * ( darcl_dr_0 + darcl_dr_1 );
    darcl_dr_0 = darcl_dr_1;
  }
/* save the derivative for the end point */
  info->dr_darcl_1 = 1.0/ darcl_dr_1;

/* scale the arclength */
  smax = info->arcl[info->n];
  for (i=2; i<=info->n; i++)
    info->arcl[i] = info->arcl[i]/smax;

/* scale the boundary condition */
  info->dr_darcl_0 = smax * info->dr_darcl_0;
  info->dr_darcl_1 = smax * info->dr_darcl_1;

/* compute the spline coefficients */
  spline( info->arcl, info->r, info->n, info->dr_darcl_0, info->dr_darcl_1,
	 info->rcoeff );
}


static void 
alloc_arcl_stretch_arrays( arcl_stretch_info *info ){
/* allocate space for all arrays in the stretch_info structure */
/* input: n = number of node points in the arclength integration. */
/* output: space for all arrays. */

  info->arcl   = vector( 1, info->n );
  info->r      = vector( 1, info->n );
  info->rcoeff = vector( 1, info->n );
}

static void
cleanup_arcl_stretching( generic_stretching *stretch_ptr ){
  arcl_stretch_info *info;

  if (stretch_ptr == NULL) return;

  info = (arcl_stretch_info *) stretch_ptr->stretch_data_ptr;

  free_arcl_stretch_arrays( info );
  free( info );
}

static void 
free_arcl_stretch_arrays( arcl_stretch_info *info ){
/* free the space for arrays in the stretch_info structure 
   input: n_tanh and all arrays
   output: none                */

  free_vector( info->arcl,   1, info->n );
  free_vector( info->r,      1, info->n );
  free_vector( info->rcoeff, 1, info->n );
}


static void 
write_arcl_stretch( int32 dir, generic_stretching *stretch_ptr ){
  arcl_stretch_info *info;

  if (stretch_ptr == NULL) return;

/* cast the generic pointer to the specific type */
  info = (arcl_stretch_info *) stretch_ptr->stretch_data_ptr;

  hput_int( info->n, "n", dir );
  hput_real( info->dr_darcl_0, "dr_darcl_0", dir );
  hput_real( info->dr_darcl_1, "dr_darcl_1", dir );

/* write in the vectors */
  hput_vector( info->arcl,   1, info->n, "arcl", dir );
  hput_vector( info->r,      1, info->n, "r", dir );
  hput_vector( info->rcoeff, 1, info->n, "rcoeff", dir );
}


static void *
copy_arcl_stretch( void *stretch_data_ptr ){
  arcl_stretch_info *info, *old_info;
  int i;
  
  old_info = (arcl_stretch_info *)stretch_data_ptr;

  info = (arcl_stretch_info *) malloc( sizeof(arcl_stretch_info) );

  info->n = old_info->n;
  info->dr_darcl_0 = old_info->dr_darcl_0;
  info->dr_darcl_1 = old_info->dr_darcl_1;

/* allocate space for the vectors */
  alloc_arcl_stretch_arrays( info );

/* copy the array elements */
  for( i=1; i<=info->n; i++){
    info->arcl[i]   = old_info->arcl[i];
    info->r[i]      = old_info->r[i];
    info->rcoeff[i] = old_info->rcoeff[i];
  }

  return (void *)info;
}
