#ifndef cubic_spline_h
#define cubic_spline_h

#define FREE_BC 1.0e30
typedef struct cubic_spline_data{
  int n, periodic;
  real *r, *x, *y, *xcoeff, *ycoeff;
  real bc_x1, bc_xn, bc_y1, bc_yn, rmax;
} cubic_spline_data;

void *
init_cubic_spline( input_output *io_ptr, generic_curve *curve_ptr,
		  char *status );
void 
read_cubic_spline( int32 dir, generic_curve *curve_ptr );

#endif
