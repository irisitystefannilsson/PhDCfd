#ifndef smooth_poly_curve_h
#define smooth_poly_curve_h

typedef struct smooth_poly_data{
  int n;
  real *r, *x, *y, *dxdr, *dydr, *vx0, *vy0;
  real x1, y1, sharpness;
} smooth_poly_data;

void *
init_smooth_poly_curve( input_output *io_ptr, generic_curve *curve_ptr, char *status );
void 
read_smooth_poly( int32 dir, generic_curve *curve_ptr );

#endif
