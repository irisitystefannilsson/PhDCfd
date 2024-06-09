typedef struct arcl_stretch_info {
  int n;
  real dr_darcl_0, dr_darcl_1;
  real *arcl, *r, *rcoeff; /* independent variable is arcl (arclength) */
} arcl_stretch_info;

void 
init_arclength_stretch( generic_stretching *stretch_ptr, generic_curve *curve_ptr );
void 
read_arcl_stretch( int32 dir, generic_stretching *stretch_ptr );
