typedef struct tanh_stretch_info {
  real a, delta, ds_dx_0, ds_dx_1;
} tanh_stretch_info;

void 
init_tanh_stretch( generic_stretching *stretch_ptr, int grid_points );
void 
read_tanh_stretch( int32 dir, generic_stretching *stretch_ptr );


