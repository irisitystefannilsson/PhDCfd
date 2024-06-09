typedef struct exp_stretch_info {
  int focus;
  real alpha;
  real ds_dx_0;
} exp_stretch_info;

void 
init_exp_stretch( generic_stretching *stretch_ptr, int grid_points );
void 
read_exp_stretch( int32 dir, generic_stretching *stretch_ptr );


