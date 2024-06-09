typedef struct linear_interp_data{
  generic_curve *curve_s0_ptr, *curve_s1_ptr;
  int flip_s0, flip_s1;
  generic_stretching *curve_s0_stretch, *curve_s1_stretch, *normal_stretch;
} linear_interp_data;

linear_interp_data *
init_linear_interp(input_output *io_ptr, generic_mapping *grid_ptr,
		   generic_curve *curve_s0_ptr, generic_curve *curve_s1_ptr);
int
read_linear_interp( int32 dir, generic_mapping *grid_ptr, generic_curve *first_c );



