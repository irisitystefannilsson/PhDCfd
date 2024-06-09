typedef struct normal_curve_data{
  generic_curve *curve_ptr;
  int is_constant;
  int flip;
  int n;
  real width_trans, width1;
  real *r, *dndr, *vn0, *un0, *width_l, *width_r, *width_jmp;

/* pointer to normal and tangential stretching structures */
  generic_stretching *tangential_stretch, *normal_stretch;

} normal_curve_data;

normal_curve_data *
init_normal_curve(input_output *io_ptr, generic_mapping *grid_ptr,
		  generic_curve *new_curve_ptr, char *status);
int
read_normal_curve( int32 dir, generic_mapping *grid_ptr, generic_curve *curve_ptr);
