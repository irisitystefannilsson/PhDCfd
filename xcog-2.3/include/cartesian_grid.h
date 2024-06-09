typedef struct cartesian_grid_info{
  real xmin, xmax, ymin, ymax;
  generic_stretching *r_stretch, *s_stretch;
}cartesian_grid_info;

/* function specifications */
cartesian_grid_info *
init_cartesian_mapping(generic_mapping *grid_ptr, real xmin, real xmax, 
		       real ymin, real ymax);
int
read_cartesian_mapping( int32 dir, generic_mapping *grid_ptr, generic_curve *first_c);
