typedef struct disc_point_info{
  int nr, ns, order;
  int r_period, s_period;

/* x(i,j) and y(i,j) are the x and y-coordinates, respectively of grid point i,j. */
  real_array_2d *x_ptr,*y_ptr;

  generic_stretching *r_stretch, *s_stretch;
}disc_point_info;

#define x_disc(i,j) compute_index_2d(info->x_ptr,i,j)
#define y_disc(i,j) compute_index_2d(info->y_ptr,i,j)

/* function specifications */
disc_point_info *
init_disc_point_map( input_output *io_ptr, generic_mapping *grid_ptr );
int
read_disc_point_map( int32 dir, generic_mapping *grid_ptr, generic_curve *first_c);



