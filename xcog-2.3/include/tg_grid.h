#include <complex.h>
#include <theodorsen_garrick.h>

typedef struct tg_grid_info{
  real xmin, xmax, ymin, ymax, s_scaling;
  real_array_2d *coordinate_list;
  tgcoefflist tgcoeff;
  generic_stretching *r_stretch, *s_stretch;
}tg_grid_info;

/* function specifications */
tg_grid_info *
init_tg_mapping(input_output *io_ptr,generic_mapping *grid_ptr);
int 
read_tg_mapping( int32 dir, generic_mapping *grid_ptr, generic_curve *first_c);
