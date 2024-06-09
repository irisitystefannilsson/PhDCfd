#ifndef hyp_grid_h
#define hyp_grid_h

#include "stretchings.h"

typedef struct{
  generic_curve *curve_ptr;
  real width, curv_factor, v_min, smooth_factor, time_factor;
  int n1, n2, n_curve, n_cut; /* total # grid points n1 = n_curve + 2*n_cut */
  int flip, both_sides, periodic, c_grid;
  real d1, d2, x_branch, y_branch;
  real_array_2d *x_, *y_;

/* boundary conditions */
  int_array_1d *hyp_bc_;

/* flag to indicate if the mapping is up to date (=1), needs to be (re-)computed (=0), */
/* or if the computation was unsuccessful (=-1). */
  int status;

/* pointer to the stretching structures */
  generic_stretching *tangential_stretch, *normal_stretch;
} hyp_grid_data;

/* prototypes */
hyp_grid_data *
init_hyp_grid(input_output *io_ptr, generic_mapping *grid_ptr, generic_curve *curve_ptr);
int
read_hyp_map(int32 dir, generic_mapping *grid_ptr, generic_curve *curve_ptr);


/* macros for evaluating the grid point coordinates */
#define x_hyp(i,j)     compute_index_2d(info->x_,i,j)
#define y_hyp(i,j)     compute_index_2d(info->y_,i,j)
#define hyp_bc(i)      compute_index_1d(info->hyp_bc_,i)

#endif
