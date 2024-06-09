#ifndef mappings_h
#define mappings_h

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>

/*#include <stupid_compiler.h>*/
#include <real.h>
#include <c_array.h>
#include <nrutil.h>
#include <spline.h>
#include <min_max.h>
#include <xplot.h>

#include "version.h"

#include "input_output.h"
#include "grid_point.h"

#include "curves.h"
#include "stretchings.h"

/* data structure for a mapping (component grid) */
typedef struct generic_mapping{

/* pointers to next and previous grids in the mapping list */
  struct generic_mapping *next_m,*prev_m; 

  int priority;             /* priority of this mapping in the composite grid */

  char *grid_name;          /* name of this component grid */
  int inverse_known;        /* inverse_known = 1 if the inverse of the mapping 
			       is defined in the grid function */

/* pointer to the function that defines the forward mapping */
  int (*grid_mapping)(grid_point *gp_ptr, void *mapping_data_ptr); 

/* pointer to the function that defines the inverse mapping 
   set to NULL if the inverse mapping is not known */
  int (*inverse_grid_mapping)(grid_point *gp_ptr, void *mapping_data_ptr);

/* pointer to the data necessary to compute the mapping */
  void *mapping_data_ptr;
/* pointer to the function that cleans the data for the mapping */
  void *(*cleanup_mapping_data)(void *mapping_data_ptr);

/* pointer to the function for interactive modification of an existing mapping
   of the present type */
  void (*set_mapping)( input_output *, struct generic_mapping *, generic_curve *first_c);

/* pointers to the functions for saving and recovering a mapping from a binary file */
  void (*write_specific_mapping)( int32 dir, void *mapping_data_ptr );

/* pointer to function for copying the specific data */
  void *(*copy_mapping_data)( void *mapping_data_ptr );

  int grid_type;            
/* special type of grid:
   1: Cartesian. The constant Jacobian will be saved in xr, xs, yr, ys. 
   2: Annular. (NOT USED ANYMORE)
   3. Half-annular. (NOT USED ANYMORE)
   4. Normal-curve mapping. Normals out from any known kind of curve.
   5. Linear interpolation mapping. Linear interpolation between two curves of any
   known kind.
   6. Interpolate bi-linearly between discrete grid points.
   7. Theodorsen-Garrick conformal mapping (suitable for wing sections)
   8. Hyperbolic mapping. Curved grid lines out from a curve. Computed by solving 
      a Hamliton-Jacobi equation.
*/

/* number of non-fictitious points in the r and s directions, respectively */
  int r_points,s_points;    
  
  int r_period,s_period;    /* r_period and s_period refer to the r and 
			       s-directions, respectively. They mean:
			       period=1: periodic
			       period=0: not periodic */
  int_array_2d *bc_ptr;      /* bc(is,id): boundary condition on side is=1,2
			       in direction id=1,2 */
  int_array_2d *curve_ptr;  /* curve(is,id): side is=1,2, direction 
			       id=1,2. curve(is,id)= a number identifying 
			       each side of the grid with a boundary curve. 
			       Used when sides of two or more grids share the
			       same boundary curve */

/* min and max x and y coordinates before rotation, translation and scaling*/
  real x_min_0, x_max_0, y_min_0, y_max_0; 

/* min and max x and y coordinates after rotation, translation and scaling */
  real x_min, x_max, y_min, y_max; 

/* plotting information */
  int grid_plot_mode;
  real xyab[2][2];
  color_type grid_color;

/* rotation, translation, scaling */
  real theta_rot, cos_theta, sin_theta, x_trans, y_trans, scaling;

/* external segments (in parameter space) on physical boundaries */
  real gap_low_r[2], gap_high_r[2], gap_low_s[2], gap_high_s[2];

} generic_mapping;

/* structure for grid quality */
typedef struct grid_quality{

  real max_ar; int i_ar, j_ar;

  real min_size_r; int i_msr, j_msr;
  real min_size_s; int i_mss, j_mss;

  real min_ortho; int i_ortho, j_ortho;
  real min_area; int i_area, j_area;

  real max_theta_r; int i_theta_r, j_theta_r;
  real max_change_r; int i_change_r, j_change_r;

  real max_theta_s; int i_theta_s, j_theta_s;
  real max_change_s; int i_change_s, j_change_s;

} grid_quality;

void 
mark_mapping_bndry(input_output *io_ptr, generic_mapping *map_ptr);
generic_mapping *
read_mapping( int32 dir, generic_curve *first_c );
void 
write_mapping( int32 dir, generic_mapping *grid_ptr );
generic_mapping *choose_mapping(input_output *io_ptr, char *name, 
				generic_curve *first_c);
void delete_generic_mapping( generic_mapping *grid_ptr );
void plot_generic_mapping(generic_mapping *grid_ptr, int draw_mode);
void 
set_boundary_condition(input_output *io_ptr, generic_mapping *grid_ptr);
void 
set_curve_label(input_output *io_ptr, generic_mapping *grid_ptr);
void 
set_gridlines(input_output *io_ptr, generic_mapping *grid_ptr);
void 
set_grid_plot_mode( input_output *io_ptr, generic_mapping *grid_ptr );
void 
show_mapping_parameters( generic_mapping *grid_ptr );
void 
reset_view_port( generic_mapping *grid_ptr );
generic_mapping *get_grid_ptr( input_output *io_ptr, generic_mapping *first_m );
void 
set_transformation(input_output *io_ptr, generic_mapping *map_ptr);
int 
forward_mapping( grid_point *gp_ptr, generic_mapping *map_ptr );
int 
backward_mapping( grid_point *gp_ptr, generic_mapping *map_ptr );
generic_mapping *
copy_mapping( char *name, generic_mapping *old_grid_ptr);
void 
map_transform( real x_old, real y_old, real *x_new, real *y_new, 
	      generic_mapping *map_ptr );
void 
inv_map_transform( real x_old, real y_old, real *x_new, real *y_new, 
		  generic_mapping *map_ptr );
void 
set_mapping( input_output *io_ptr, generic_mapping *map_ptr, generic_curve *first_c );
void 
get_grid_quality( grid_quality *quality, generic_mapping *map_ptr );

/* mappings known to xcog */
#include "cartesian_grid.h"
#include "normal_curve.h"
#include "linear_interp_grid.h"
#include "disc_point_mapping.h"
#include "tg_grid.h"
#include "hyperbolic.h"

/* define some macros for computing indices of multi-dimensional arrays */
#define bc(i,j)      compute_index_2d(grid_ptr->bc_ptr,i,j)
#define curve(is,id) compute_index_2d(grid_ptr->curve_ptr,is,id)

#endif
