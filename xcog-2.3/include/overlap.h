#ifndef overlap_h
#define overlap_h

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
#include <linked_list.h>

#include "version.h"

#include "input_output.h"
#include "grid_point.h"

/* typedefs */

/* data structure for an interpolation point */
typedef struct interp_point{

/* pointers to the previous interpolation point 
   for the first point in the list, prev == NULL */
  struct interp_point *prev;

/* a flag to indicate if the interpolation point is still in use */
  int active;

/* the indices of the interpolation point in this grid */
  int i_point, j_point;
/* the pointer to the other grid */
  struct component_grid *grid_loc;
/* r and s coordinates for the interpolation point in the other grids
   coordinate system */
  real r_loc, s_loc;
/* the indices of the lower left gridpoint in the other grid. The 
   interpolation stencil will use the points 
   i_loc <= i <= i_loc + interp_width - 1, 
   j_loc <= j <= j_loc + interp_width - 1 */
  int i_loc, j_loc;
} interp_point;


/* data structure for a bad grid point */
typedef struct bad_point{
  struct bad_point *prev;
  int i, j;
  int type; 
/* the value of type can be one of the following:
   BOX: Bad discretization point, 
   CROSS: Dead interpolation location,
   FILLED_CIRCLE: Non-explicit interpolation location,
   TRIANGLE: Bad interpolation point because of dead interpolation location,
   INV_TRIANGLE: Bad interpolation point because of non-explicit interpolation location.
*/
} bad_point;


/* data structure for describing a boundary of a component grid */
typedef struct boundary_info{
  struct boundary_info *part1, *part2;
  int start, stop, i_const, j_const;
  real xmin, xmax, ymin, ymax;
} boundary_info;


/* data structure for global description of the boundary */
typedef struct boundary_part{
  int curve_label, n_points, top, inside;
  real_array_2d *coordinates;
} boundary_part;

/* data structure for an overlapping grid */
typedef struct overlapping_grid{

/* number of grids */
  int n_grids;

  struct component_grid *first,*last; /* pointers to the first and the last 
					  component grids */

  int extra;          /* number of fictitious grid points outside each boundary.*/
  int extra_period;   /* number of extra grid points at a periodic boundary. */
/* 
   Observe that extra_period grid points are added at the lower bound but
   extra_period-1 grid points are added at the upper bound. See the above 
   descriptions of range and dim for details.
*/
  int disc_width;     /* width of interior discretisation stencil. Must be odd. */
  int normal_width;   /* width of boundary discretization stencil in the 
			 direction normal to the boundary */
  int tangent_width;  /* width of boundary discretization stencil in the 
			 direction tangential to the boundary */
  int corner_width;   /* width of the stencil close to a corner */
  int interp_width;   /* width of the interpolation stencil. Must be odd. */
  char interp_type;   /* type of interpolation, 
			 ='e': explicit (uncoupled) interpolation 
			 ='i': implicit (coupled) interpolation   */

/* trimming style: minimize overlap = 0, maximize overlap = 1 */
  int trim_style;

  int correct_mismatch;

/* linked list for the global boundary info */
  linked_list *boundary_list;

  real x_min, x_max, y_min, y_max; /* min and max x and y coordinates */

  real max_n_dist;

  char *title;

} overlapping_grid;
  

/* data structure for a component grid */
typedef struct component_grid{

  struct component_grid *next,*prev; /* pointers to next and previous grids in 
				   the composite grid list */

  int priority;             /* priority of this grid in the overlapping grid */

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

/* number of non-fictitious points in the r and s directions, respectively */
  int r_points,s_points;    
  
/* dimensions of arrays in the r and s directions, respectively */
  int r_dim,s_dim;
  
/*  Hence, if r_period == 0: r_dim = r_points + 2*extra 
    and if    r_period == 1: r_dim = r_points + 2*extra_period - 1
    See the description of range below for more details.
*/
  int_array_2d *range_ptr;  /* range of non-fictitious points in the arrays
			       range(is,id) refers to side is=1,2 and direction
			       id=1,2. */
/* If r_period == 0, s_period == 0:
   range(1,1)=extra+1, 
   range(2,1)=r_points+extra, 
   range(1,2)=extra+1,
   range(2,2)=s_points+extra. 

   If r_period == 1:
   range(1,1) = extra_period + 1,
   range(2,1) = r_points + extra_period,

   If s_period == 1:
   range(1,2) = extra_period + 1,
   range(2,2) = s_points + extra_period,
*/

  real r_step,s_step;       /* grid sizes in the r and s-directions, 
			       respectively. r_step = 1/(r_points-1) */

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

  int_array_1d *left_interp_ptr, *right_interp_ptr
    , *lower_interp_ptr, *upper_interp_ptr; 


  int_array_2d *flag_ptr;   /* flag(i,j) is the flag of gridpoint i,j. 
			       flag(i,j)=grid_number: interior point,
			       flag(i,j)=-k: interpolate from grid k
			       flag(i,j)=0: unused point. */
  real_array_2d *x_ptr,*y_ptr;
                            /* x(i,j) and y(i,j) are the x and y-coordinates, 
			       respectively of grid point i,j. */
  int n_interp;             /* number of interpolation points in this 
			       component grid */
  interp_point *last_interp; /* pointer to the last interpolation point. 
				We save the interpolation points on a stack. */

  bad_point *last_bad_point; /* pointer to the stack of bad grid points */

  real x_min, x_max, y_min, y_max; /* min and max x and y coordinates */

  int grid_type;            
/* special type of grid:
   1: Cartesian.
   2: Annular. (NOT USED ANYMORE)
   3. Half-annular. (NOT USED ANYMORE)
   4. Normal-curve mapping. Normals out from any known kind of curve.
   5. Linear interpolation mapping. Linear interpolation between two curves of any
   known kind.
   6. Bi-linear or bi-qubic interpolation between discrete grid points.
   7. Theodorssen-Garrick mapping.
   8. Hyperbolic mapping. If the initial curve is periodic, the mapping 
      will get an O- or a C-grid topology. See the topology info below.
*/

  color_type grid_color;

/* rotation, translation, scaling */
  real cos_theta, sin_theta, x_trans, y_trans, scaling;

/* boundary subdivision information */
  boundary_info *left, *right, *lower, *upper;

/* external segments (in parameter space) on physical boundaries */
  real gap_low_r[2], gap_high_r[2], gap_low_s[2], gap_high_s[2];

/* topology = 0 when the grid mapping is one-to-one or periodic */
/* topology = 1 for a c-grid */
/* branch_direction=1 when the branch cut is aligned with the r-direction and */
/* branch_direction=2 when the branch cut is aligned with the s-direction */
  int topology, branch_direction; 
/* for a c-grid topology, we save the index of the first and last grid points */
/* on the real boundary */
  int first_point, last_point;

} component_grid;

/* define some macros for computing indices of multi-dimensional arrays */
#define bc(i,j)    compute_index_2d(grid_ptr->bc_ptr,i,j)
#define range(i,j) compute_index_2d(grid_ptr->range_ptr,i,j)
#define curve(is,id) compute_index_2d(grid_ptr->curve_ptr,is,id)
#define left_interp(i) compute_index_1d(grid_ptr->left_interp_ptr,i)
#define right_interp(i) compute_index_1d(grid_ptr->right_interp_ptr,i)
#define lower_interp(i) compute_index_1d(grid_ptr->lower_interp_ptr,i)
#define upper_interp(i) compute_index_1d(grid_ptr->upper_interp_ptr,i)

#define flag(i,j)  compute_index_2d(grid_ptr->flag_ptr,i,j)
#define other_flag(i,j)  compute_index_2d(other_ptr->flag_ptr,i,j)


#define x(i,j)     compute_index_2d(grid_ptr->x_ptr,i,j)
#define y(i,j)     compute_index_2d(grid_ptr->y_ptr,i,j)

/* prototype */
int 
overlap(input_output *io_ptr, overlapping_grid *comp_ptr, int redo, 
	int show_phys, int show_holes);
overlapping_grid *
new_overlapping_grid( int n_grids );
void 
delete_overlapping_grid( overlapping_grid *over_ptr );
component_grid *
new_component_grid( char *name);
boundary_part *
new_boundary_part(int curve_label);
boundary_part *
delete_boundary_part( boundary_part *this_part );
void 
insert_grid( component_grid *grid_ptr, overlapping_grid *over_ptr );
void 
bi_linear(int i_loc, int j_loc, grid_point *gp_ptr, component_grid *grid_ptr);
void 
bi_cubic(int i_loc, int j_loc, grid_point *gp_ptr, component_grid *grid_ptr);
void 
plot_overlapping( overlapping_grid *comp_ptr, int plot_mode );
void
plot_global_boundary( boundary_part *boundary_ptr, int plot_mode );
boundary_part *
find_curve( linked_list *boundary_list, int curve_label );
int 
forward_grid_mapping( grid_point *gp_ptr, component_grid *map_ptr );
int 
backward_grid_mapping( grid_point *gp_ptr, component_grid *map_ptr );
void 
subdivide_boundaries( component_grid *grid_ptr );

#endif
