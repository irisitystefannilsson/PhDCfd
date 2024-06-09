/* data structures for xcog */
#ifndef xcog_h
#define xcog_h

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
#include "mappings.h"
#include "overlap.h"

/* data structure for a composite grid */
typedef struct xcog_data{

  int extra;          /* number of fictitious grid points outside each boundary.*/

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

  real max_n_dist;

/* min and max x and y coordinates */
  real x_min, x_max, y_min, y_max; 

/* number of mappings */
  int n_grids;

/* pointers to the first and the last mappings */
  struct generic_mapping *first_m,*last_m; 

/* number of curves defined */
  int n_curves;

/* pointer to the first and last curve in the geometry description */
  struct generic_curve *first_c, *last_c;

/* plotting information */
  int comp_plot_mode, grid_plot_mode, curve_plot_mode;
  real xyab[2][2];

  char *title;

  int local_interp;

  int show_physical_boundary, show_holes, correct_mismatch;
} xcog_data;
  

void print_mapping_list(xcog_data *comp_ptr);
void print_curve_list(xcog_data *comp_ptr);
void overlap_parameters(input_output *io_ptr, xcog_data *comp_ptr);
void set_priority( input_output *io_ptr, xcog_data *comp_ptr );
void show_composite_parameters(xcog_data *comp_ptr);
xcog_data *new_xcog_data( void );
void set_global_boundaries( xcog_data *comp_ptr );
void show_all( xcog_data *comp_ptr );

void insert_curve( generic_curve *curve_ptr, xcog_data *comp_ptr);
void 
save_curves_mappings( int32 root, xcog_data *comp_ptr );
void 
read_curves_mappings( int32 root, xcog_data *comp_ptr );

char *unique_curve_name( xcog_data *comp_ptr, char *name );
char *unique_mapping_name( xcog_data *comp_ptr, char *name );
void insert_mapping( generic_mapping *grid_ptr, xcog_data *comp_ptr );

void show_overlap_param( xcog_data *comp_ptr );
void welcome_to_xcog( void );
void general_help( void );
void tutorial_intro( void );
void run_tutorial( input_output *io_ptr );
void print_help( char *brief_help );

void save_xcog_data( input_output *io_ptr, overlapping_grid *comp_ptr );
overlapping_grid *compute_overlap(input_output *io_ptr, xcog_data *comp_ptr, 
				  overlapping_grid *over_ptr, int *overlap_ok);
void set_global_plot_mode(input_output *io_ptr, xcog_data *comp_ptr,
			  overlapping_grid *over_ptr);
void component_grid_list( input_output *io_ptr, xcog_data *xcog_data_ptr );
void 
save_grid_plot3d(generic_mapping *map_ptr, FILE *fp, int i_min, int i_max, 
		 int j_min, int j_max );

#endif
