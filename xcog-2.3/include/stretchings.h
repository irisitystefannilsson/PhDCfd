#ifndef stretchings_h
#define stretchings_h

#include "curves.h"

typedef struct generic_stretching{

  int stretch_type;
/* stretch_type: */
/* 1 : arclength stretching */
/* 2 : exponential stretching */
/* 3 : layer stretching */
/* 4 : hyperbolic tangent stretching */

  void *stretch_data_ptr;
  void (*uniform_to_stretch)( real uniform, struct generic_stretching *stretch_ptr, 
			     real *stretch, real *d_stretch, real *d2_stretch );
  void (*write_stretch_data)( int32 dir, struct generic_stretching *stretch_ptr );
  void (*set_stretching)( input_output *io_ptr, 
			 struct generic_stretching *stretch_ptr,
			 generic_curve *curve_ptr,
			 int *grid_points);
  void (*cleanup_stretching)( struct generic_stretching *stretch_ptr );
  void *(*copy_stretching)( void *stretch_data_ptr );
} generic_stretching;


generic_stretching *choose_stretching( input_output *io_ptr, 
				      generic_curve *curve_ptr, int *grid_points );
char *stretching_name( generic_stretching *stretch_ptr );
void plot_stretch( generic_stretching *stretch_ptr, generic_curve *curve_ptr, 
		  int grid_points, int plot_mode );
generic_stretching *
read_stretching( int32 dir );
void 
write_stretch( int32 dir, generic_stretching *stretch_ptr );
generic_stretching *delete_generic_stretching( generic_stretching *stretch_ptr );
generic_stretching *copy_stretching( generic_stretching *old_stretch_ptr );
int set_stretching( input_output *io_ptr, generic_stretching *stretch_ptr,
		    generic_curve *curve_ptr, int *grid_points );
void uniform_to_stretch( real uniform, generic_stretching *stretch_ptr, 
			real *stretch, real *d_stretch, real *d2_stretch );

/* stretching */
#include "layer_stretch.h"
#include "arcl_stretch.h"
#include "exponential_stretch.h"
#include "tanh_stretch.h"

#endif
