#ifndef curves_h
#define curves_h

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
#include "hdf_stuff.h"
#include "input_output.h"

/* data structure for reporting a point on a curve */
typedef struct curve_point{
  real r, x, xr, xrr, y, yr, yrr;
} curve_point;


typedef struct generic_curve{

/* pointers to previous and next curves in the list */
  struct generic_curve *prev, *next;

/* name */
  char *curve_name;

/* type */
  int curve_type; /*  curve_type     note
		      0              straight line
		      1              cubic spline
		      2              smoothed polygon 
                      3              circular arc */

/* periodicity flag */
  int periodic;

/* min and max x and y coordinates before rotation, translation and scaling*/
  real x_min_0, x_max_0, y_min_0, y_max_0; 

/* min and max x and y coordinates after rotation, translation and scaling */
  real x_min, x_max, y_min, y_max; 

/* pointer to the specific curve data */
  void *curve_data_ptr;

/* pointer to the specific curve evaluation function */
  int (*curve_mapping)(curve_point *, void *);

/* pointer to the specific curve cleanup function */
  void *(*cleanup_curve_data)(void *);

/* pointer to the function for interactive modification of an existing curve
   of the present type */
  void (*set_curve)( input_output *, struct generic_curve * );

/* pointers to the functions for saving and recovering a curve from a binary file */
  void (*write_specific_curve)( int32 dir, void *curve_data_ptr );

/* pointer to the function that plots the control points */
  void (*plot_control_points)( struct generic_curve *curve_ptr );

/* pointer to function that copies and normalizes the parameter values at all */
/* node points of the curve. The pointer to the vector is returned by the function */
/* and the size of the vector is returned in *n. */
  real *(*copy_nodes)( struct generic_curve *curve_ptr, int *n );

/* pointer to function for copying the specific data */
  void *(*copy_curve_data)( void *curve_data_ptr );

/* is this curve used in the definition of a grid ? */
  int used_by_grid;

/* plot information */
  int plot_mode;
  real xyab[2][2];
  color_type curve_color;

/* step size for plotting it */
  real r_step;

/* number of tickmarks to be plotted */
  int n_ticks;

/* rotation, translation, scaling */
  real theta_rot, cos_theta, sin_theta, x_trans, y_trans, scaling;

} generic_curve;

/* curves known to xcog */
#include "straight_line.h"
#include "smooth_poly_curve.h"
#include "cubic_spline.h"
#include "circular_arc.h"

generic_curve *
read_curve( int32 curve_dir );
void 
write_curve( int32 curve_dir, generic_curve *curve_ptr );
generic_curve *choose_curve( input_output *io_ptr, char *name );
void delete_generic_curve( generic_curve *curve_ptr);
generic_curve *new_generic_curve(char *name);
void plot_curve( generic_curve *curve_ptr, int plot_mode );
void set_curve_plot_mode( input_output *io_ptr, generic_curve *curve_ptr );
void show_curve_parameters( generic_curve *curve_ptr );
generic_curve *get_curve_ptr( input_output *io_ptr, generic_curve *first_c);
int get_2curve_ptr( input_output *io_ptr, generic_curve *first_c, 
		   generic_curve **curve_1_ptr, generic_curve **curve_2_ptr, 
		   int level);
int curve_function( curve_point *cp_ptr, generic_curve *curve_ptr );
void curve_view_port( generic_curve *curve_ptr );
void set_curve_trans(input_output *io_ptr, generic_curve *curve_ptr);
void rot_trans_scale( real x_old, real y_old, real *x_new, real *y_new, 
		     generic_curve *curve_ptr );
void inv_rot_trans_scale( real x_old, real y_old, real *x_new, real *y_new, 
			 generic_curve *curve_ptr );
generic_curve *copy_curve( char *name, generic_curve *old_curve_ptr);
void set_curve( input_output *io_ptr, generic_curve *curve_ptr );
real *copy_nodes( generic_curve *curve_ptr, int *n_nodes );
void plot_control_points( generic_curve *curve_ptr );

#endif
