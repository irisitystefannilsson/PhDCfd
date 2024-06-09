#include "xcog.h"

void load_curves_mappings( int fd, composite_grid *comp_ptr ){
  generic_curve *curve_ptr;
  generic_mapping *grid_ptr;

/* save the number of curves */
  read( fd, (char *) &(comp_ptr->n_curves), sizeof(int) );
  for (curve_ptr = comp_ptr->first_c; curve_ptr != NULL;
       curve_ptr = curve_ptr->next){
/* save the generic curve data */
    read_curve( fd, curve_ptr );
/* save the curve-type dependent data */
    (*curve_ptr->read_specific_curve)( fd, curve_ptr->curve_data_ptr );
  }

/* save the number of mappings */
  read( fd, (char *) &(comp_ptr->n_grids), sizeof(int) );
  for (grid_ptr = comp_ptr->first_m; grid_ptr != NULL;
       grid_ptr = grid_ptr->next_m){
/* save the generic mapping data */
    read_mapping( fd, grid_ptr );
/* save the mapping-type dependent data */
    (*grid_ptr->read_specific_mapping)( fd, grid_ptr->mapping_data_ptr );
  }

}


void read_curve( int fd, generic_curve *curve_ptr ){
  int name_length;

  name_length = strlen( curve_ptr->curve_name );
  read( fd, (char *) &name_length, sizeof(int) );
  read( fd, curve_ptr->curve_name, name_length * sizeof(char) );
/* add the terminating null character to the string */
  curve_ptr->curve_name[name_length] = '\0';
  read( fd, (char *) &(curve_ptr->curve_type), sizeof(int) );
  read( fd, (char *) &(curve_ptr->periodic), sizeof(int) );
  read( fd, (char *) &(curve_ptr->x_min), sizeof(real) );
  read( fd, (char *) &(curve_ptr->x_max), sizeof(real) );
  read( fd, (char *) &(curve_ptr->y_min), sizeof(real) );
  read( fd, (char *) &(curve_ptr->y_max), sizeof(real) );
  read( fd, (char *) &(curve_ptr->used_by_grid), sizeof(int) );
  read( fd, (char *) &(curve_ptr->plot_mode), sizeof(int) );
  read( fd, (char *) curve_ptr->xyab, 4*sizeof(real) );
  read( fd, (char *) &(curve_ptr->r_step), sizeof(real) );
  read( fd, (char *) &(curve_ptr->n_ticks), sizeof(int) );

/* assign the appropriate curve-specific pointers */
  switch( curve_ptr->curve_type ){

  case 0:
/* straight line */
    curve_ptr->curve_data_ptr = malloc( sizeof(straight_line_data) );
    curve_ptr->curve_mapping = straight_line;
    curve_ptr->cleanup_curve_data = cleanup_straight_line;
    curve_ptr->set_curve = set_straight_line;
    curve_ptr->read_specific_curve = read_straight_line;
    curve_ptr->write_specific_curve = write_straight_line;

    break;

  case 1:
/* cubic spline */
    curve_ptr->curve_data_ptr = malloc( sizeof(cubic_spline_data) );
    curve_ptr->curve_mapping = cubic_spline_curve;
    curve_ptr->cleanup_curve_data = cleanup_cubic_spline;
    curve_ptr->set_curve = set_cubic_spline;
    curve_ptr->read_specific_curve = read_cubic_spline;
    curve_ptr->write_specific_curve = write_cubic_spline;
    break;

  case 2:
/* smoothed polygon */
    curve_ptr->curve_data_ptr = malloc( sizeof(smooth_poly_data ) );
    curve_ptr->curve_mapping = smooth_poly_curve;
    curve_ptr->cleanup_curve_data = cleanup_smooth_poly_curve;
    curve_ptr->set_curve = set_smooth_polygon;
    curve_ptr->read_specific_curve = read_smooth_poly;
    curve_ptr->write_specific_curve = write_smooth_poly;
    break;

  case 3:
/* circular arc */
    curve_ptr->curve_data_ptr = malloc( sizeof(circular_arc_data ) );
    curve_ptr->curve_mapping = circular_arc;
    curve_ptr->cleanup_curve_data = cleanup_circular_arc;
    curve_ptr->set_curve = set_circular_arc;
    curve_ptr->read_specific_curve = read_circular_arc;
    curve_ptr->write_specific_curve = write_circular_arc;
    break;

  default:
    ;
  }
}

void read_mapping( int fd, generic_mapping *grid_ptr ){
  int name_length;

  read( fd, (char *) &(grid_ptr->priority), sizeof(int) );
  read( fd, (char *) &name_length, sizeof(int) );
  read( fd, grid_ptr->grid_name, name_length * sizeof(char) );
/* add the terminating null character to the string */
  grid_ptr->grid_name[name_length] = '\0';
  read( fd, (char *) &(grid_ptr->inverse_known), sizeof(int) );
  read( fd, (char *) &(grid_ptr->grid_type), sizeof(int) );
  read( fd, (char *) &(grid_ptr->r_points), sizeof(int) );
  read( fd, (char *) &(grid_ptr->s_points), sizeof(int) );
  read( fd, (char *) &(grid_ptr->r_dim), sizeof(int) );
  read( fd, (char *) &(grid_ptr->s_dim), sizeof(int) );
  read_int_array_2d( grid_ptr->range_ptr, fd );
  read( fd, (char *) &(grid_ptr->r_step), sizeof(real) );
  read( fd, (char *) &(grid_ptr->s_step), sizeof(real) );
  read( fd, (char *) &(grid_ptr->r_period), sizeof(int) );
  read( fd, (char *) &(grid_ptr->s_period), sizeof(int) );
  read_int_array_2d( grid_ptr->bc_ptr, fd );
  read_int_array_2d( grid_ptr->curve_ptr, fd );
  read( fd, (char *) &(grid_ptr->x_min), sizeof(real) );
  read( fd, (char *) &(grid_ptr->x_max), sizeof(real) );
  read( fd, (char *) &(grid_ptr->y_min), sizeof(real) );
  read( fd, (char *) &(grid_ptr->y_max), sizeof(real) );
  read( fd, (char *) &(grid_ptr->grid_plot_mode), sizeof(int) );
  read( fd, (char *) grid_ptr->xyab, 4*sizeof(real) );

}
