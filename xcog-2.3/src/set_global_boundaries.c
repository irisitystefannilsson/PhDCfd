#include "xcog.h"

void set_global_boundaries( xcog_data *xcog_data_ptr ){
  generic_mapping *grid_ptr;
  generic_curve *curve_ptr;

  if (xcog_data_ptr->n_curves < 1 && xcog_data_ptr->n_grids < 1){
    xcog_data_ptr->x_min = 0.0;
    xcog_data_ptr->x_max = 1.0;
    xcog_data_ptr->y_min = 0.0;
    xcog_data_ptr->y_max = 1.0;
  }
  else{

/* initialize bounding box */
    xcog_data_ptr->x_min =  1.0e10;
    xcog_data_ptr->x_max = -1.0e10;
    xcog_data_ptr->y_min =  1.0e10;
    xcog_data_ptr->y_max = -1.0e10;

/* check all grids */
    for (grid_ptr = xcog_data_ptr->first_m; grid_ptr != NULL;
	 grid_ptr = grid_ptr->next_m){
      if (xcog_data_ptr->x_min > grid_ptr->x_min)
	xcog_data_ptr->x_min = grid_ptr->x_min;
      if (xcog_data_ptr->x_max < grid_ptr->x_max)
	xcog_data_ptr->x_max = grid_ptr->x_max;
      if (xcog_data_ptr->y_min > grid_ptr->y_min)
	xcog_data_ptr->y_min = grid_ptr->y_min;
      if (xcog_data_ptr->y_max < grid_ptr->y_max)
	xcog_data_ptr->y_max = grid_ptr->y_max;
    }

/* check all curves */
    for (curve_ptr = xcog_data_ptr->first_c; curve_ptr != NULL;
	 curve_ptr = curve_ptr->next){
      if (xcog_data_ptr->x_min > curve_ptr->x_min)
	xcog_data_ptr->x_min = curve_ptr->x_min;
      if (xcog_data_ptr->x_max < curve_ptr->x_max)
	xcog_data_ptr->x_max = curve_ptr->x_max;
      if (xcog_data_ptr->y_min > curve_ptr->y_min)
	xcog_data_ptr->y_min = curve_ptr->y_min;
      if (xcog_data_ptr->y_max < curve_ptr->y_max)
	xcog_data_ptr->y_max = curve_ptr->y_max;
    }
  }
}

void
show_all( xcog_data *xcog_data_ptr ){
  xcog_data_ptr->xyab[0][0] = xcog_data_ptr->x_min - 
    0.05*(xcog_data_ptr->x_max - xcog_data_ptr->x_min);
  xcog_data_ptr->xyab[0][1] = xcog_data_ptr->x_max + 
    0.05*(xcog_data_ptr->x_max - xcog_data_ptr->x_min);
  xcog_data_ptr->xyab[1][0] = xcog_data_ptr->y_min - 
    0.05*(xcog_data_ptr->y_max - xcog_data_ptr->y_min);
  xcog_data_ptr->xyab[1][1] = xcog_data_ptr->y_max + 
    0.05*(xcog_data_ptr->y_max - xcog_data_ptr->y_min);
}

