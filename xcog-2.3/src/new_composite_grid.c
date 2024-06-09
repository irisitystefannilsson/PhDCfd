#include "xcog.h"

xcog_data *
new_xcog_data( void ){
  xcog_data *xcog_data_ptr;

  if ((xcog_data_ptr = (xcog_data *)malloc(sizeof(xcog_data))) == NULL)
    printf("memory error in new_xcog_data\n");
  
/* no grids initially */
  xcog_data_ptr->n_grids = 0;
/*  xcog_data_ptr->first   = NULL;
  xcog_data_ptr->last    = NULL;*/

/* number of ghost points */
  xcog_data_ptr->extra = 0;

/* discretization width 3 */
  xcog_data_ptr->disc_width = 3;    /* must be odd */
  xcog_data_ptr->normal_width = 3;
  xcog_data_ptr->tangent_width = 3; /* must be odd */
  xcog_data_ptr->corner_width = 3;
    
/* interpolation width 2: bi-linear interpolation */
  xcog_data_ptr->interp_width = 2;
    
/* set the interpolation type */
  xcog_data_ptr->interp_type = 'e';
    
/* default trimming style is to minimize the overlap */
  xcog_data_ptr->trim_style = 0;

/* default additional tolerance for two curves to be regarded as the same curve */
  xcog_data_ptr->max_n_dist = 0.0;

/* the default is NOT to look at the physical boundary and NOT to look */
/* at the used grid points after the holes have been cut */
  xcog_data_ptr->show_physical_boundary = FALSE;
  xcog_data_ptr->show_holes = FALSE;
  xcog_data_ptr->correct_mismatch = TRUE;

/* initialize bounding box */
  xcog_data_ptr->x_min = 0.0;
  xcog_data_ptr->x_max = 1.0;
  xcog_data_ptr->y_min = 0.0;
  xcog_data_ptr->y_max = 1.0;

/* initialize the pointers to the first and last mappings */
  xcog_data_ptr->first_m = NULL;
  xcog_data_ptr->last_m  = NULL;

/* initialize the pointers to the first and last curves */
  xcog_data_ptr->first_c = NULL;
  xcog_data_ptr->last_c  = NULL;
/* no curves initially */
  xcog_data_ptr->n_curves = 0;

/* plotting stuff */
  xcog_data_ptr->comp_plot_mode = 0;
  xcog_data_ptr->grid_plot_mode = 0;
  xcog_data_ptr->curve_plot_mode = 0;

  xcog_data_ptr->title = NULL;

  xcog_data_ptr->xyab[0][0] = 0.0;
  xcog_data_ptr->xyab[0][1] = 1.0;
  xcog_data_ptr->xyab[1][0] = 0.0;
  xcog_data_ptr->xyab[1][1] = 1.0;

  xcog_data_ptr->local_interp = 0;

  return xcog_data_ptr;
}




