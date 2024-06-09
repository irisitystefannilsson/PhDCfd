#include "xcog.h"

composite_grid *init_composite_grid( void ){
  composite_grid *comp_ptr;

  comp_ptr = (composite_grid *)malloc(sizeof(composite_grid));
  
/* no grids initially */
  comp_ptr->n_grids = 0;
  comp_ptr->first   = NULL;
  comp_ptr->last    = NULL;

/* number of ghost points */
  comp_ptr->extra = 0;

/* discretization width 3 */
  comp_ptr->disc_width = 3;    /* must be odd */
  comp_ptr->normal_width = 3;
  comp_ptr->tangent_width = 3; /* must be odd */
  comp_ptr->corner_width = 3;
    
/* interpolation width 3 */
  comp_ptr->interp_width = 3;  /* must be odd */
    
/* set the interpolation type */
  comp_ptr->interp_type = 'e';
    
/* number of extra points at periodic boundaries */
  comp_ptr->extra_period = int_max((comp_ptr->disc_width - 1)/2, 
				   (comp_ptr->interp_width - 1)/2);

/* initialize bounding box */
  comp_ptr->x_min =  1.0e10;
  comp_ptr->x_max = -1.0e10;
  comp_ptr->y_min =  1.0e10;
  comp_ptr->y_max = -1.0e10;

/* initialize the pointers to the first and last curves */
  comp_ptr->first_c = NULL;
  comp_ptr->last_c  = NULL;
/* no curves initially */
  comp_ptr->n_curves = 0;

  return comp_ptr;
}




