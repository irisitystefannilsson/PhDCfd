#include "overlap.h"

#define R0 ((i_loc-range(1,1))   * grid_ptr->r_step)
#define R1 ((i_loc+1-range(1,1)) * grid_ptr->r_step)
#define ALPHA_0(r)   ((R1 - r)/grid_ptr->r_step)
#define D_ALPHA_0(r) (-1.0/grid_ptr->r_step)
#define ALPHA_1(r)   ((r - R0)/grid_ptr->r_step)
#define D_ALPHA_1(r) (1.0/grid_ptr->r_step)

#define S0 ((j_loc-range(1,2))  * grid_ptr->s_step)
#define S1 ((j_loc+1-range(1,2))* grid_ptr->s_step)
#define BETA_0(s)   ((S1 - s)/grid_ptr->s_step)
#define D_BETA_0(s) (-1.0/grid_ptr->s_step)
#define BETA_1(s)   ((s - S0)/grid_ptr->s_step)
#define D_BETA_1(s) (1.0/grid_ptr->s_step)

#define LIN_X(r,s) (ALPHA_0(r) * BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		    ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		    ALPHA_0(r) * BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		    ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) )
#define LIN_X_R(r,s) (D_ALPHA_0(r) * BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		      D_ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		      D_ALPHA_0(r) * BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		      D_ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) )
#define LIN_X_S(r,s) (ALPHA_0(r) * D_BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		      ALPHA_1(r) * D_BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		      ALPHA_0(r) * D_BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		      ALPHA_1(r) * D_BETA_1(s) * x(i_loc+1,j_loc+1) )

#define LIN_Y(r,s) (ALPHA_0(r) * BETA_0(s) * y(i_loc  ,j_loc  ) +  \
		    ALPHA_1(r) * BETA_0(s) * y(i_loc+1,j_loc  ) +  \
		    ALPHA_0(r) * BETA_1(s) * y(i_loc  ,j_loc+1) +  \
		    ALPHA_1(r) * BETA_1(s) * y(i_loc+1,j_loc+1) )
#define LIN_Y_R(r,s) (D_ALPHA_0(r) * BETA_0(s) * y(i_loc  ,j_loc  ) +  \
		      D_ALPHA_1(r) * BETA_0(s) * y(i_loc+1,j_loc  ) +  \
		      D_ALPHA_0(r) * BETA_1(s) * y(i_loc  ,j_loc+1) +  \
		      D_ALPHA_1(r) * BETA_1(s) * y(i_loc+1,j_loc+1) )
#define LIN_Y_S(r,s) (ALPHA_0(r) * D_BETA_0(s) * y(i_loc  ,j_loc  ) +  \
		      ALPHA_1(r) * D_BETA_0(s) * y(i_loc+1,j_loc  ) +  \
		      ALPHA_0(r) * D_BETA_1(s) * y(i_loc  ,j_loc+1) +  \
		      ALPHA_1(r) * D_BETA_1(s) * y(i_loc+1,j_loc+1) )

void bi_linear(int i_loc, int j_loc, grid_point *gp_ptr, 
		      component_grid *grid_ptr){
  gp_ptr->x  = LIN_X(gp_ptr->r, gp_ptr->s);
  gp_ptr->y  = LIN_Y(gp_ptr->r, gp_ptr->s);
  gp_ptr->xr = LIN_X_R(gp_ptr->r, gp_ptr->s);
  gp_ptr->xs = LIN_X_S(gp_ptr->r, gp_ptr->s);
  gp_ptr->yr = LIN_Y_R(gp_ptr->r, gp_ptr->s);
  gp_ptr->ys = LIN_Y_S(gp_ptr->r, gp_ptr->s);
}
