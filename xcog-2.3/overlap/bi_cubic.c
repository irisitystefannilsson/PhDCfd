#include "overlap.h"

#define R0 ((i_loc-range(1,1))   * grid_ptr->r_step)
#define R1 ((i_loc+1-range(1,1)) * grid_ptr->r_step)
#define R2 ((i_loc+2-range(1,1)) * grid_ptr->r_step)
#define R3 ((i_loc+3-range(1,1)) * grid_ptr->r_step)
#define DR3 (grid_ptr->r_step*grid_ptr->r_step*grid_ptr->r_step)
#define ALPHA_0(r)   ((R3 - r)*(R2 - r)*(R1 - r)/(6*DR3))
#define ALPHA_1(r)   ((R3 - r)*(R2 - r)*(r - R0)/(2*DR3))
#define ALPHA_2(r)   ((R3 - r)*(r - R1)*(r - R0)/(2*DR3))
#define ALPHA_3(r)   ((r - R2)*(r - R1)*(r - R0)/(6*DR3))
#define D_ALPHA_0(r) (( -(R2 - r)*(R1 - r) - (R3 - r)*(R1 - r) - (R3 - r)*(R2 - r) ) \
		      /(6*DR3))
#define D_ALPHA_1(r) (( -(R2 - r)*(r - R0) - (R3 - r)*(r - R0) + (R3 - r)*(R2 - r) ) \
		      /(2*DR3))
#define D_ALPHA_2(r) (( -(r - R1)*(r - R0) + (R3 - r)*(r - R0) + (R3 - r)*(r - R1) ) \
		      /(2*DR3))
#define D_ALPHA_3(r) (( +(r - R1)*(r - R0) + (r - R2)*(r - R0) + (r - R2)*(r - R1) ) \
		      /(6*DR3))

#define S0 ((j_loc-range(1,2))  * grid_ptr->s_step)
#define S1 ((j_loc+1-range(1,2))* grid_ptr->s_step)
#define S2 ((j_loc+2-range(1,2)) * grid_ptr->s_step)
#define S3 ((j_loc+3-range(1,2)) * grid_ptr->s_step)
#define DS3 (grid_ptr->s_step*grid_ptr->s_step*grid_ptr->s_step)
#define BETA_0(s)   ((S3 - s)*(S2 - s)*(S1 - s)/(6*DS3))
#define BETA_1(s)   ((S3 - s)*(S2 - s)*(s - S0)/(2*DS3))
#define BETA_2(s)   ((S3 - s)*(s - S1)*(s - S0)/(2*DS3))
#define BETA_3(s)   ((s - S2)*(s - S1)*(s - S0)/(6*DS3))
#define D_BETA_0(s) (( -(S2 - s)*(S1 - s) - (S3 - s)*(S1 - s) - (S3 - s)*(S2 - s) ) \
		      /(6*DS3))
#define D_BETA_1(s) (( -(S2 - s)*(s - S0) - (S3 - s)*(s - S0) + (S3 - s)*(S2 - s) ) \
		      /(2*DS3))
#define D_BETA_2(s) (( -(s - S1)*(s - S0) + (S3 - s)*(s - S0) + (S3 - s)*(s - S1) ) \
		      /(2*DS3))
#define D_BETA_3(s) (( +(s - S1)*(s - S0) + (s - S2)*(s - S0) + (s - S2)*(s - S1) ) \
		      /(6*DS3))

#define CUB_X(r,s) (ALPHA_0(r) * BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		    ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		    ALPHA_2(r) * BETA_0(s) * x(i_loc+2,j_loc  ) +  \
		    ALPHA_3(r) * BETA_0(s) * x(i_loc+3,j_loc  ) +  \
		    ALPHA_0(r) * BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		    ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) +  \
		    ALPHA_2(r) * BETA_1(s) * x(i_loc+2,j_loc+1) +  \
		    ALPHA_3(r) * BETA_1(s) * x(i_loc+3,j_loc+1) +  \
		    ALPHA_0(r) * BETA_2(s) * x(i_loc  ,j_loc+2) +  \
		    ALPHA_1(r) * BETA_2(s) * x(i_loc+1,j_loc+2) +  \
		    ALPHA_2(r) * BETA_2(s) * x(i_loc+2,j_loc+2) +  \
		    ALPHA_3(r) * BETA_2(s) * x(i_loc+3,j_loc+2) +  \
		    ALPHA_0(r) * BETA_3(s) * x(i_loc  ,j_loc+3) +  \
		    ALPHA_1(r) * BETA_3(s) * x(i_loc+1,j_loc+3) +  \
		    ALPHA_2(r) * BETA_3(s) * x(i_loc+2,j_loc+3) +  \
		    ALPHA_3(r) * BETA_3(s) * x(i_loc+3,j_loc+3) )
#define CUB_X_R(r,s) (D_ALPHA_0(r) * BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		      D_ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		      D_ALPHA_2(r) * BETA_0(s) * x(i_loc+2,j_loc  ) +  \
		      D_ALPHA_3(r) * BETA_0(s) * x(i_loc+3,j_loc  ) +  \
		      D_ALPHA_0(r) * BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		      D_ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) +  \
		      D_ALPHA_2(r) * BETA_1(s) * x(i_loc+2,j_loc+1) +  \
		      D_ALPHA_3(r) * BETA_1(s) * x(i_loc+3,j_loc+1) +  \
		      D_ALPHA_0(r) * BETA_2(s) * x(i_loc  ,j_loc+2) +  \
		      D_ALPHA_1(r) * BETA_2(s) * x(i_loc+1,j_loc+2) +  \
		      D_ALPHA_2(r) * BETA_2(s) * x(i_loc+2,j_loc+2) +  \
		      D_ALPHA_3(r) * BETA_2(s) * x(i_loc+3,j_loc+2) +  \
		      D_ALPHA_0(r) * BETA_3(s) * x(i_loc  ,j_loc+3) +  \
		      D_ALPHA_1(r) * BETA_3(s) * x(i_loc+1,j_loc+3) +  \
		      D_ALPHA_2(r) * BETA_3(s) * x(i_loc+2,j_loc+3) +  \
		      D_ALPHA_3(r) * BETA_3(s) * x(i_loc+3,j_loc+3) )
#define CUB_X_S(r,s) (ALPHA_0(r) * D_BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		      ALPHA_1(r) * D_BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		      ALPHA_2(r) * D_BETA_0(s) * x(i_loc+2,j_loc  ) +  \
		      ALPHA_3(r) * D_BETA_0(s) * x(i_loc+3,j_loc  ) +  \
		      ALPHA_0(r) * D_BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		      ALPHA_1(r) * D_BETA_1(s) * x(i_loc+1,j_loc+1) +  \
		      ALPHA_2(r) * D_BETA_1(s) * x(i_loc+2,j_loc+1) +  \
		      ALPHA_3(r) * D_BETA_1(s) * x(i_loc+3,j_loc+1) +  \
		      ALPHA_0(r) * D_BETA_2(s) * x(i_loc  ,j_loc+2) +  \
		      ALPHA_1(r) * D_BETA_2(s) * x(i_loc+1,j_loc+2) +  \
		      ALPHA_2(r) * D_BETA_2(s) * x(i_loc+2,j_loc+2) +  \
		      ALPHA_3(r) * D_BETA_2(s) * x(i_loc+3,j_loc+2) +  \
		      ALPHA_0(r) * D_BETA_3(s) * x(i_loc  ,j_loc+3) +  \
		      ALPHA_1(r) * D_BETA_3(s) * x(i_loc+1,j_loc+3) +  \
		      ALPHA_2(r) * D_BETA_3(s) * x(i_loc+2,j_loc+3) +  \
		      ALPHA_3(r) * D_BETA_3(s) * x(i_loc+3,j_loc+3) )


#define CUB_Y(r,s) (ALPHA_0(r) * BETA_0(s) * y(i_loc  ,j_loc  ) +  \
		    ALPHA_1(r) * BETA_0(s) * y(i_loc+1,j_loc  ) +  \
		    ALPHA_2(r) * BETA_0(s) * y(i_loc+2,j_loc  ) +  \
		    ALPHA_3(r) * BETA_0(s) * y(i_loc+3,j_loc  ) +  \
		    ALPHA_0(r) * BETA_1(s) * y(i_loc  ,j_loc+1) +  \
		    ALPHA_1(r) * BETA_1(s) * y(i_loc+1,j_loc+1) +  \
		    ALPHA_2(r) * BETA_1(s) * y(i_loc+2,j_loc+1) +  \
		    ALPHA_3(r) * BETA_1(s) * y(i_loc+3,j_loc+1) +  \
		    ALPHA_0(r) * BETA_2(s) * y(i_loc  ,j_loc+2) +  \
		    ALPHA_1(r) * BETA_2(s) * y(i_loc+1,j_loc+2) +  \
		    ALPHA_2(r) * BETA_2(s) * y(i_loc+2,j_loc+2) +  \
		    ALPHA_3(r) * BETA_2(s) * y(i_loc+3,j_loc+2) +  \
		    ALPHA_0(r) * BETA_3(s) * y(i_loc  ,j_loc+3) +  \
		    ALPHA_1(r) * BETA_3(s) * y(i_loc+1,j_loc+3) +  \
		    ALPHA_2(r) * BETA_3(s) * y(i_loc+2,j_loc+3) +  \
		    ALPHA_3(r) * BETA_3(s) * y(i_loc+3,j_loc+3) )
#define CUB_Y_R(r,s) (D_ALPHA_0(r) * BETA_0(s) * y(i_loc  ,j_loc  ) +  \
		      D_ALPHA_1(r) * BETA_0(s) * y(i_loc+1,j_loc  ) +  \
		      D_ALPHA_2(r) * BETA_0(s) * y(i_loc+2,j_loc  ) +  \
		      D_ALPHA_3(r) * BETA_0(s) * y(i_loc+3,j_loc  ) +  \
		      D_ALPHA_0(r) * BETA_1(s) * y(i_loc  ,j_loc+1) +  \
		      D_ALPHA_1(r) * BETA_1(s) * y(i_loc+1,j_loc+1) +  \
		      D_ALPHA_2(r) * BETA_1(s) * y(i_loc+2,j_loc+1) +  \
		      D_ALPHA_3(r) * BETA_1(s) * y(i_loc+3,j_loc+1) +  \
		      D_ALPHA_0(r) * BETA_2(s) * y(i_loc  ,j_loc+2) +  \
		      D_ALPHA_1(r) * BETA_2(s) * y(i_loc+1,j_loc+2) +  \
		      D_ALPHA_2(r) * BETA_2(s) * y(i_loc+2,j_loc+2) +  \
		      D_ALPHA_3(r) * BETA_2(s) * y(i_loc+3,j_loc+2) +  \
		      D_ALPHA_0(r) * BETA_3(s) * y(i_loc  ,j_loc+3) +  \
		      D_ALPHA_1(r) * BETA_3(s) * y(i_loc+1,j_loc+3) +  \
		      D_ALPHA_2(r) * BETA_3(s) * y(i_loc+2,j_loc+3) +  \
		      D_ALPHA_3(r) * BETA_3(s) * y(i_loc+3,j_loc+3) )
#define CUB_Y_S(r,s) (ALPHA_0(r) * D_BETA_0(s) * y(i_loc  ,j_loc  ) +  \
		      ALPHA_1(r) * D_BETA_0(s) * y(i_loc+1,j_loc  ) +  \
		      ALPHA_2(r) * D_BETA_0(s) * y(i_loc+2,j_loc  ) +  \
		      ALPHA_3(r) * D_BETA_0(s) * y(i_loc+3,j_loc  ) +  \
		      ALPHA_0(r) * D_BETA_1(s) * y(i_loc  ,j_loc+1) +  \
		      ALPHA_1(r) * D_BETA_1(s) * y(i_loc+1,j_loc+1) +  \
		      ALPHA_2(r) * D_BETA_1(s) * y(i_loc+2,j_loc+1) +  \
		      ALPHA_3(r) * D_BETA_1(s) * y(i_loc+3,j_loc+1) +  \
		      ALPHA_0(r) * D_BETA_2(s) * y(i_loc  ,j_loc+2) +  \
		      ALPHA_1(r) * D_BETA_2(s) * y(i_loc+1,j_loc+2) +  \
		      ALPHA_2(r) * D_BETA_2(s) * y(i_loc+2,j_loc+2) +  \
		      ALPHA_3(r) * D_BETA_2(s) * y(i_loc+3,j_loc+2) +  \
		      ALPHA_0(r) * D_BETA_3(s) * y(i_loc  ,j_loc+3) +  \
		      ALPHA_1(r) * D_BETA_3(s) * y(i_loc+1,j_loc+3) +  \
		      ALPHA_2(r) * D_BETA_3(s) * y(i_loc+2,j_loc+3) +  \
		      ALPHA_3(r) * D_BETA_3(s) * y(i_loc+3,j_loc+3) )

void bi_cubic(int i_loc, int j_loc, grid_point *gp_ptr, 
		      component_grid *grid_ptr){
  gp_ptr->x  = CUB_X(gp_ptr->r, gp_ptr->s);
  gp_ptr->y  = CUB_Y(gp_ptr->r, gp_ptr->s);
  gp_ptr->xr = CUB_X_R(gp_ptr->r, gp_ptr->s);
  gp_ptr->xs = CUB_X_S(gp_ptr->r, gp_ptr->s);
  gp_ptr->yr = CUB_Y_R(gp_ptr->r, gp_ptr->s);
  gp_ptr->ys = CUB_Y_S(gp_ptr->r, gp_ptr->s);
}
