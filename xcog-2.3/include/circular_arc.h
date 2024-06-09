typedef struct circular_arc_data{
  real pi, theta_0, theta_1, xcenter, ycenter, radius;
} circular_arc_data;

circular_arc_data *
init_circular_arc(generic_curve *curve_ptr, real theta_0, real theta_1, real radius,
		  real xcenter, real ycenter );
void 
read_circular_arc( int32 dir, generic_curve *curve_ptr );
