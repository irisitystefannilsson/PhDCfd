typedef struct straight_line_data{
  real x_begin, x_end, y_begin, y_end;
} straight_line_data;

straight_line_data *
init_straight_line(generic_curve *curve_ptr, real x_begin, real y_begin, 
		   real x_end, real y_end );
void 
read_straight_line(int32 dir, generic_curve *curve_ptr);

