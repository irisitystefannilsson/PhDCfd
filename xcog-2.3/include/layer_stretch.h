typedef struct layer_stretch_info {
  int n_tanh, n_nodes;
  real scale;
  real *strength, *focus, *u0, *str_sharp;
  real *str_node, *uni_node;
} layer_stretch_info;

void 
init_layer_stretch( generic_stretching *stretch_ptr, generic_curve *curve_ptr );
void 
read_layer_stretch( int32 dir, generic_stretching *stretch_ptr );
