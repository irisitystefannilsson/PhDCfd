#include "xcog.h"

void print_mapping_list( xcog_data *xcog_data_ptr ){
  generic_mapping *grid_ptr;

  printf("\n");
  if (xcog_data_ptr->first_m == NULL)
    printf("The mapping list is empty.\n");
  else
    printf("List of all mappings presently known to xcog.\n");
  for (grid_ptr = xcog_data_ptr->first_m; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next_m)
    show_mapping_parameters( grid_ptr );

}


void print_curve_list( xcog_data *xcog_data_ptr ){
  generic_curve *curve_ptr;

  if (xcog_data_ptr->first_c == NULL)
    printf("The curve list is empty.\n");
  else
    printf("List of all curves presently known in the geometry.\n");
  for (curve_ptr = xcog_data_ptr->first_c; curve_ptr != NULL; 
       curve_ptr = curve_ptr->next)
    show_curve_parameters( curve_ptr );

}

