#include "xcog.h"

char *unique_curve_name( xcog_data *xcog_data_ptr, char *name ){
  int i, name_length;
  char *new_name;
  generic_curve *curve_ptr;

  name_length = strlen(name) + 7;
  new_name = (char *) malloc( (name_length+1) * sizeof(char) );

  i = 1;
  do{
    i++;
    sprintf(new_name, "%s_%i", name, i);
    for (curve_ptr = xcog_data_ptr->first_c; curve_ptr != NULL &&
         strcmp(curve_ptr->curve_name, new_name) != 0; curve_ptr = curve_ptr->next);
  } while( curve_ptr != NULL );

  return new_name;
}

char *unique_mapping_name( xcog_data *xcog_data_ptr, char *name ){
  int i, name_length;
  char *new_name;
  generic_mapping *grid_ptr;

  name_length = strlen(name) + 7;
  new_name = (char *) malloc( (name_length+1) * sizeof(char) );

  i = 1;
  do{
    i++;
    sprintf(new_name, "%s_%i", name, i);
    for (grid_ptr = xcog_data_ptr->first_m; grid_ptr != NULL &&
         strcmp(grid_ptr->grid_name, new_name) != 0; grid_ptr = grid_ptr->next_m);
  } while( grid_ptr != NULL );

  return new_name;
}

