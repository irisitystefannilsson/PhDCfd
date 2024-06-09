#include "xcog.h"

void insert_curve( generic_curve *curve_ptr, xcog_data *xcog_data_ptr){
/* insert the new curve into the linked list */
  if (xcog_data_ptr->first_c != NULL)
    xcog_data_ptr->first_c->prev = curve_ptr;
  curve_ptr->next = xcog_data_ptr->first_c;
  xcog_data_ptr->first_c = curve_ptr;
  if (xcog_data_ptr->last_c == NULL)
    xcog_data_ptr->last_c = curve_ptr;
/* increase the count of mappings */
  xcog_data_ptr->n_curves++;
}
