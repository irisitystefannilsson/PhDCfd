#include "xcog.h"

void 
read_curves_mappings( int32 root, xcog_data *xcog_data_ptr ){
  generic_curve *curve_ptr, *new_curve_ptr;
  generic_mapping *grid_ptr, *new_grid_ptr;
  char *new_curve_name, *new_grid_name, *interp_string, curve_name[80]
    , mapping_name[80];
  int i, n_curves, n_grids;
  int32 curve_dir, mapping_dir, this_curve_dir, this_map_dir, xcog_dir;

/* get the curve directory */
  if ((curve_dir = locate_dir("curves", root)) != -1){
/* read the number of curves */
    hget_int( &n_curves, "n_curves", curve_dir );
    for( i=1; i<=n_curves; i++ ){

      sprintf(curve_name, "curve-%i", i);
/* locate the directory */
      if ( (this_curve_dir = locate_dir(curve_name, curve_dir)) != -1){
/* read the generic and specific curve data */
	if ( (new_curve_ptr = read_curve( this_curve_dir )) ){

/* is the name already in use? */
	  for (curve_ptr = xcog_data_ptr->first_c; curve_ptr != NULL;
	       curve_ptr = curve_ptr->next){
	    if (strcmp(curve_ptr->curve_name, new_curve_ptr->curve_name) == 0){
/* the name is taken, rename curve_ptr */
	      new_curve_name = unique_curve_name(xcog_data_ptr, 
						 new_curve_ptr->curve_name );
	      free( curve_ptr->curve_name );
	      curve_ptr->curve_name = new_curve_name;
	      printf("Warning: The old curve `%s' was renamed to `%s' because of a name "
		     "conflict\n",
		     new_curve_ptr->curve_name, new_curve_name);
	    }
	  } /* end for all current curves */

/* insert the curve in list */
	  insert_curve( new_curve_ptr, xcog_data_ptr );
	} /* end if read_curve */
/* end access to the directory for the present curve */
	Vdetach( this_curve_dir );
      } /* end if locate_dir */

    } /* end for all curves */

/* end access to the directory for all curves */
    Vdetach( curve_dir );
  }

/* get the mappings directory */
  if ((mapping_dir = locate_dir("mappings", root)) != -1){
/* read the number of mappings */
    hget_int( &n_grids, "n_grids", mapping_dir );
    for( i=1; i<=n_grids; i++ ){ 

      sprintf(mapping_name, "mapping-%i", i);
/* locate the directory */
      if ( (this_map_dir = locate_dir(mapping_name, mapping_dir)) != -1){

/* create a new mapping and initialize it from a file */
	if ( (new_grid_ptr = read_mapping( this_map_dir, xcog_data_ptr->first_c )) ){

/* is the mapping name already in use? */
	  for (grid_ptr = xcog_data_ptr->first_m; grid_ptr != NULL; 
	       grid_ptr = grid_ptr->next_m){ 
	    if (strcmp(grid_ptr->grid_name, new_grid_ptr->grid_name) == 0){ 
/* the name is taken, rename grid_ptr */
	      new_grid_name = unique_mapping_name(xcog_data_ptr, 
						  new_grid_ptr->grid_name);
	      free( grid_ptr->grid_name ); 
	      grid_ptr->grid_name = new_grid_name; 
	      printf("Warning: The old mapping `%s' was renamed to `%s' because of a " 
		     "name conflict\n", 
		     new_grid_ptr->grid_name, new_grid_name); 
	    } 
	  } /* end for all existing mappings */

/* insert new_grid_ptr in the mapping list */
	  insert_mapping( new_grid_ptr, xcog_data_ptr );
	} /* end if read_mapping */

/* end access to the directory for the present curve */
	Vdetach( this_map_dir );
      } /* end if locate_dir */
    } /* end for */

/* end access to the directory for all mappings */
    Vdetach( mapping_dir );
  } /* end if locate_dir */

/* get the xcog_parameters directory */
  if ((xcog_dir = locate_dir("xcog_parameters", root)) != -1){
/* read the overlap parameters */
    hget_int( &(xcog_data_ptr->disc_width), "disc_width", xcog_dir );
    hget_int( &(xcog_data_ptr->normal_width), "normal_width", xcog_dir ); 
    hget_int( &(xcog_data_ptr->tangent_width), "tangent_width", xcog_dir );
    hget_int( &(xcog_data_ptr->corner_width), "corner_width", xcog_dir ); 
    hget_int( &(xcog_data_ptr->interp_width), "interp_width", xcog_dir ); 
    hget_int( &(xcog_data_ptr->extra), "extra", xcog_dir ); 

    interp_string = hget_string( "interp_type", xcog_dir );
    if (!strcmp(interp_string, "explicit"))
      xcog_data_ptr->interp_type = 'e';
    else
      xcog_data_ptr->interp_type = 'i';
    free(interp_string);

    hget_real( &(xcog_data_ptr->max_n_dist), "max_n_dist", xcog_dir ); 
    hget_int(  &(xcog_data_ptr->show_physical_boundary), "show_physical_boundary", 
	     xcog_dir);
    if ( !hget_int(  &(xcog_data_ptr->show_holes), "show_holes", xcog_dir) )
      xcog_data_ptr->show_holes = FALSE;
    if ( !hget_int(  &(xcog_data_ptr->correct_mismatch), "correct_mismatch", xcog_dir) )
      xcog_data_ptr->correct_mismatch = TRUE;
    hget_int(  &(xcog_data_ptr->trim_style), "trim_style", xcog_dir );

/* plotting stuff */
    hget_int( &(xcog_data_ptr->comp_plot_mode), "comp_plot_mode", xcog_dir );
    hget_int( &(xcog_data_ptr->grid_plot_mode), "grid_plot_mode", xcog_dir );
    hget_int( &(xcog_data_ptr->curve_plot_mode),"curve_plot_mode", xcog_dir );

    hget_real( &(xcog_data_ptr->xyab[0][0]), "xyab00", xcog_dir);
    hget_real( &(xcog_data_ptr->xyab[1][0]), "xyab10", xcog_dir);
    hget_real( &(xcog_data_ptr->xyab[0][1]), "xyab01", xcog_dir);
    hget_real( &(xcog_data_ptr->xyab[1][1]), "xyab11", xcog_dir);
  
/* title */
    xcog_data_ptr->title = hget_string( "title", xcog_dir );

/* local interpolation */
    hget_int( &(xcog_data_ptr->local_interp), "local_interp", xcog_dir );

/* terminate access to the xcog_parameter directory */
    Vdetach( xcog_dir );
  } /* end if locate dir */
}


void insert_mapping( generic_mapping *grid_ptr, xcog_data *xcog_data_ptr ){
/* insert the new mapping into the linked list */
  if (xcog_data_ptr->first_m != NULL)
    xcog_data_ptr->first_m->prev_m = grid_ptr;
  grid_ptr->next_m = xcog_data_ptr->first_m;
  xcog_data_ptr->first_m = grid_ptr;
  if (xcog_data_ptr->last_m == NULL)
    xcog_data_ptr->last_m = grid_ptr;
/* increase the count of mappings */
  xcog_data_ptr->n_grids++;
}


