#include "xcog.h"
#include "hdf_stuff.h"

void 
save_curves_mappings( int32 root, xcog_data *xcog_data_ptr ){
  generic_curve *curve_ptr;
  generic_mapping *grid_ptr;
  int counter;
  int32 curve_dir, mapping_dir, xcog_dir, one_dir;
  char curve_name[80], mapping_name[80];

/* make a directory for the curves */
  if ((curve_dir = create_dir("curves", "curves", root)) != -1){
/* save the number of curves */
    hput_int(xcog_data_ptr->n_curves, "n_curves", curve_dir);
/* save each specific curve */
    counter = 0;
    for (curve_ptr = xcog_data_ptr->first_c; curve_ptr != NULL;
	 curve_ptr = curve_ptr->next){
      sprintf(curve_name, "curve-%i", ++counter);
      if ((one_dir = create_dir(curve_name, "one curve", curve_dir)) != -1){
/* save the generic and specific curve data */
	write_curve( one_dir, curve_ptr );
/* release the sub-directory */
	Vdetach(one_dir);
      }
      else
	printf("Error: save_curves_mappings: unable to create a directory for %s\n",
	       curve_name);
    } /* end for all curves */
/* release the directory for all curves */
    Vdetach(curve_dir);
  } /* end if create dir */
  else
    printf("Error: save_curves_mappings: unable to create a directory for "
	   "the curves\n");

/* make a directory for the mappings */
  if ((mapping_dir = create_dir("mappings", "mappings", root)) != -1){
/* save the number of mappings */
    hput_int(xcog_data_ptr->n_grids, "n_grids", mapping_dir);
/* save each specific mapping */
    counter = 0;
    for (grid_ptr = xcog_data_ptr->first_m; grid_ptr != NULL;
	 grid_ptr = grid_ptr->next_m){
      sprintf(mapping_name, "mapping-%i", ++counter);
      if ((one_dir = create_dir(mapping_name, "one mapping", mapping_dir)) != -1){
/* save the generic and specific part of the data */
	write_mapping( one_dir, grid_ptr );
/* release the sub-directory */
	Vdetach(one_dir);
      }
      else
	printf("Error: save_curves_mappings: unable to create a directory for %s\n",
	       mapping_name);
    } /* end for all mappings */
/* release the directory for all mappings */
    Vdetach(mapping_dir);
  } /* end if create dir */
  else
    printf("Error: save_curves_mappings: unable to create a directory for "
	   "the mappings\n");

/* make a directory for the xcog parameters */
  if ((xcog_dir = create_dir("xcog_parameters", "xcog_parameters", root)) != -1){
/* save the overlap parameters */
    hput_int(xcog_data_ptr->disc_width, "disc_width", xcog_dir);
    hput_int(xcog_data_ptr->normal_width, "normal_width", xcog_dir);
    hput_int(xcog_data_ptr->tangent_width, "tangent_width", xcog_dir);
    hput_int(xcog_data_ptr->corner_width, "corner_width", xcog_dir);
    hput_int(xcog_data_ptr->interp_width, "interp_width", xcog_dir);
    hput_string((xcog_data_ptr->interp_type == 'e')? "explicit" : "implicit", 
		"interp_type", xcog_dir);
    hput_int(xcog_data_ptr->extra, "extra", xcog_dir);
    hput_int(xcog_data_ptr->trim_style, "trim_style", xcog_dir);

    hput_real(xcog_data_ptr->max_n_dist, "max_n_dist", xcog_dir);

    hput_int(xcog_data_ptr->show_physical_boundary, "show_physical_boundary", xcog_dir);
    hput_int(xcog_data_ptr->show_holes, "show_holes", xcog_dir);
    hput_int(xcog_data_ptr->correct_mismatch, "correct_mismatch", xcog_dir);

/* plotting stuff */
    hput_int(xcog_data_ptr->comp_plot_mode, "comp_plot_mode", xcog_dir);
    hput_int(xcog_data_ptr->grid_plot_mode, "grid_plot_mode", xcog_dir);
    hput_int(xcog_data_ptr->curve_plot_mode, "curve_plot_mode", xcog_dir);

/* bounding box */
    hput_real(xcog_data_ptr->xyab[0][0], "xyab00", xcog_dir);
    hput_real(xcog_data_ptr->xyab[1][0], "xyab10", xcog_dir);
    hput_real(xcog_data_ptr->xyab[0][1], "xyab01", xcog_dir);
    hput_real(xcog_data_ptr->xyab[1][1], "xyab11", xcog_dir);
  
/* title */
    hput_string(xcog_data_ptr->title, "title", xcog_dir);

/* local interpolation */
    hput_int(xcog_data_ptr->local_interp, "local_interp", xcog_dir);

/* release the directory */
    Vdetach(xcog_dir);
  }
  else
    printf("Error: save_curves_mappings: unable to create a directory for "
	   "the xcog parameters\n");
}


