#include "xcog.h"
#include "hdf_stuff.h"

/* static prototypes */
static char *
print_color(color_type grid_color);
static void 
delete_mapping( xcog_data *xcog_data_ptr, generic_mapping *grid_ptr );
static void 
delete_curve( xcog_data *xcog_data_ptr, generic_curve *curve_ptr );
static xcog_data *
reset_xcog( xcog_data *xcog_data_ptr, overlapping_grid *over_ptr );
static color_type
choose_color(input_output *io_ptr);
/* end static prototypes */


void 
main(int argc, char **argv){
  char *token, new_grid[80], *title[3], *file_name, *name, *statefile=NULL
    , *comfile=NULL;
  int icom, quit, replot, overlap_ok, f_number, tmp_fd, i, arg_len;
  int level=0, use_graphics=TRUE, i_min, i_max, j_min, j_max;
  int32 root;
  input_output io, *io_ptr=&io;
  xcog_data *xcog_data_ptr;
  generic_mapping *grid_ptr, *other_ptr, *new_grid_ptr;
  generic_curve *curve_ptr, *other_curve_ptr, *new_curve_ptr;
  overlapping_grid *over_ptr=NULL;
  FILE *fp;

  real xy[][2] = {{0.0, 1.0}, {0.67, 0.77}, {0.67, 0.37}};
  real width[]  = {0.65, 0.32, 0.32};
  real height[] = {0.6, 0.35, 0.35};

#include "xcog_com.h"

/* read the command line arguments */
  for (i = 1; i < argc; i++) {
    arg_len = strlen(argv[i]);
    if (!strncmp(argv[i], "-nographics", arg_len)) {
      use_graphics = FALSE;
    } 
    else if (!strncmp(argv[i], "-command", arg_len)){
      if (++i >= argc){
	printf("follow the -command option with a filename\n");
	exit(1);
      }
/* copy the filename for later use */
      comfile = (char *) malloc( (strlen(argv[i])+1) * sizeof(char) );
      strcpy( comfile, argv[i] );
    }
    else if (!strncmp(argv[i], "-state", arg_len)){
      if (++i >= argc){
	printf("follow the -state option with a filename\n");
	exit(1);
      }
/* copy the filename for later use */
      statefile = (char *) malloc( (strlen(argv[i])+1) * sizeof(char) );
      strcpy( statefile, argv[i] );
    }
    else{
      printf("Usage:\n"
	     "%s [-c[ommand] command_file] [-s[tate] state_file] [-n[ographics]]\n", 
	     argv[0]);
      exit(1);
    }
  }

/* start the interactive package */
  title[0] = "Composite grid (Xcog version " VERSION ")";
  title[1] = "Mapping (Xcog version " VERSION ")";
  title[2] = "Curve (Xcog version " VERSION ")";
  PL_init(io_ptr, title, xy, width, height, 3, use_graphics);

/* initialize the composite grid */
  xcog_data_ptr = new_xcog_data();

/* turn on the scale in all windows */
  PL_scale( 0, 1);
  PL_scale( 1, 1);
  PL_scale( 2, 1);

/* welcome message */
  welcome_to_xcog( );

  overlap_ok = 0;
  quit = 0;
  replot = 1;

/* try opening the file with the commands */
  if (comfile){
    if ((fp = open_this_ascii_file( comfile, 'r', 1, 1)) != NULL){
/* close the previous input file */
      if (io_ptr->read_command != stdin)
	fclose( io_ptr->read_command );
/* assign the new file as the input stream */
      io_ptr->read_command = fp;
    }
  }

/* try to read a saved state */
  if (statefile){
    if ( (root = open_hdf_file( statefile, 'r')) > 0 ){
      printf("Reading curves and mappings from the data-base... ");
      read_curves_mappings( root, xcog_data_ptr );
      printf("done\n");
      close_hdf_file( root );
/* compute new bounding box */
      set_global_boundaries( xcog_data_ptr );
      show_all( xcog_data_ptr );
      PL_window(0, xcog_data_ptr->xyab);
/* turn on grid boundaries and all curves */
      xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
      xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
/* update the plotting window */
      replot = TRUE;
    }
  }
  
/* open a log file where all commands will be saved */
/* iterate on finding a unique file name xcog-log-1 ... xcog-log-n */
  for (f_number = 1; ; f_number++){
    sprintf(new_grid, "xcog-log-%-i.com", f_number);
/* make sure the file does not exist */
    if ( (tmp_fd = open(new_grid, O_WRONLY, 0666)) == -1) break;
    close( tmp_fd );
  }
  if ( (fp = open_this_ascii_file( new_grid, 'w', 1, 0 )) != NULL){
    printf("All commands are beeing saved in the file `%s'\n", new_grid);
/* close the previous copy file */
    if (io.copy_command != NULL)
      fclose( io.copy_command );
/* assign the new file as the copy stream */
    io.copy_command = fp;
/* copy the file-name */
    io.copy_file_name = (char *) malloc( (strlen(new_grid)+1) * sizeof(char) );
    io.copy_file_name = strcpy( io.copy_file_name, new_grid );
  }


  do{

    if (replot){
/* plot the mapping */
      PL_erase(0);
      PL_start_plot(0);
/* plot the curves */
      for (curve_ptr = xcog_data_ptr->first_c; curve_ptr != NULL;
	   curve_ptr = curve_ptr->next)
	plot_curve( curve_ptr, xcog_data_ptr->curve_plot_mode );
/* plot the mappings */
      for (grid_ptr = xcog_data_ptr->first_m; grid_ptr != NULL;
	   grid_ptr = grid_ptr->next_m)
	plot_generic_mapping( grid_ptr, xcog_data_ptr->grid_plot_mode );
/* plot the composite grid */
      plot_overlapping( over_ptr, xcog_data_ptr->comp_plot_mode );
      PL_end_plot();
    }

    icom = get_command( io_ptr, "xcog>", COMMAND, LAST_COM, 
		       LEVEL, SAVE_ON_COPY, ARGUMENT);

/* the default is not to update the plotting window */
    replot = 0;

    switch (icom) {
    case MAKE_MAPPING:
/* make mapping */
      sprintf( new_grid, "grid-%i", xcog_data_ptr->n_grids+1);
      token = get_word( io_ptr, "New mapping name: ", new_grid, 1 );
/* check for uniqueness */
      for (other_ptr = xcog_data_ptr->first_m; other_ptr != NULL; 
	   other_ptr = other_ptr->next_m)
	if (strcmp(token, other_ptr->grid_name) == 0) break;
      if ( other_ptr != NULL ){
	printf("A mapping with that name already exits.\n");
      }
      else{
/* copy the name */
	name = (char *) malloc( (strlen(token)+1)*sizeof(char) );
	name = strcpy( name, token );
/* set the specific type of mapping */
	if ((grid_ptr = choose_mapping( io_ptr, name, xcog_data_ptr->first_c)) 
	    != NULL ){
/* insert in list */
	  insert_mapping( grid_ptr, xcog_data_ptr );
/* assign the default priority */
	  grid_ptr->priority = xcog_data_ptr->n_grids;
/* turn on grid boundaries */
	  xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
	  xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	  xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	  set_global_boundaries( xcog_data_ptr );
	  show_all( xcog_data_ptr );
	  PL_window(0, xcog_data_ptr->xyab);
/* update the plotting window */
	  replot = 1;
	}
/* free the space for the name */
	free( name );

      }
      break;

    case CHANGE_MAPPING:
/* change mapping */
      if (xcog_data_ptr->first_m == NULL){
	printf("The mapping list is empty.\n");
      }
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
/* change mapping */
	set_mapping( io_ptr, grid_ptr, xcog_data_ptr->first_c );
/* turn on grid boundaries */
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	set_global_boundaries( xcog_data_ptr );
	show_all( xcog_data_ptr );
	PL_window(0, xcog_data_ptr->xyab);
/* update the plotting window */
	replot = 1;
      }
      break;

    case COPY_MAPPING:
/* copy mapping */
      if (xcog_data_ptr->first_m == NULL){
	printf("The mapping list is empty.\n");
      }
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
	sprintf( new_grid, "grid-%i", xcog_data_ptr->n_grids+1);
	token = get_word( io_ptr, "New mapping name: ", new_grid, 1 );
/* check for uniqueness */
	for (other_ptr = xcog_data_ptr->first_m; other_ptr != NULL; 
	     other_ptr = other_ptr->next_m)
	  if (strcmp(token, other_ptr->grid_name) == 0) break;
	if ( other_ptr != NULL ){
	  printf("A mapping with that name already exits.\n");
	}
	else{
/* copy the name */
	  name = (char *) malloc( (strlen(token)+1)*sizeof(char) );
	  name = strcpy( name, token );
/* copy the mapping `grid_ptr' to `new_grid_ptr' and insert it into the list */
	  new_grid_ptr = copy_mapping( name, grid_ptr );
	  free( name );
/* insert in list */
	  insert_mapping( new_grid_ptr, xcog_data_ptr );
/* assign the default priority */
	  new_grid_ptr->priority = xcog_data_ptr->n_grids;
/* ask for rotations, translations or scalings */
	  set_transformation( io_ptr, new_grid_ptr );
/* turn on grid boundaries */
	  xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
	  xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	  xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	  set_global_boundaries( xcog_data_ptr );
	  show_all( xcog_data_ptr );
	  PL_window(0, xcog_data_ptr->xyab);
/* update the plotting window */
	  replot = 1;
	}
      }
      break;

    case TRANSFORM_MAPPING:
/* tranform mapping */
      if (xcog_data_ptr->first_m == NULL){
	printf("The mapping list is empty.\n");
      }
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
	set_transformation( io_ptr, grid_ptr );
/* turn on grid boundaries */
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	set_global_boundaries( xcog_data_ptr );
	show_all( xcog_data_ptr );
	PL_window(0, xcog_data_ptr->xyab);
	replot = 1;
      }
      break;

    case RENAME_MAPPING: 
/* rename mapping */
      if (xcog_data_ptr->first_m == NULL){
	printf("The mapping list is empty.\n");
      }
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
	token = get_word( io_ptr, "new mapping name: ", grid_ptr->grid_name, 1 ); 
/*  check for uniqueness  */
	for (other_ptr = xcog_data_ptr->first_m; other_ptr != NULL; 
	     other_ptr = other_ptr->next_m) 
	  if (strcmp(token, other_ptr->grid_name) == 0) break; 
	if ( other_ptr != NULL ){ 
	  printf("A mapping with that name already exits.\n"); 
	} 
	else{ 
/*  free the space for the previous name  */
	  free(grid_ptr->grid_name); 
/*  make space for the new name and assign it  */
	  grid_ptr->grid_name = (char *) malloc( (strlen(token)+1)*sizeof(char) );
	  grid_ptr->grid_name = strcpy( grid_ptr->grid_name, token ); 
	} 
      }
      break; 

    case DELETE_MAPPING:
/* delete mapping */
      if (xcog_data_ptr->first_m == NULL)
	printf("The mapping list is empty.\n");
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
	delete_mapping( xcog_data_ptr, grid_ptr);
/* turn on grid boundaries */
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	set_global_boundaries( xcog_data_ptr );
	show_all( xcog_data_ptr );
	PL_window(0, xcog_data_ptr->xyab);
/* update the plotting window */
	replot = 1;
      }
      break;

    case COLOR_MAPPING:
/* change color of a mapping */
      if (xcog_data_ptr->first_m == NULL)
	printf("The mapping list is empty.\n");
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
	printf("The mapping `%s' is presently %s\n", grid_ptr->grid_name, 
	       print_color(grid_ptr->grid_color));
	grid_ptr->grid_color = choose_color(io_ptr);
/* update the plotting window */
	replot = TRUE;
      }
      break;

    case MARK_MAPPING_BNDRY:
      if (xcog_data_ptr->first_m == NULL)
	printf("The mapping list is empty.\n");
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
	mark_mapping_bndry( io_ptr, grid_ptr );
/* turn off composite grid */
	xcog_data_ptr->comp_plot_mode = 0;
/* turn on grid boundaries */
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* update the plotting window */
	replot = TRUE;
      }
      break;

    case MAKE_CURVE:
/* make curve */
      sprintf( new_grid, "curve-%i", xcog_data_ptr->n_curves+1);
      token = get_word( io_ptr, "New curve name: ", new_grid, 1 );
/* check for uniqueness */
      for (curve_ptr = xcog_data_ptr->first_c; curve_ptr != NULL; 
	   curve_ptr = curve_ptr->next)
	if (strcmp(token, curve_ptr->curve_name) == 0) break;
      if ( curve_ptr != NULL ){
	printf("A curve with that name already exits.\n");
      }
      else{
/* copy the name */
	name = (char *) malloc( (strlen(token)+1)*sizeof(char) );
	name = strcpy( name, token );
/* set the specific type of curve */
	if ((curve_ptr=choose_curve( io_ptr, name)) != NULL ){
/* insert in curve list */
	  insert_curve( curve_ptr, xcog_data_ptr );
/* turn on grid boundaries */
	  xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
	  xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	  xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	  set_global_boundaries( xcog_data_ptr );
	  show_all( xcog_data_ptr );
	  PL_window(0, xcog_data_ptr->xyab);
/* update the plotting window */
	  replot = 1;
	}
	free( name );
      }
      break;

    case CHANGE_CURVE:
/* change curve */
      if (xcog_data_ptr->first_c == NULL)
	printf("The curve list is empty.\n");
      else if ((curve_ptr = get_curve_ptr( io_ptr, xcog_data_ptr->first_c )) != NULL){
	set_curve( io_ptr, curve_ptr );
/* turn on grid boundaries */
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	set_global_boundaries( xcog_data_ptr );
	show_all( xcog_data_ptr );
	PL_window(0, xcog_data_ptr->xyab);
/* update the plotting window */
	replot = 1;
      }
      break;

    case COPY_CURVE:
/* copy curve */
      if (xcog_data_ptr->first_c == NULL)
	printf("The curve list is empty.\n");
      else if ((curve_ptr = get_curve_ptr( io_ptr, xcog_data_ptr->first_c )) != NULL){
	sprintf( new_grid, "curve-%i", xcog_data_ptr->n_curves+1);
	token = get_word( io_ptr, "New curve name: ", new_grid, 1 );
/* check for uniqueness */
	for (other_curve_ptr = xcog_data_ptr->first_c; other_curve_ptr != NULL; 
	     other_curve_ptr = other_curve_ptr->next)
	  if (strcmp(token, other_curve_ptr->curve_name) == 0) break;
	if ( other_curve_ptr != NULL ){
	  printf("A curve with that name already exits.\n");
	}
	else{
/* copy the name */
	  name = (char *) malloc( (strlen(token)+1)*sizeof(char) );
	  name = strcpy( name, token );
/* copy the curve `curve_ptr' to `new_curve_ptr' and insert it into the list */
	  new_curve_ptr = copy_curve( name, curve_ptr );
	  free( name );
/* insert in list */
	  insert_curve( new_curve_ptr, xcog_data_ptr );
/* ask for rotations, translations or scalings */
	  set_curve_trans( io_ptr, new_curve_ptr );
/* turn on grid boundaries */
	  xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
	  xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	  xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	  set_global_boundaries( xcog_data_ptr );
	  show_all( xcog_data_ptr );
	  PL_window(0, xcog_data_ptr->xyab);
/* update the plotting window */
	  replot = 1;
	}
      }
      break;

    case TRANSFORM_CURVE:
/* transform curve */
      if (xcog_data_ptr->first_c == NULL)
	printf("The curve list is empty.\n");
      else if ((curve_ptr = get_curve_ptr( io_ptr, xcog_data_ptr->first_c )) != NULL){
	set_curve_trans( io_ptr, curve_ptr );
/* turn on grid boundaries */
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	set_global_boundaries( xcog_data_ptr );
	show_all( xcog_data_ptr );
	PL_window(0, xcog_data_ptr->xyab);
/* update the plotting window */
	replot = 1;
      }
      break;

    case RENAME_CURVE:
/* rename curve */
      if (xcog_data_ptr->first_c == NULL)
	printf("There are no curves to rename.\n");
      else if ((curve_ptr = get_curve_ptr( io_ptr, xcog_data_ptr->first_c )) != NULL){
	token = get_word( io_ptr, "new curve name: ", curve_ptr->curve_name, 1 );
	/* check for uniqueness */
	for (other_curve_ptr = xcog_data_ptr->first_c; other_curve_ptr != NULL; 
	     other_curve_ptr = other_curve_ptr->next)
	  if (strcmp(token, other_curve_ptr->curve_name) == 0) break;
	if ( other_curve_ptr != NULL ){
	  printf("A curve with that name already exits.\n");
	}
	else{
	  /* free the space for the previous name */
	  free(curve_ptr->curve_name);
	  /* make space for the new name and assign it */
	  curve_ptr->curve_name = (char *) malloc( (strlen(token)+1) * sizeof(char) );
	  curve_ptr->curve_name = strcpy( curve_ptr->curve_name, token );
	}
      }
      break;

    case DELETE_CURVE:
/* delete curve */
      if (xcog_data_ptr->first_c == NULL)
	printf("The curve list is empty.\n");
      else if ((curve_ptr = get_curve_ptr( io_ptr, xcog_data_ptr->first_c )) != NULL){
	if (curve_ptr->used_by_grid){
	  printf("This curve can't be deleted because it is used in the "
		 "definition\n");
	  printf("of a mapping.\n");
	}
	else{
	  delete_curve( xcog_data_ptr, curve_ptr);
/* turn on grid boundaries */
	  xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
	  xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	  xcog_data_ptr->comp_plot_mode = 0;
/* compute new bounding box */
	  set_global_boundaries( xcog_data_ptr );
	  show_all( xcog_data_ptr );
	  PL_window(0, xcog_data_ptr->xyab);
/* update the plotting window */
	  replot = 1;
	}
      }
      break;

    case BOUNDARY_CONDITION:
/* boundary condition */
      if (xcog_data_ptr->first_m == NULL){
	printf("The mapping list is empty.\n");
      }
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
	set_boundary_condition(io_ptr, grid_ptr);
/* redraw the figure if grid_plot_mode == 0 or if grid_plot_mode & 16 */
	replot = (xcog_data_ptr->grid_plot_mode == 0 || 
		  xcog_data_ptr->grid_plot_mode & 16)? 1 : 0;
/* turn on grid boundaries */
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	xcog_data_ptr->comp_plot_mode = 0;
      }
      break;

    case CURVE_LABEL:
/* curve label */
      if (xcog_data_ptr->first_m == NULL){
	printf("The mapping list is empty.\n");
      }
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) != NULL){
	set_curve_label(io_ptr, grid_ptr);
/* redraw the figure if grid_plot_mode == 0 or if grid_plot_mode & 32 */
/* or if grid_plot_mode & 1 */
	replot = (xcog_data_ptr->grid_plot_mode == 0 || 
		  xcog_data_ptr->grid_plot_mode & 32 ||
		  xcog_data_ptr->grid_plot_mode & 1 )? 1 : 0;
/* turn on grid boundaries */
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	xcog_data_ptr->comp_plot_mode = 0;
      }
      break;

    case GRID_LINES:
/* grid lines */
      if (xcog_data_ptr->first_m == NULL){
	printf("The mapping list is empty.\n");
      }
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
	set_gridlines(io_ptr, grid_ptr);
/* redraw the figure if grid_plot_mode == 0 or if grid_plot_mode & 4 or */
/* if grid_plot_mode & 8*/
	replot = (xcog_data_ptr->grid_plot_mode == 0 || 
		  xcog_data_ptr->grid_plot_mode & 4 ||
		  xcog_data_ptr->grid_plot_mode & 8)? 1 : 0;
/* turn on grid boundaries */
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
	xcog_data_ptr->comp_plot_mode = 0;
      }
      break;


    case GRID_QUALITY:
/* grid quality */
      if (xcog_data_ptr->first_m == NULL){
	printf("The mapping list is empty.\n");
      }
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){

/* plot the mapping with the directional arrows, all grid lines and quality */
/* information which is computed by a call to grid_quality from plot_generic_mapping */
	PL_window(1, grid_ptr->xyab);
	PL_erase(1);
	PL_start_plot(1);
	plot_generic_mapping( grid_ptr, 2+8+256 );
	PL_end_plot();
      }
      break;

    case OVERLAP_PARAMETERS:
/* overlap parameters */
      overlap_parameters( io_ptr, xcog_data_ptr);
      break;

    case LIST_OF_COMPONENTS:
/* change component grid list */
      component_grid_list( io_ptr, xcog_data_ptr );
/* turn on grid boundaries */
      xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
/* turn off composite grid */
      xcog_data_ptr->comp_plot_mode = 0;
      replot = TRUE;
      break;

    case COMPUTE_OVERLAP:
/* compute overlap */
      over_ptr = compute_overlap( io_ptr, xcog_data_ptr, over_ptr, &overlap_ok );
      if (overlap_ok == 1){
/* show the composite grid and the interpolation points */
	xcog_data_ptr->comp_plot_mode = 1 + 2;
/* turn off all grid lines */
	xcog_data_ptr->curve_plot_mode = 0;
	xcog_data_ptr->grid_plot_mode = 0;
/* update the plotting window */
	replot = 1;
/* print as hint */
	printf("\nThe overlapping grid algorithm succeeded.\n");
	printf("A CIRCLE indicates an interpolation point.\n"); 
      }
      else if (overlap_ok == 0){
/* show the composite grid and the bad points */
	xcog_data_ptr->comp_plot_mode = 1 + 4;
/* turn off all grid lines */
	xcog_data_ptr->curve_plot_mode = 0;
	xcog_data_ptr->grid_plot_mode = 0;
/* update the plotting window */
	replot = 1;
    printf("\n"
"Error: The overlapping grid algorithm failed. Follow the following steps\n"
"to improve the situation:\n\n"
"o Make sure that all parts of the physical domain are covered by grids.\n"
"o Make sure there are sufficiently many grid points in the overlap regions.\n"
"o Make sure that a curve-label has been assigned to the sides of the grids \n"
"  that are aligned with the physical boundary. Sides that are completely aligned \n"
"  with the boundary should be given a positive curve-label, and sides that \n"
"  are partly on the boundary, and partly inside of another grid, should be \n"
"  given a negative curve-label. The absolute value of the curve-label should \n"
"  be the same for all sides of all grids that are aligned with the same part \n"
"  of the boundary. Hence, if the domain is simply connected, only one \n"
"  curve-label value should be used. For a doubly connected domain, two \n"
"  values should be used, and so on.\n"
"\n");
/* print a hint */
	printf("\n"
"The grid points that didn't satisfy the requirements for a discretization point or\n"
"an interpolation point are marked in the following way:\n\n"
" BOX               Bad discretization point.\n"
" CROSS             Dead interpolation location.\n"
" FILLED CIRCLE     Non-explicit interpolation location.\n"
" TRIANGLE          Bad interpolation point because of a dead interpolation location.\n"
" DOWNWARD TRIANGLE Bad interpolation point because of a non-explicit interpolation\n"
"                   location.\n");
      }
      else{
/* cancel: just update the plotting window */
	replot = 1;
      }

      break;

    case RE_COMPUTE:
/* re-compute overlap */
      if (over_ptr != NULL) 
	overlap_ok = overlap(io_ptr, over_ptr, TRUE, FALSE, FALSE); 
/* the holes are not */
/* cut again, so it does not make sense to look at the physical boundary again */
      else{
	printf("This command requires an existing overlapping grid, but no such grid "
	       "is currently defined.\n");
	break;
      }

      if (overlap_ok == 1){
/* show the composite grid and the interpolation points */
	xcog_data_ptr->comp_plot_mode = 1 + 2;
/* turn off all grid lines */
	xcog_data_ptr->curve_plot_mode = 0;
	xcog_data_ptr->grid_plot_mode = 0;
/* update the plotting window */
	replot = 1;
/* print as hint */
	printf("\nThe overlapping grid algorithm succeeded.\n");
	printf("A CIRCLE indicates an interpolation point.\n"); 
      }
      else if (overlap_ok == 0){
/* show the composite grid and the bad points */
	xcog_data_ptr->comp_plot_mode = 1 + 4;
/* turn off all grid lines */
	xcog_data_ptr->curve_plot_mode = 0;
	xcog_data_ptr->grid_plot_mode = 0;
/* update the plotting window */
	replot = 1;
    printf("\n"
"Error: The overlapping grid algorithm failed. Follow the following steps\n"
"to improve the situation:\n\n"
"o Make sure that all parts of the physical domain are covered by grids.\n"
"o Make sure there are sufficiently many grid points in the overlap regions.\n"
"o Make sure the boundary condition is non-zero on all physical boundaries.\n"
"o Make sure the curve-label has the same non-zero value when the same physical\n"
"  boundary is represented by more than one grid.\n"
"o Make sure that all external corners are identified.\n\n");
/* print a hint */
	printf("\n"
"The grid points that didn't satisfy the requirements for a discretization point or\n"
"an interpolation point are marked in the following way:\n\n"
" BOX               Bad discretization point.\n"
" CROSS             Dead interpolation location.\n"
" FILLED CIRCLE     Non-explicit interpolation location.\n"
" TRIANGLE          Bad interpolation point because of a dead interpolation location.\n"
" DOWNWARD TRIANGLE Bad interpolation point because of a non-explicit interpolation\n"
"                   location.\n");
      }
      else{
/* cancel: just update the plotting window */
	replot = 1;
      }

      break;

    case SAVE_COMPONENT:
/* save a component grid in PLOT3D format */
      if (xcog_data_ptr->first_m == NULL){
	printf("The mapping list is empty.\n");
      }
      else if ((grid_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) 
	       != NULL){
/* open an ascii file and save the component grid */
	if ( (fp = open_ascii_file(io_ptr, "component grid file name: ",
				   "test.plot3d", &file_name, 'w', 0, 
				   SAVE_ON_COPY[SAVE_COMPONENT])) != NULL){
	  i_min = int_max( 1, int_min(grid_ptr->r_points, 
				      get_int(io_ptr, "Starting i-index: ", 
					      1, LEVEL+1) ) );
	  i_max = int_max( 1, int_min(grid_ptr->r_points, 
				      get_int(io_ptr, "Ending i-index: ", 
					      grid_ptr->r_points, LEVEL+1) ) );
	  j_min = int_max( 1, int_min(grid_ptr->s_points, 
				      get_int(io_ptr, "Starting j-index: ", 
					      1, LEVEL+1) ) );
	  j_max = int_max( 1, int_min(grid_ptr->s_points, 
				      get_int(io_ptr, "Ending j-index: ", 
					      grid_ptr->s_points, LEVEL+1) ) );
	  save_grid_plot3d( grid_ptr, fp, i_min, i_max, j_min, j_max );
	  fclose( fp );
	}
      }
      break;

    case SAVE_OVERLAPPING:
/* save overlapping */
      if ( overlap_ok != 0 && over_ptr != NULL ){
	save_xcog_data( io_ptr, over_ptr );
      }
      else
	printf("There is no valid overlapping grid to save.\n");

      break;

    case SHOW:
/* show */
      if (io.copy_command != NULL)
	printf("All commands are currently beeing saved on file `%s'.\n",
	       io.copy_file_name);
      else
	printf("The commands are not beeing saved on file.\n");
/* parameters */
      print_mapping_list(xcog_data_ptr);
      print_curve_list(xcog_data_ptr);

      break;

    case PLOT_MODE:
/* plot mode */
      set_global_plot_mode( io_ptr, xcog_data_ptr, over_ptr );
      break;

    case READ_COMMAND: 
/* read command-file */
      if ( (fp = open_ascii_file( io_ptr, "command file name: ",
				 "intro.com", &file_name, 'r', 1, 
				 SAVE_ON_COPY[READ_COMMAND])) != NULL){
/* close the previous input file */
	if (io.read_command != stdin)
	  fclose( io.read_command );
/* assign the new file as the input stream */
	io.read_command = fp;
      }

      break;

    case START_SAVING:
/* start saving-commands */
      

      if (io.copy_command != NULL){
      }

/* try to open the new copy file */
      if ( (fp = open_ascii_file( io_ptr, "Save commands in file named: ", 
				 "commands.com", &file_name, 'w', 1,
				 SAVE_ON_COPY[START_SAVING])) != NULL){
/* close the previous copy file */
	if (io.copy_command != NULL){
/* remind the user of the file name */
	  printf("Closing the previous command file `%s'\n", io.copy_file_name);
/* close the file */
	  fclose( io.copy_command );
	  io.copy_command = NULL;
/* free the space occupied by the file-name */
	  free( io.copy_file_name );
	}

/* assign the new file as the copy stream */
	io.copy_command = fp;
/* copy the file-name */
	if ((io.copy_file_name = (char *) 
	     malloc( (strlen(file_name)+1) * sizeof(char) )) == NULL)
	  printf("memory error in xcog, copy filename\n");
	io.copy_file_name = strcpy( io.copy_file_name, file_name );
      }

      break;

    case STOP_SAVING:
/* stop saving-commands */
      if (io.copy_command != NULL){
/* remind the user of the file name */
	printf("Closing command file `%s'\n", io.copy_file_name);
/* close the file */
	fclose( io.copy_command );
	io.copy_command = NULL;
/* free the space occupied by the file-name */
	free( io.copy_file_name );
      }

      break;

    case RESET_XCOG:
/* reset */
      xcog_data_ptr = reset_xcog( xcog_data_ptr, over_ptr );
      xcog_data_ptr->grid_plot_mode = 0;
      xcog_data_ptr->curve_plot_mode = 0;
      overlap_ok = 0;
      over_ptr = NULL;
/* compute new bounding box */
      set_global_boundaries( xcog_data_ptr );
      show_all( xcog_data_ptr );
      PL_window(0, xcog_data_ptr->xyab);
      PL_erase(0);
      PL_erase(1);
      PL_erase(2);
/* turn on the scale in all windows */
      PL_scale( 0, 1);
      PL_scale( 1, 1);
      PL_scale( 2, 1);
/* update the plotting window */
      replot = 1;

      break;

    case READ_STATE:
/* read saved state */
      if ((xcog_data_ptr->n_grids > 0 || xcog_data_ptr->n_curves > 0) &&
	  get_yes_no( io_ptr, "Do you want to reset xcog before reading the "
		     "saved state>", level )){
	xcog_data_ptr = reset_xcog( xcog_data_ptr, over_ptr );
	xcog_data_ptr->grid_plot_mode = 0;
	xcog_data_ptr->curve_plot_mode = 0;
	overlap_ok = 0;
	over_ptr = NULL;
      }
      name = get_word( io_ptr, "File name for reading curves and mappings: ",
		      "test.xcog", 2);
      if ( (root = open_hdf_file( name, 'r' )) > 0 ){
	read_curves_mappings( root, xcog_data_ptr );
	close_hdf_file( root );
/* update the view port according to the new bounding box */
	PL_window(0, xcog_data_ptr->xyab);
/* turn on grid boundaries and all curves because right now, we don't save the */
/* component grid, so if only comp_plot != 0, nothing will be shown on screen */
	xcog_data_ptr->grid_plot_mode = xcog_data_ptr->grid_plot_mode | 1;
	xcog_data_ptr->curve_plot_mode = xcog_data_ptr->curve_plot_mode | 1;
/* update the plotting window */
	replot = 1;
      }

      break;

    case SAVE_STATE:
/* save current state */
      name = get_word( io_ptr, "File name for saving all curves and mappings: ", 
		      "test.xcog", 2);
      if ( (root = open_hdf_file( name, 'i')) > 0){
	save_curves_mappings( root, xcog_data_ptr );
	close_hdf_file( root );
      }
      else{
	printf("Unable to open the database file %s\n", name);
      }

      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );

      break;

    case TUTORIAL:
      run_tutorial( io_ptr );
      break;

    case QUIT:
/* end */
/* stop saving-commands */
      if (io.copy_command != NULL){
/* remind the user of the file name */
	printf("Closing command file `%s'\n", io.copy_file_name);
/* close the file */
	fclose( io.copy_command );
	io.copy_command = NULL;
/* free the space occupied by the file-name */
	free( io.copy_file_name );
      }
/* close the connection to X */
      PL_reset(io_ptr);
/* clean the heap */
      xcog_data_ptr = reset_xcog( xcog_data_ptr, over_ptr );
      free(xcog_data_ptr);

      printf("Good bye\n");
      quit = 1;

      break;

    default:
      ;
    }
  }
  while(!quit);

}


static xcog_data *
reset_xcog( xcog_data *xcog_data_ptr, overlapping_grid *over_ptr ){
  generic_mapping *grid_ptr, *next_grid_ptr;
  generic_curve *curve_ptr, *next_curve_ptr;

/* delete the overlapping grid */
  delete_overlapping_grid( over_ptr );

/* delete all mappings */
  grid_ptr = xcog_data_ptr->first_m;
  while (grid_ptr != NULL){
    next_grid_ptr = grid_ptr->next_m;
/* this routine clobbers the heap */
    delete_mapping( xcog_data_ptr, grid_ptr);
    grid_ptr = next_grid_ptr;
  }

/* delete all curves */
  curve_ptr = xcog_data_ptr->first_c;
  while (curve_ptr != NULL){
    next_curve_ptr = curve_ptr->next;
    delete_curve( xcog_data_ptr, curve_ptr);
    curve_ptr = next_curve_ptr;
  }

/* delete the structure itself */
  free( xcog_data_ptr );

/* make a new composite grid */
  return new_xcog_data();
}


static void 
delete_mapping( xcog_data *xcog_data_ptr, generic_mapping *grid_ptr ){
/* update the list of mappings */
  if (grid_ptr->prev_m != NULL)
    grid_ptr->prev_m->next_m = grid_ptr->next_m;
  if (grid_ptr->next_m != NULL)
    grid_ptr->next_m->prev_m = grid_ptr->prev_m;
  if (xcog_data_ptr->first_m == grid_ptr)
    xcog_data_ptr->first_m = grid_ptr->next_m;
  if (xcog_data_ptr->last_m == grid_ptr)
    xcog_data_ptr->last_m = grid_ptr->prev_m;

/* decrease the count of mappings */
  xcog_data_ptr->n_grids--;

  delete_generic_mapping( grid_ptr );
}


static void 
delete_curve( xcog_data *xcog_data_ptr, generic_curve *curve_ptr ){
/* update the list */
  if (curve_ptr->prev != NULL)
    curve_ptr->prev->next = curve_ptr->next;
  if (curve_ptr->next != NULL)
    curve_ptr->next->prev = curve_ptr->prev;
  if (xcog_data_ptr->first_c == curve_ptr)
    xcog_data_ptr->first_c = curve_ptr->next;
  if (xcog_data_ptr->last_c == curve_ptr)
    xcog_data_ptr->last_c = curve_ptr->prev;

/* decrease the count of mappings */
  xcog_data_ptr->n_curves--;

  delete_generic_curve( curve_ptr );
}


static char *
print_color(color_type grid_color){
  enum{RED, GREEN, BLUE, BROWN, PURPLE, BLACK, GREY};
  real min_dist, dist[7], colors[7][3] = {
    {255.0, 0.0, 0.0}, {50.0, 205.0, 50.0}, {0.0, 0.0, 255.0},
    {190.0, 81.0, 6.0}, {240.0, 26.0, 207.0}, {0.0, 0.0, 0.0},
    {155.0, 166.0, 187.0}};
  int i, j, color_index;
/* normalize */
  for (i=0; i<7; i++)
    for (j=0; j<3; j++)
      colors[i][j] /= 255.0;

/* compute the distance */
  for (i=0; i<7; i++){
    dist[i] = fabs(grid_color.r - colors[i][0]) + 
      fabs(grid_color.g - colors[i][1]) + 
	fabs(grid_color.b - colors[i][2]);
  }

/* search for the smallest distance */
  min_dist = dist[0]; color_index = 0;
  for (i=1; i<7; i++)
    if (dist[i] < min_dist){
      min_dist = dist[i];
      color_index = i;
    }

/* return the appropriate string */
  switch(color_index){
  case RED:
    return "red";
  case GREEN:
    return "green";
  case BLUE:
    return "blue";
  case BROWN:
    return "brown";
  case PURPLE:
    return "purple";
  case BLACK:
    return "black";
  case GREY:
    return "grey";
  default:
    return "unknown color";
  } /* end switch */
}

static color_type
choose_color(input_output *io_ptr){
  color_type color;
#include "choose_color_com.h"
  switch (get_command(io_ptr, "Choose color>", COMMAND, LAST_COM, LEVEL, 
		     SAVE_ON_COPY, ARGUMENT)){
  case RED:
/* red */
    color.r = 255.0;
    color.g = 0.0;
    color.b = 0.0;
    break;

  case GREEN:
/* green */
    color.r = 50.0;
    color.g = 205.0;
    color.b = 50.0;
    break;

  case BLUE:
/* blue */
    color.r = 0.0;
    color.g = 0.0;
    color.b = 255.0;
    break;

  case BROWN:
/* brown */
    color.r = 190.0;
    color.g = 81.0;
    color.b = 6.0;
    break;

  case PURPLE:
/* purple */
    color.r = 240.0;
    color.g = 26.0;
    color.b = 207.0;
    break;

  case BLACK:
/* black */
    color.r = 0.0;
    color.g = 0.0;
    color.b = 0.0;
    break;

  case GREY:
  default:
/* grey */
    color.r = 155.0;
    color.g = 166.0;
    color.b = 187.0;
    break;

  }

  color.r = color.r/255.0;
  color.g = color.g/255.0;
  color.b = color.b/255.0;

  return color;
}
