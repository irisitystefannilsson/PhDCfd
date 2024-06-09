#include "xcog.h"

/* static prototypes */
static void 
copy_overlapping( overlapping_grid *over_ptr, xcog_data *xcog_data_ptr );
static void 
copy_grid( component_grid *grid_ptr, generic_mapping *map_ptr, int local_interp );
static generic_mapping *
find_mapping_priority(xcog_data *xcog_data_ptr, int priority );
static void 
init_grid_arrays(component_grid *grid_ptr, generic_mapping *map_ptr,
		 overlapping_grid *over_ptr);
/* end static prototypes */


void 
component_grid_list( input_output *io_ptr, xcog_data *xcog_data_ptr ){
  generic_mapping *map_ptr;
  int priority, replot;
  const int plot_mode=1;
  char *prompt;
  int level=1, quit, icom, n_grids, max_priority;

#include "component_grid_list_com.h"

/* find the max priority */
  max_priority = 0;
  for (map_ptr = xcog_data_ptr->first_m; map_ptr != NULL; 
       map_ptr = map_ptr->next_m){
    if (map_ptr->priority > max_priority) max_priority = map_ptr->priority;
  }

/* Insert all the grids with positive priority */

  n_grids = 0;
  for (priority = 1; priority <= max_priority; priority++){
    map_ptr = find_mapping_priority( xcog_data_ptr, priority );
    if (map_ptr != NULL){
      n_grids++;
      map_ptr->priority = n_grids;
    }
  }

/* ask the user for modifications of the list, i.e. the priority of each mapping */

  prompt = "component grid list>";

  quit = 0;
  replot = (xcog_data_ptr->grid_plot_mode == plot_mode && 
	    n_grids == xcog_data_ptr->n_grids)? 0 : 1;
  do{

    if (replot){
/* plot the mapping */
      PL_erase(0);
      PL_start_plot(0);
      for (map_ptr = xcog_data_ptr->first_m; map_ptr != NULL;
           map_ptr = map_ptr->next_m)
/* only plot the grids with positive priority */
	if (map_ptr->priority > 0)
	  plot_generic_mapping( map_ptr, plot_mode );
      PL_end_plot();
    }

    icom = get_command(io_ptr, prompt, command, ncom, level,
		       save_on_copy, argument);

/* default is plot the screen after the command */
    replot = 1;
    switch (icom) {

    case 0: 
/* insert all mappings */
      if (n_grids == xcog_data_ptr->n_grids){
	printf("All mappings are already in the list.\n");
	replot = 0;
      }
      else{
	priority = xcog_data_ptr->n_grids;
	for (map_ptr = xcog_data_ptr->first_m; map_ptr != NULL; 
	     map_ptr = map_ptr->next_m){
	  map_ptr->priority = priority;
	  priority--;
	} 
	n_grids = xcog_data_ptr->n_grids;
      }
      
      break;

    case 1:
/* insert first */
      if ( (map_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) != NULL){
/* check if this grid is already in the list */
	if ( map_ptr->priority > 0 )
 	  printf("That mapping is already in the list\n"); 
 	else{ 
/* insert that mapping first */
	  n_grids++;
	  map_ptr->priority = n_grids;
	}
      }
      break;

    case 2:
/* remove-from-list */
      if ( (map_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) != NULL){
/* check if this grid is in the list */
	if ( map_ptr->priority == 0 )
 	  printf("That mapping is not in the list\n"); 
 	else{ 
/* lower the priority of all mappings with priority higher than this one */
	  for (priority = map_ptr->priority+1; priority <= n_grids; priority++){
	    find_mapping_priority( xcog_data_ptr, priority )->priority--;
	  }
/* delete that mapping from the list*/
	  n_grids--;
	  map_ptr->priority = 0;
	}
      }

      break;

    case 3:
/* move-first */
      if ( (map_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) != NULL){
	if (map_ptr->priority == n_grids)
	  printf("That mapping is already first\n");
	else if (map_ptr->priority == 0)
	  printf("That mapping is not in the list\n");
	else{
/* update the list */
	  for (priority=map_ptr->priority+1; priority<=n_grids; priority++)
	    find_mapping_priority( xcog_data_ptr, priority )->priority--;
/* move the present element last */
	  map_ptr->priority = n_grids;
	}
      }
      replot = 0;
      break;

    case 4:
/* move-last */
      if ( (map_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) != NULL){
	if (map_ptr->priority == 1)
	  printf("That mapping is already last\n");
	else if (map_ptr->priority == 0)
	  printf("That mapping is not in the list\n");
	else{
/* update the list */
	  for (priority=map_ptr->priority-1; priority>=1; priority--)
	    find_mapping_priority( xcog_data_ptr, priority )->priority++;
/* move the present element last */
	  map_ptr->priority = 1;
	}
      }
      replot = 0;
      break;

    case 5:
/* move up */
      if ( (map_ptr = get_grid_ptr( io_ptr, xcog_data_ptr->first_m )) != NULL){
	if (map_ptr->priority == n_grids)
	  printf("That mapping is already first\n");
	else if (map_ptr->priority == 0)
	  printf("That mapping is not in the list\n");
	else{
/* swap this and the one above it */
	  find_mapping_priority( xcog_data_ptr, map_ptr->priority+1 )->priority--;
	  map_ptr->priority++;
	}
      }
      replot = 0;
      break;

    case 6:
/* show mapping list */
      if (n_grids == 0)
	printf("The mapping list is empty\n");
      else{
	printf("List of mappings to be used in the composite grid.\n");
	for (priority = n_grids; priority >= 1; priority--){
	  map_ptr = find_mapping_priority( xcog_data_ptr, priority );
	  printf("%i: %s\n", map_ptr->priority, map_ptr->grid_name);
	}
	printf("\n");
      }
      replot = 0;
      break;

    case 7:
/* help */
      while ( (icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				  command, ncom, level+1, NULL, NULL)) == -1 );
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case 8:
/* exit */
      quit = 1;
      break;

    default:
      replot = 0;
      ;

    } 

  } while(!quit);

  return;
}

overlapping_grid *
compute_overlap(input_output *io_ptr, xcog_data *xcog_data_ptr, 
		overlapping_grid *over_ptr, int *overlap_ok){
  generic_mapping *map_ptr;
  int priority, total_grid_points, total_n_interp, total_hole_points, hole_points
    , i, j;
  component_grid *grid_ptr;
  int n_grids, max_priority;

/* Insert all the grids with positive priority */

/* find the max priority */
  max_priority = 0;
  for (map_ptr = xcog_data_ptr->first_m; map_ptr != NULL; 
       map_ptr = map_ptr->next_m){
    if (map_ptr->priority > max_priority) max_priority = map_ptr->priority;
  }

  n_grids = 0;
  for (priority = 1; priority <= max_priority; priority++){
    map_ptr = find_mapping_priority( xcog_data_ptr, priority );
    if (map_ptr != NULL){
      n_grids++;
      map_ptr->priority = n_grids;
    }
  }

/* check if there are any mappings left */
  if (n_grids == 0){
    printf("There are no mappings left to construct an overlapping grid from.\n");
    *overlap_ok = -1;
    return over_ptr;
  }

/* Copy the data to the overlapping grid structure */
  printf("Evaluating the mapping functions at the grid points and initializing the\n"
	 "data structure...");
  fflush(stdout);

/* remove any old overlapping grid */
  if (over_ptr != NULL) delete_overlapping_grid( over_ptr );
/* make a new structure */
  over_ptr = new_overlapping_grid( n_grids );
  copy_overlapping( over_ptr, xcog_data_ptr );

  for (priority = 1; priority <= n_grids; priority++){
 
    map_ptr = find_mapping_priority( xcog_data_ptr, priority );
/* create the component grid structure */
    grid_ptr = new_component_grid( map_ptr->grid_name );
/* insert it into the list */
    insert_grid( grid_ptr, over_ptr );
/* copy sizes and pointers from the mapping */
    copy_grid( grid_ptr, map_ptr, xcog_data_ptr->local_interp );
/* allocate space for the large arrays and initialize them */
    init_grid_arrays( grid_ptr, map_ptr, over_ptr );
  }


  printf("Done\n");

/* compute the overlap */
/* overlap_ok == 1 indicates sucess, 0 indicates error and -1 means cancel */
  *overlap_ok = overlap(io_ptr, over_ptr, FALSE, xcog_data_ptr->show_physical_boundary, 
			xcog_data_ptr->show_holes );

/* print grid stats */
  if (*overlap_ok == 1){
    total_grid_points = 0;
    total_n_interp = 0;
    total_hole_points = 0;
    printf("\n# grid points  # interpolation points  # hole points   name \n");
    for (grid_ptr = over_ptr->first; grid_ptr != NULL;
	 grid_ptr = grid_ptr->next){
/* count the number of hole-points by looking at the flag array */
      hole_points = 0;
      for (i=range(1,1); i<=range(2,1); i++)
	for (j=range(1,2); j<=range(2,2); j++)
	  if (flag(i,j) == 0) hole_points++;
/* print the result for this grid */
      printf("%-14i %-23i %-15i %s\n", 
	     grid_ptr->r_points * grid_ptr->s_points, 
	     grid_ptr->n_interp, hole_points, grid_ptr->grid_name); 
      total_grid_points += grid_ptr->r_points * grid_ptr->s_points;
      total_n_interp += grid_ptr->n_interp;
      total_hole_points += hole_points;
    }
    printf("_______________________________________________________________________\n");
    printf("%-14i %-23i %-15i %s\n", total_grid_points, total_n_interp, 
	   total_hole_points, "overlapping grid");
  }
  return over_ptr;
}

static generic_mapping *
find_mapping_priority( xcog_data *xcog_data_ptr, int priority ){
  generic_mapping *map_ptr;
/* find the mapping with that priority */
  for (map_ptr = xcog_data_ptr->last_m; map_ptr != NULL; 
       map_ptr = map_ptr->prev_m)
    if (map_ptr->priority == priority) break;
  return map_ptr;
}

static void 
copy_overlapping( overlapping_grid *over_ptr, xcog_data *xcog_data_ptr ){

  over_ptr->extra = xcog_data_ptr->extra;

  over_ptr->disc_width = xcog_data_ptr->disc_width;
  over_ptr->normal_width = xcog_data_ptr->normal_width;
  over_ptr->tangent_width = xcog_data_ptr->tangent_width;
  over_ptr->corner_width = xcog_data_ptr->corner_width;

  over_ptr->interp_width = xcog_data_ptr->interp_width;
  over_ptr->interp_type = xcog_data_ptr->interp_type;

  over_ptr->trim_style = xcog_data_ptr->trim_style;

  over_ptr->x_min = xcog_data_ptr->x_min;
  over_ptr->x_max = xcog_data_ptr->x_max;
  over_ptr->y_min = xcog_data_ptr->y_min;
  over_ptr->y_max = xcog_data_ptr->y_max;

  over_ptr->max_n_dist = xcog_data_ptr->max_n_dist;
  over_ptr->correct_mismatch = xcog_data_ptr->correct_mismatch;
}


static void 
copy_grid( component_grid *grid_ptr, generic_mapping *map_ptr, int local_interp ){
  int i, j;

#define g_range(i, j) compute_index_2d( grid_ptr->range_ptr, i, j)
#define g_bc(i, j)    compute_index_2d( grid_ptr->bc_ptr, i, j)
#define g_curve(i, j) compute_index_2d( grid_ptr->curve_ptr, i, j)
#define g_cut(i, j)   compute_index_2d( grid_ptr->cut_ptr, i, j)

#define m_range(i, j) compute_index_2d( map_ptr->range_ptr, i, j)
#define m_bc(i, j)    compute_index_2d( map_ptr->bc_ptr, i, j)
#define m_curve(i, j) compute_index_2d( map_ptr->curve_ptr, i, j)
#define m_cut(i, j)   compute_index_2d( map_ptr->cut_ptr, i, j)

  grid_ptr->grid_type = map_ptr->grid_type;

  grid_ptr->priority = map_ptr->priority;

  grid_ptr->inverse_known = map_ptr->inverse_known;
  grid_ptr->grid_mapping = map_ptr->grid_mapping;
  grid_ptr->inverse_grid_mapping = map_ptr->inverse_grid_mapping;
  grid_ptr->mapping_data_ptr = map_ptr->mapping_data_ptr;

  grid_ptr->r_points = map_ptr->r_points;
  grid_ptr->s_points = map_ptr->s_points;

  grid_ptr->r_period = map_ptr->r_period;
  grid_ptr->s_period = map_ptr->s_period;

  grid_ptr->r_step = 1.0/(grid_ptr->r_points - 1);
  grid_ptr->s_step = 1.0/(grid_ptr->s_points - 1);

  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++){
      g_bc(i,j) = m_bc(i,j);
      g_curve(i,j) = m_curve(i,j);
    }

/* don't use the forwards mapping if local_interp is on */
  if (local_interp)
    grid_ptr->grid_mapping = NULL;

  grid_ptr->grid_color.r = map_ptr->grid_color.r;
  grid_ptr->grid_color.g = map_ptr->grid_color.g;
  grid_ptr->grid_color.b = map_ptr->grid_color.b;

/* rotation, translation and scaling */
  grid_ptr->cos_theta = map_ptr->cos_theta;
  grid_ptr->sin_theta = map_ptr->sin_theta;
  grid_ptr->x_trans = map_ptr->x_trans;
  grid_ptr->y_trans = map_ptr->y_trans;
  grid_ptr->scaling = map_ptr->scaling;

/* external gaps */
  for (i=0; i<2; i++){
    grid_ptr->gap_low_r[i]  = map_ptr->gap_low_r[i];
    grid_ptr->gap_high_r[i] = map_ptr->gap_high_r[i];
    grid_ptr->gap_low_s[i]  = map_ptr->gap_low_s[i];
    grid_ptr->gap_high_s[i] = map_ptr->gap_high_s[i];
  }

#undef g_range
#undef g_bc
#undef g_curve

#undef m_range
#undef m_bc
#undef m_curve
}

static void 
init_grid_arrays(component_grid *grid_ptr, generic_mapping *map_ptr, 
		 overlapping_grid *over_ptr){
  int i,j,ierr=0;
  grid_point gp;
  hyp_grid_data *info;

/* number of extra points at periodic boundaries */
  over_ptr->extra_period = int_max((over_ptr->disc_width - 1)/2, 
				   (over_ptr->interp_width - 1)/2);
/* the next line saves the run if disc_width = 1 and interp_width = 1 */
  over_ptr->extra_period = int_max(1, over_ptr->extra_period);

  if (grid_ptr->r_period){
    range(1,1) = over_ptr->extra_period + 1; 
    range(2,1) = grid_ptr->r_points + over_ptr->extra_period;
    grid_ptr->r_dim = grid_ptr->r_points + 2*over_ptr->extra_period - 1;
  }
  else{
    range(1,1) = over_ptr->extra + 1; 
    range(2,1) = grid_ptr->r_points + over_ptr->extra;
    grid_ptr->r_dim = grid_ptr->r_points + 2*over_ptr->extra;
  }

  if (grid_ptr->s_period){
    range(1,2) = over_ptr->extra_period + 1; 
    range(2,2) = grid_ptr->s_points + over_ptr->extra_period;
    grid_ptr->s_dim = grid_ptr->s_points + 2*over_ptr->extra_period - 1;
  }
  else{
    range(1,2) = over_ptr->extra + 1; 
    range(2,2) = grid_ptr->s_points + over_ptr->extra;
    grid_ptr->s_dim = grid_ptr->s_points + 2*over_ptr->extra;
  }

  grid_ptr->r_step = 1/((real) (grid_ptr->r_points-1));
  grid_ptr->s_step = 1/((real) (grid_ptr->s_points-1));


/* hyperbolic c-grids */
  if (grid_ptr->grid_type == 8){
    info = (hyp_grid_data *) map_ptr->mapping_data_ptr;
    if (info->c_grid){
      grid_ptr->topology = 1;
/* right now, the branch-cut is along the r-direction */
      grid_ptr->branch_direction = 1;
      grid_ptr->first_point = range(1,1) + info->n_cut;
      grid_ptr->last_point  = grid_ptr->first_point + info->n_curve - 1;
    }
  }

/* define the boundary classification arrays */
  grid_ptr->left_interp_ptr = create_int_array_1d( grid_ptr->s_dim );
  grid_ptr->right_interp_ptr = create_int_array_1d( grid_ptr->s_dim );
  grid_ptr->lower_interp_ptr = create_int_array_1d( grid_ptr->r_dim );
  grid_ptr->upper_interp_ptr = create_int_array_1d( grid_ptr->r_dim );

/* assign default values */
  for (i=1; i<= grid_ptr->r_dim; i++){
    lower_interp(i) = 0;
    upper_interp(i) = 0;
  }
  for (j=1; j<= grid_ptr->s_dim; j++){
    left_interp(j)  = 0;
    right_interp(j) = 0;
  }

/* define the multi-dimensional arrays */

  grid_ptr->x_ptr = create_real_array_2d(grid_ptr->r_dim,grid_ptr->s_dim);
  grid_ptr->y_ptr = create_real_array_2d(grid_ptr->r_dim,grid_ptr->s_dim);

  grid_ptr->flag_ptr = create_int_array_2d(grid_ptr->r_dim,grid_ptr->s_dim);

/* initialize bounding box */
  grid_ptr->x_min =  1.e10;
  grid_ptr->x_max = -1.e10;
  grid_ptr->y_min =  1.e10;
  grid_ptr->y_max = -1.e10;

/* assign x, y, xr, yr, xs and ys. */
  for (i = 1; i<= grid_ptr->r_dim; i++){
    for (j= 1; j<= grid_ptr->s_dim; j++){
      gp.r = (i-range(1,1)) * grid_ptr->r_step;
      gp.s = (j-range(1,2)) * grid_ptr->s_step;
      ierr = forward_mapping(&gp, map_ptr); 
      x(i,j)  = gp.x;
      y(i,j)  = gp.y;
      if (gp.x < grid_ptr->x_min) grid_ptr->x_min = gp.x;
      if (gp.x > grid_ptr->x_max) grid_ptr->x_max = gp.x;
      if (gp.y < grid_ptr->y_min) grid_ptr->y_min = gp.y;
      if (gp.y > grid_ptr->y_max) grid_ptr->y_max = gp.y;
    }
  }
/* error check */
  if (ierr != 0)
    printf("define_grid: error in grid for kg = %d\n",grid_ptr->priority);

/* update min and max for x and y in the overlapping grid */
  if (grid_ptr->x_min < over_ptr->x_min) over_ptr->x_min = grid_ptr->x_min;
  if (grid_ptr->x_max > over_ptr->x_max) over_ptr->x_max = grid_ptr->x_max;
  if (grid_ptr->y_min < over_ptr->y_min) over_ptr->y_min = grid_ptr->y_min;
  if (grid_ptr->y_max > over_ptr->y_max) over_ptr->y_max = grid_ptr->y_max;

/* ensure perfect periodicity if r_period = 1 or s_period = 1 */
  if (grid_ptr->r_period ==1)
    for (j=1; j<= grid_ptr->s_dim; j++){
      x(range(2,1),j) = x(range(1,1),j);
      y(range(2,1),j) = y(range(1,1),j);
    }

  if (grid_ptr->s_period ==1)
    for (i=1; i<= grid_ptr->r_dim; i++){
      x(i,range(2,2)) = x(i,range(1,2));
      y(i,range(2,2)) = y(i,range(1,2));
    }

/* ensure perfect periodicity for c-grids */
  if (grid_ptr->topology == 1)
    {
      if (grid_ptr->branch_direction == 1)
	{
/* branch cut */
	  for (i=grid_ptr->last_point; i<=grid_ptr->r_dim; i++)
	    {
	      x(i,range(1,2)) = x(grid_ptr->r_dim+1-i, range(1,2));
	      y(i,range(1,2)) = y(grid_ptr->r_dim+1-i, range(1,2));
	    }
/* ghost points */
	  for (j=1; j<=over_ptr->extra; j++)
	    {
	      for (i=1; i<=grid_ptr->first_point; i++)
		{
		  x(i,range(1,2)-j) = x(grid_ptr->r_dim+1-i, range(1,2)+j);
		  y(i,range(1,2)-j) = y(grid_ptr->r_dim+1-i, range(1,2)+j);
		}
	      for (i=grid_ptr->last_point; i<=grid_ptr->r_dim; i++)
		{
		  x(i,range(1,2)-j) = x(grid_ptr->r_dim+1-i, range(1,2)+j);
		  y(i,range(1,2)-j) = y(grid_ptr->r_dim+1-i, range(1,2)+j);
		}
	    } /* end for all ghost points */
	}
    }

} /* end define grid */


