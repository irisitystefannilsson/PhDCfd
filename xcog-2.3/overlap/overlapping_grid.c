#include "overlap.h"

/* static prototypes */
static void 
delete_grid( component_grid *grid_ptr );
static interp_point *
delete_interpolation_list( component_grid *grid_ptr );
static bad_point *
delete_bad_list( component_grid *grid_ptr );
static boundary_info *
new_boundary_info( int i_const, int j_const, int start, int stop, 
		  component_grid *grid_ptr );
static void 
delete_boundary_info( boundary_info *info );
/* static prototypes */


overlapping_grid *
new_overlapping_grid( int n_grids ){
  overlapping_grid *over_ptr;

  over_ptr = (overlapping_grid *) malloc(sizeof(overlapping_grid));
  
  over_ptr->n_grids = n_grids;

/* pointers to the first and last component grids */
  over_ptr->first   = NULL;
  over_ptr->last    = NULL;

/* number of ghost points */
  over_ptr->extra = 0;

/* number of extra points at periodic boundaries */
  over_ptr->extra_period = 0;

/* discretization width quantities */
  over_ptr->disc_width = 0;
  over_ptr->normal_width = 0;
  over_ptr->tangent_width = 0;
  over_ptr->corner_width = 0;
    
/* interpolation width */
  over_ptr->interp_width = 0;
    
/* interpolation type */
  over_ptr->interp_type = 'e';
    
/* default trimming style is to minimize the overlap */
  over_ptr->trim_style = 0;

/* correct for boundary mismatch */
  over_ptr->correct_mismatch = TRUE;

/* new head for the list of boundary info */
  over_ptr->boundary_list = new_linked_list();

/* bounding box */
  over_ptr->x_min =  1.e10;
  over_ptr->x_max = -1.e10;
  over_ptr->y_min =  1.e10;
  over_ptr->y_max = -1.e10;

/* tolerance for two curves to be regarded as the same curve */
  over_ptr->max_n_dist = 0.0;

/* title */
  over_ptr->title = NULL;

  return over_ptr;
}

void 
delete_overlapping_grid( overlapping_grid *over_ptr ){
  component_grid *grid_ptr, *next_grid_ptr;
  linked_list_member * this_boundary;

  if (over_ptr == NULL) return;

/* delete all mappings */
  grid_ptr = over_ptr->first;
  while (grid_ptr != NULL){
    next_grid_ptr = grid_ptr->next;
    delete_grid( grid_ptr);
    grid_ptr = next_grid_ptr;
  }

/* free the global boundary info */
  for (this_boundary = over_ptr->boundary_list->first; this_boundary != NULL;
       this_boundary = this_boundary->next){
    delete_boundary_part( (boundary_part *) this_boundary->data );
  }
/* free the list itself */
  over_ptr->boundary_list = delete_linked_list( over_ptr->boundary_list );

/* free the title */
  if (over_ptr->title != NULL)
    free( over_ptr->title );

/* free the structure itself */
  free( over_ptr );
}

component_grid *
new_component_grid( char *name){
  component_grid *grid_ptr;
  int i, j;

  grid_ptr = (component_grid *) malloc(sizeof(component_grid));
  
/* initialize all fields */
  grid_ptr->prev = NULL;
  grid_ptr->next = NULL;
  grid_ptr->priority = 0;

  grid_ptr->grid_name = (char *) malloc( (strlen(name)+1)*sizeof(char) );
  grid_ptr->grid_name = strcpy( grid_ptr->grid_name, name );

/* save the pointers to the data and the functions for
   evaluating the mapping */
  grid_ptr->inverse_known = 0;
  grid_ptr->grid_mapping = NULL;
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->mapping_data_ptr = NULL;

/* set the grid type and the periodicity flag */
  grid_ptr->grid_type= 0; /* unknown grid type */

/* set the number of grid lines */
  grid_ptr->r_points = 0;
  grid_ptr->s_points = 0;

  grid_ptr->r_dim = 0;
  grid_ptr->s_dim = 0;
  grid_ptr->range_ptr     = create_int_array_2d(2,2);

  grid_ptr->r_step = 0.0;
  grid_ptr->s_step = 0.0;

  grid_ptr->r_period = 0;
  grid_ptr->s_period = 0;

  grid_ptr->bc_ptr        = create_int_array_2d(2,2);
  grid_ptr->curve_ptr     = create_int_array_2d(2,2);

  grid_ptr->left_interp_ptr = NULL;
  grid_ptr->right_interp_ptr = NULL;
  grid_ptr->lower_interp_ptr = NULL;
  grid_ptr->upper_interp_ptr = NULL;

  grid_ptr->flag_ptr = NULL;
  grid_ptr->x_ptr = NULL;
  grid_ptr->y_ptr = NULL;

  grid_ptr->n_interp = 0;
  grid_ptr->last_interp = NULL;

  grid_ptr->last_bad_point = NULL;

  grid_ptr->x_min = 0.0;
  grid_ptr->x_max = 1.0;
  grid_ptr->y_min = 0.0;
  grid_ptr->y_max = 1.0;

/* initialize to black */
  grid_ptr->grid_color.r = 0.0;
  grid_ptr->grid_color.g = 0.0;
  grid_ptr->grid_color.b = 0.0;

  grid_ptr->cos_theta = 1.0;
  grid_ptr->sin_theta = 0.0;
  grid_ptr->x_trans = 0.0;
  grid_ptr->y_trans = 0.0;
  grid_ptr->scaling = 1.0;

/* initialize the 2-d arrays */
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++){
      bc(i,j) = 0;
      range(i,j) = 0;
      curve(i,j) = 0;
    }

/* initialize the boundary information */
  grid_ptr->left  = NULL;
  grid_ptr->right = NULL;
  grid_ptr->lower = NULL;
  grid_ptr->upper = NULL;

/* external gaps */
  grid_ptr->gap_low_r[0] =  2.0;
  grid_ptr->gap_low_r[1] = -1.0;
  grid_ptr->gap_high_r[0] =  2.0;
  grid_ptr->gap_high_r[1] = -1.0;
  grid_ptr->gap_low_s[0] =  2.0;
  grid_ptr->gap_low_s[1] = -1.0;
  grid_ptr->gap_high_s[0] =  2.0;
  grid_ptr->gap_high_s[1] = -1.0;

/* initialize the topology info */
  grid_ptr->topology = 0;
  grid_ptr->branch_direction = 0;

/* initialize the c-grid info */
  grid_ptr->first_point = 0;
  grid_ptr->last_point  = 0;

  return grid_ptr;
}

static void 
delete_grid( component_grid *grid_ptr ){

  if (grid_ptr == NULL) return;

/* deallocate the name */
  free( grid_ptr->grid_name );

/* free the space of the mapping */
  grid_ptr->range_ptr = delete_int_array_2d( grid_ptr->range_ptr );
  grid_ptr->bc_ptr = delete_int_array_2d( grid_ptr->bc_ptr );
  grid_ptr->curve_ptr = delete_int_array_2d( grid_ptr->curve_ptr );

  grid_ptr->left_interp_ptr = delete_int_array_1d( grid_ptr->left_interp_ptr );
  grid_ptr->right_interp_ptr = delete_int_array_1d( grid_ptr->right_interp_ptr );
  grid_ptr->lower_interp_ptr = delete_int_array_1d( grid_ptr->lower_interp_ptr );
  grid_ptr->upper_interp_ptr = delete_int_array_1d( grid_ptr->upper_interp_ptr );

  grid_ptr->x_ptr = delete_real_array_2d( grid_ptr->x_ptr );
  grid_ptr->y_ptr = delete_real_array_2d( grid_ptr->y_ptr );
  grid_ptr->flag_ptr = delete_int_array_2d( grid_ptr->flag_ptr );

/* free the list of interpolation points */
  grid_ptr->last_interp = delete_interpolation_list( grid_ptr );

/* free the list of bad points */
  grid_ptr->last_bad_point = delete_bad_list( grid_ptr );

/* free the boundary information trees */
  delete_boundary_info( grid_ptr->left  );
  delete_boundary_info( grid_ptr->right );
  delete_boundary_info( grid_ptr->lower );
  delete_boundary_info( grid_ptr->upper );

/* free the structure itself */
  free( grid_ptr );
}

boundary_part *
new_boundary_part(int curve_label){
  boundary_part *boundary_ptr;

  boundary_ptr = (boundary_part *) malloc( sizeof(boundary_part) );
  boundary_ptr->curve_label = curve_label;
  boundary_ptr->n_points = 0;
  boundary_ptr->top = 0;
  boundary_ptr->inside = -1; /* undecided */
  boundary_ptr->coordinates = NULL;

  return boundary_ptr;
}

boundary_part *
delete_boundary_part( boundary_part *boundary_ptr ){
  if (boundary_ptr){
/* free the array of coordinates */
    delete_real_array_2d( boundary_ptr->coordinates );
/* free the structure itself */
    free( boundary_ptr );
  }
  return NULL;
}

void 
subdivide_boundaries( component_grid *grid_ptr ){

/* create the top-level of the information tree */
  grid_ptr->left  = new_boundary_info( range(1,1), 0, range(1,2), range(2,2), 
				      grid_ptr );
  grid_ptr->right = new_boundary_info( range(2,1), 0, range(1,2), range(2,2), 
				      grid_ptr );
  grid_ptr->lower = new_boundary_info( 0, range(1,2), range(1,1), range(2,1), 
				      grid_ptr );
  grid_ptr->upper = new_boundary_info( 0, range(2,2), range(1,1), range(2,1), 
				      grid_ptr );
}

static boundary_info *
new_boundary_info( int i_const, int j_const, int start, int stop, 
		  component_grid *grid_ptr ){
  boundary_info *info=NULL;
  int i, j, mid;

/* check some things */
  if ( !(i_const==0 || j_const==0) ){
    printf("ERROR in new_boundary_info. Either i_const == 0 or j_const == 0,\n"
	   "but i_const=%i and j_const=%i\n", info->i_const, info->j_const);
    exit(1);
  }

  info = (boundary_info *) malloc( sizeof(boundary_info) );
  
  info->i_const = i_const;
  info->j_const = j_const;
  info->start = start;
  info->stop  = stop;

/* initialize the sub-trees */
  if (stop - start < 4){
    info->part1 = NULL;
    info->part2 = NULL;

/* initialize the bounding box so it becomes easy to compute it */
    info->xmin = 1.0e10;
    info->xmax =-1.0e10;
    info->ymin = 1.0e10;
    info->ymax =-1.0e10;
    
    if ( info->i_const == 0 ){
      j = info->j_const;
      for (i = start; i<= stop; i++){
	if (x(i,j) < info->xmin) info->xmin = x(i,j);
	if (x(i,j) > info->xmax) info->xmax = x(i,j);
	if (y(i,j) < info->ymin) info->ymin = y(i,j);
	if (y(i,j) > info->ymax) info->ymax = y(i,j);
      }
    }
    else if ( info->j_const == 0 ){
      i = info->i_const;
      for (j = start; j<= stop; j++){
	if (x(i,j) < info->xmin) info->xmin = x(i,j);
	if (x(i,j) > info->xmax) info->xmax = x(i,j);
	if (y(i,j) < info->ymin) info->ymin = y(i,j);
	if (y(i,j) > info->ymax) info->ymax = y(i,j);
      }
    }
  }
  else{
    mid = (start + stop)/2;
    info->part1 = new_boundary_info( i_const, j_const, start, mid, grid_ptr );
    info->part2 = new_boundary_info( i_const, j_const, mid,  stop, grid_ptr );

/* take the bounding box from the sub-trees */
    info->xmin = real_min( info->part1->xmin, info->part2->xmin );
    info->xmax = real_max( info->part1->xmax, info->part2->xmax );
    info->ymin = real_min( info->part1->ymin, info->part2->ymin );
    info->ymax = real_max( info->part1->ymax, info->part2->ymax );
  }

   
/* done */
  return info;
}

static void 
delete_boundary_info( boundary_info *info ){

/* are we at the bottom of the tree yet? */
  if (info == NULL) return;

/* delete the information tree recursively */
  delete_boundary_info( info->part1 );
  delete_boundary_info( info->part2 );

/* free the structure itself */
  free( info );

/* done */
}

boundary_part *
find_curve( linked_list *boundary_list, int curve_label ){
  linked_list_member *this_boundary;
  boundary_part *boundary_ptr=NULL;
/* go through the list and check if we already have a node for the current curve_label */
  for (this_boundary = boundary_list->first; this_boundary != NULL;
       this_boundary = this_boundary->next){
    boundary_ptr = this_boundary->data;
    if (boundary_ptr->curve_label == curve_label) break;
  } /* end for all present parts of the boundary */

/* did we find a link with matching curve_label? */
  if (this_boundary)
    return boundary_ptr;
  else
    return NULL;
}

#define boundary_coord(i,j) compute_index_2d(boundary_ptr->coordinates, i, j)
void
plot_global_boundary( boundary_part *boundary_ptr, int plot_mode ){
  int i;
  color_type black={0.0, 0.0, 0.0};
  const int h_center = 1, v_center = -1, font_number = 3;
  char label[80];
/* make sure there is something to plot */
  if (boundary_ptr == NULL) return;
  if (boundary_ptr->top <= 0) return;

/* set color */
  PL_color( black );
  PL_poll();
/* draw a polygon through all node points */
  if (plot_mode & 1){
    PL_move( boundary_coord(1,1), boundary_coord(2,1) );
    for (i=2; i<= boundary_ptr->top; i++)
      PL_draw( boundary_coord(1,i), boundary_coord(2,i) );
/* draw the final segment */
    PL_draw( boundary_coord(1,1), boundary_coord(2,1) );

  }
/* plot labeled markers at all boundary points */
  if (plot_mode & 2){
    for (i=1; i<= boundary_ptr->top; i++){
      PL_marker( boundary_coord(1,i), boundary_coord(2,i), CIRCLE );
/*       sprintf(label, "%i", i); */
/*       PL_plot_string(label, boundary_coord(1,i), boundary_coord(2,i),  */
/* 		     h_center, v_center, font_number); */
    }
  }
/* plot the curve-label next to the curve */
  if (plot_mode & 4){
    i=1;
    sprintf(label, "curve-label=%i", boundary_ptr->curve_label); 
    PL_plot_string(label, boundary_coord(1,i), boundary_coord(2,i),  
		   h_center, v_center, font_number); 
  }
} /* end plot_global_boundary */
#undef boundary_coord


void 
plot_overlapping(overlapping_grid *over_ptr, int plot_mode ){
  int i,j,h_center,v_center,font_number, prev_real;
  char c_flag[20];
  component_grid *grid_ptr;
  bad_point *bad_ptr;
  interp_point *interp_ptr;
/*
   plot_mode controls what is plotted according to:
   plot_mode = 1  composite gridpoints.
   plot_mode = 2  CIRCLE all interpolation points.
   plot_mode = 4  mark the bad points according to their type.
   plot_mode = 8  flag array at each non-fictious point.
   plot_mode = 16 title
   plot_mode = 32 physical boundary points and the physical part of mixed boundaries

   Different options are combined by or-ing, i.e. plot_mode = 1 | 2 gives both
   option 1 and 2. 
*/

/* only plot the grid if there is one! */
  if (over_ptr == NULL) return;

  for (grid_ptr = over_ptr->first; grid_ptr != NULL;
       grid_ptr = grid_ptr->next){
    
/* set color */
    PL_color( grid_ptr->grid_color );
/* draw the non-fictitious points with flag !=0 in solid */
    if ((plot_mode & 1) && grid_ptr->flag_ptr != NULL && 
	grid_ptr->x_ptr != NULL && grid_ptr->y_ptr != NULL){
/*    PL_shade(1.0); */
/* iso i-lines */
      for (i = range(1,1); i <= range(2,1); i++){
	j = range(1,2);
	PL_poll();
	PL_move(x(i,j),y(i,j));
	prev_real = ( flag(i,j) != 0 );
	for (j = range(1,2)+1; j <= range(2,2); j++){
	  if (flag(i,j) != 0 && prev_real){
	    PL_draw(x(i,j),y(i,j));
	    prev_real = 1;
	  }
	  else if (flag(i,j) != 0){
	    prev_real = 1;
	    PL_move(x(i,j),y(i,j));
	  }
	  else
	    prev_real = 0; 
	}
      }
/* iso j-lines */
      for (j = range(1,2); j <= range(2,2); j++){
	i = range(1,1);
	PL_poll();
	PL_move(x(i,j),y(i,j));
	prev_real = ( flag(i,j) != 0 );
	for (i = range(1,1)+1; i<= range(2,1); i++){
	  if (flag(i,j) != 0 && prev_real){
	    PL_draw(x(i,j),y(i,j));
	    prev_real = 1;
	  }
	  else if (flag(i,j) != 0){
	    prev_real = 1;
	    PL_move(x(i,j),y(i,j));
	  }
	  else
	    prev_real = 0; 
	}
      }
    }
    
    if ( (plot_mode & 2) && grid_ptr->x_ptr != NULL &&
	grid_ptr->y_ptr != NULL){
      h_center = 1; v_center = 1; font_number = 3;
/* mark all interpolation points with an CIRCLE */
      PL_poll();
      for (interp_ptr = grid_ptr->last_interp;
	   interp_ptr != NULL; interp_ptr = interp_ptr->prev){
/* don't plot inactive interpolation points */
	if (interp_ptr->active){
/* set the color to equal that of the donor grid */
	  PL_color( interp_ptr->grid_loc->grid_color );
	  i = interp_ptr->i_point; j = interp_ptr->j_point;
	  PL_marker(x(i,j), y(i,j), CIRCLE);
	}
      }
    }
    
/* mark all bad points */
    if ( (plot_mode & 4) ){
      PL_poll();
      for (bad_ptr = grid_ptr->last_bad_point; bad_ptr != NULL;
	   bad_ptr = bad_ptr->prev){
	PL_marker(x(bad_ptr->i, bad_ptr->j), 
		  y(bad_ptr->i, bad_ptr->j), 
		  bad_ptr->type);
      }
    }
    
/* plot the value of grid_ptr->flag(i,j) on the location of each grid point */
    if ((plot_mode & 8) && grid_ptr->flag_ptr != NULL && 
	grid_ptr->x_ptr != NULL && grid_ptr->y_ptr != NULL){
      h_center = 1; v_center = 1; font_number = 3;
      for (i = range(1,1); i <= range(2,1); i++){
	PL_poll();
	for (j = range(1,2); j <= range(2,2); j++){
/* only plot the used, non-discretization points */
	  if (flag(i,j) != 0 && flag(i,j) != grid_ptr->priority){
	    sprintf(c_flag, "%2i", flag(i,j));
	    PL_plot_string(c_flag, x(i,j), y(i,j), h_center, v_center, font_number);
	  }
	}
      }
    }
    
/* physical boundaries */
    if ( plot_mode & 32 ){
      PL_poll();
/* left */
      if (abs(curve(1,1)) > 0){
	j = range(1,2);
	i = range(1,1);
	PL_move(x(i,j),y(i,j));
	prev_real = ( left_interp(j) > 0 );
	for (j = range(1,2)+1; j <= range(2,2); j++){
	  if (prev_real && left_interp(j) > 0 ){
	    PL_draw(x(i,j),y(i,j));
	    prev_real = 1;
	  }
	  else if (left_interp(j) > 0){
	    prev_real = 1;
	    PL_move(x(i,j),y(i,j));
	  }
	  else
	    prev_real = 0; 
	} /* end for j */
      } /* end physical or mixed boundary */

/* right */
      if (abs(curve(2,1)) > 0){
	j = range(1,2);
	i = range(2,1);
	PL_move(x(i,j),y(i,j));
	prev_real = ( right_interp(j) > 0 );
	for (j = range(1,2)+1; j <= range(2,2); j++){
	  if (prev_real && right_interp(j) > 0){
	    PL_draw(x(i,j),y(i,j));
	    prev_real = 1;
	  }
	  else if ( right_interp(j) > 0 ){
	    prev_real = 1;
	    PL_move(x(i,j),y(i,j));
	  }
	  else
	    prev_real = 0; 
	} /* end for j */
      } /* end physical or mixed boundary */

/* lower */
      if (abs(curve(1,2)) > 0){
	j = range(1,2);
	i = range(1,1);
	PL_move(x(i,j),y(i,j));
	prev_real = ( lower_interp(i) > 0 );
	for (i = range(1,1)+1; i <= range(2,1); i++){
	  if (prev_real && lower_interp(i) > 0){
	    PL_draw(x(i,j),y(i,j));
	    prev_real = 1;
	  }
	  else if (lower_interp(i) > 0){
	    prev_real = 1;
	    PL_move(x(i,j),y(i,j));
	  }
	  else
	    prev_real = 0; 
	} /* end for i */
      } /* end physical or mixed boundary */

/* upper */
      if (abs(curve(2,2)) > 0){
	j = range(2,2);
	i = range(1,1);
	PL_move(x(i,j),y(i,j));
	prev_real = ( upper_interp(i) > 0 );
	for (i = range(1,1)+1; i <= range(2,1); i++){
	  if ( prev_real && upper_interp(i) > 0 ){
	    PL_draw(x(i,j),y(i,j));
	    prev_real = 1;
	  }
	  else if ( upper_interp(i) > 0 ){
	    prev_real = 1;
	    PL_move(x(i,j),y(i,j));
	  }
	  else
	    prev_real = 0; 
	} /* end for i */
      } /* end physical or mixed boundary */
    } /* end plot_mode & 32 */

  } /* end for all grids */

/* title */
  if ( plot_mode & 16 ){
    PL_title( 0, over_ptr->title );
  }

} /* end grid_plot */




static interp_point *
delete_interpolation_list( component_grid *grid_ptr ){
  interp_point *interp_ptr, *next_victim;

  interp_ptr = grid_ptr->last_interp;
  while (interp_ptr != NULL){
    next_victim = interp_ptr->prev;
    free(interp_ptr);
    interp_ptr= next_victim;
  }
  grid_ptr->n_interp = 0;

  return NULL;
}

static bad_point *
delete_bad_list( component_grid *grid_ptr ){
  bad_point *bad_point_ptr, *next_victim;

  bad_point_ptr = grid_ptr->last_bad_point;
  while (bad_point_ptr != NULL){
    next_victim = bad_point_ptr->prev;
    free(bad_point_ptr);
    bad_point_ptr = next_victim;
  }
  
  return NULL;
}


void 
insert_grid( component_grid *grid_ptr, overlapping_grid *over_ptr ){
/* insert the component grid pointed to by `grid_ptr' into the linked list */
  if (over_ptr->first != NULL)
    over_ptr->first->prev = grid_ptr;
  grid_ptr->next = over_ptr->first;
  over_ptr->first = grid_ptr;
  if (over_ptr->last == NULL)
    over_ptr->last = grid_ptr;
}


int 
forward_grid_mapping( grid_point *gp_ptr, component_grid *grid_ptr ){
  grid_point gp_old;
  int msg;

  gp_old.r = gp_ptr->r;
  gp_old.s = gp_ptr->s;

/* take care of any errors in the specific mapping routine and pass them on */
  if ( (msg=(*grid_ptr->grid_mapping)( &gp_old, grid_ptr->mapping_data_ptr )) != 0)
    return msg;

/* the transformation is done in the order: rotation - translation - scaling */
  gp_ptr->x = grid_ptr->scaling* (grid_ptr->x_trans + 
				 gp_old.x * grid_ptr->cos_theta -
				 gp_old.y * grid_ptr->sin_theta);

  gp_ptr->y = grid_ptr->scaling* (grid_ptr->y_trans + 
				 gp_old.x * grid_ptr->sin_theta +
				 gp_old.y * grid_ptr->cos_theta);

/* jacobian */
  gp_ptr->xr = grid_ptr->scaling* (gp_old.xr * grid_ptr->cos_theta -
				  gp_old.yr * grid_ptr->sin_theta);

  gp_ptr->xs = grid_ptr->scaling* (gp_old.xs * grid_ptr->cos_theta -
				  gp_old.ys * grid_ptr->sin_theta);
				 
  gp_ptr->yr = grid_ptr->scaling* (gp_old.xr * grid_ptr->sin_theta +
				  gp_old.yr * grid_ptr->cos_theta);

  gp_ptr->ys = grid_ptr->scaling* (gp_old.xs * grid_ptr->sin_theta +
				  gp_old.ys * grid_ptr->cos_theta);
				 
  return 0;
}


int 
backward_grid_mapping( grid_point *gp_ptr, component_grid *grid_ptr ){
  int status;
  real x_new, y_new;

  x_new = gp_ptr->x;
  y_new = gp_ptr->y;

/* invert the rotation - translation - scaling */
  gp_ptr->x = (x_new/grid_ptr->scaling - grid_ptr->x_trans) * grid_ptr->cos_theta +
    (y_new/grid_ptr->scaling - grid_ptr->y_trans) * grid_ptr->sin_theta;

  gp_ptr->y =-(x_new/grid_ptr->scaling - grid_ptr->x_trans) * grid_ptr->sin_theta +
    (y_new/grid_ptr->scaling - grid_ptr->y_trans) * grid_ptr->cos_theta;

/* invert the un-transformed (x,y) coordinates */
  status = (*grid_ptr->inverse_grid_mapping)( gp_ptr, grid_ptr->mapping_data_ptr );

/* restore gp_ptr->x, gp_ptr->y */
  gp_ptr->x = x_new;
  gp_ptr->y = y_new;

  return status;
}

