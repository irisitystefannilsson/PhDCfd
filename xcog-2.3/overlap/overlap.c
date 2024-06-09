#include "overlap.h"

/* static prototypes */
static int
close_to_boundary(real xp, real yp, int curve_label, overlapping_grid *over_ptr);
static void
narrow_interp_width(overlapping_grid *over_ptr);
static void
bridge_the_gap(int i0, int j0, real x0, real y0,
	       int i1, int j1, real x1, real y1, 
	       component_grid *hole_ptr, int cutter_curve, overlapping_grid *over_ptr);
static int
mark_one_cell(real x_c, real y_c, int cutter_curve, component_grid *hole_ptr, 
	      overlapping_grid *comp_ptr, int *i_cell, int *j_cell);
static int
inside_global_domain( real xp, real yp, int curve_label, overlapping_grid *over_ptr );
static int
inside_cutter_curve(real xp, real yp, boundary_part *boundary_ptr);
static int 
step_1(overlapping_grid *comp_ptr);
static int 
step_2(overlapping_grid *comp_ptr);
static int 
step_3(input_output *io_ptr, overlapping_grid *comp_ptr);
static int 
step_4(overlapping_grid *comp_ptr);
static int 
step_5(overlapping_grid *comp_ptr);
static int 
step_6(overlapping_grid *comp_ptr);
static int 
step_7(overlapping_grid *comp_ptr);
static int
normal_distance(grid_point *gp_ptr, int cutter_curve, component_grid *grid_ptr, 
		overlapping_grid *comp_ptr);
static int 
needed_by_disc(int i_point, int j_point, component_grid *grid_ptr,
	       overlapping_grid *comp_ptr);
static int 
needed_by_interp(int i_point, int j_point, int priority,
		 overlapping_grid *comp_ptr);
static void 
remove_deactive_interp(component_grid *grid_ptr);
static void 
set_flag(component_grid *grid_ptr, int flag_value, int extra_period);
static int 
explicit(input_output *io_ptr, overlapping_grid *comp_ptr);
static int 
disc_point(int i, int j, component_grid *grid_ptr, 
	   overlapping_grid *comp_ptr);
static int 
interp_from_higher(real x, real y, int curve, 
		   interp_point *new_interp_ptr, 
		   component_grid *other_grid_ptr, overlapping_grid *comp_ptr);
static int 
interp_from_lower(real x, real y, int curve, 
		  interp_point *new_interp_ptr,
		  component_grid *own_grid_ptr, 
		  component_grid *other_grid_ptr, overlapping_grid *comp_ptr);
static int 
good_interp_loc(interp_point *new_interp_ptr,
		component_grid *grid_ptr,
		overlapping_grid *comp_ptr);
static interp_point *
create_interp_point(interp_point *template_ptr, 
		    interp_point *prev_ptr, component_grid *grid_ptr);
static void 
change_sign(int i_loc, int j_loc, int iw1, component_grid *grid_ptr);
static int 
invert_mapping(grid_point *gp_ptr, int i_cell, int j_cell, int curve, 
	       component_grid *grid_ptr, overlapping_grid *comp_ptr);
static int 
get_curve(int i, int j, component_grid *grid_ptr, 
	  overlapping_grid *comp_ptr);
static void 
insert_bad_point(int i, int j, int type, component_grid *grid_ptr);
static int 
newton(component_grid *grid_ptr, grid_point *goal_ptr, real r0, real s0);
static void 
closest_boundary_point( component_grid *grid_ptr, real x0, real y0, 
		       int job, int *ip, int *jp, 
		       real *n_dist, real *parameter );
static int 
inside_region(real xp, real yp, component_grid *grid_ptr);
static int 
inside_sub_region(real xp, real yp, int i_min, int i_max, int j_min, 
		  int j_max, component_grid *grid_ptr);
static int 
locate_cell(real xp, real yp, int *i_coord, int *j_coord, 
	    int i_cell, int j_cell, component_grid *grid_ptr);
static int 
boundary_interp( real xp, real yp, component_grid *this_grid, 
		overlapping_grid *comp_ptr);
static int 
local_newton(component_grid *grid_ptr, grid_point *goal_ptr,
	     int interp_width, int i_close, int j_close,
	     real r0, real s0);
static int 
ray_intersect( real xp, real yp, boundary_info *info, 
	      component_grid *grid_ptr );
/* end static prototypes */



#include <sys/types.h>
#include <sys/times.h>

int 
overlap(input_output *io_ptr, overlapping_grid *comp_ptr, int redo, 
	int show_phys, int show_holes){
  linked_list_member *this_boundary;
  boundary_part *boundary_ptr;
  int step_3_ok, icom, col_flag, actual_iw=0;
  struct tms time0, time1, time2, time3, time4, time5, time6;
  char * break_com[3], *file_name;
  const int LEVEL=1;
  FILE *fp;
  break_com[0] = "proceed";
  break_com[1] = "break";
  break_com[2] = "postscript-copy";

  times(&time0);

  if (comp_ptr->interp_width == 2 && comp_ptr->interp_type == 'i'){
    actual_iw = 2;
    comp_ptr->interp_width = 3;
  }

  if (!redo) step_1(comp_ptr);
  times(&time1);


  if (show_phys){
      PL_erase(0);
      PL_start_plot(0);
      for (this_boundary = comp_ptr->boundary_list->first; this_boundary != NULL;
	   this_boundary = this_boundary->next){
	boundary_ptr = this_boundary->data;
/* draw a polygon between the boundary grid points and print the corresponding */
/* curve_label */
	plot_global_boundary( boundary_ptr, 1|4); 
      } /* end for each part of the boundary */
      PL_end_plot();
/* break or proceed? */    
      while ((icom=get_command(io_ptr, "You are looking at the physical boundary>", 
			       break_com, 2, LEVEL, NULL, NULL)) != 0)
	{
	  switch (icom)
	    {
	    case 1:		/* exit to calling routine */
	      return FALSE;
	      break;
	    case 2:		/* postscript copy */
	      if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
					 "boundary.ps", &file_name, 'w', 0, 
					 TRUE) ) != NULL)
		{
		  col_flag = get_yes_no( io_ptr, "Do you want to save color "
					"information>", LEVEL );
		  PL_postscript( 0, fp, col_flag );
		  fclose( fp );
		}
	      break;
	    default:
	      ;
	    } /* end switch */
	} /* end while */
  } /* end if show_phys */

  if (!redo) step_2(comp_ptr);
  times(&time2);

/* plot the physical and mixed boundaries */
#ifdef PLOT_MIXED
  PL_erase(0);
  PL_start_plot(0);
  plot_overlapping(comp_ptr, 32);
  PL_end_plot();
/* break or proceed? */    
  while ((icom=get_command(io_ptr, "Physical boundaries for the purpose "
			   "of interpolation> ", 
			   break_com, 2, LEVEL, NULL, NULL)) != 0)
    {
      switch (icom)
	{
	case 1:			/* exit to calling routine */
	  return FALSE;
	  break;
	case 2:			/* postscript copy */
	  if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
				     "hole.ps", &file_name, 'w', 0, 
				     FALSE) ) != NULL)
	    {
	      col_flag = get_yes_no( io_ptr, "Do you want to save color "
				    "information>", LEVEL );
	      PL_postscript( 0, fp, col_flag );
	      fclose( fp );
	    }
	  break;
	default:
	  ;
	} /* end switch */
    } /* end while */
#endif

/* plot the overlapping grid */
  if (show_holes){
    PL_erase(0);
    PL_start_plot(0);
    plot_overlapping( comp_ptr, 1 ); /* used points */
    PL_end_plot();

/* break or proceed? */    
    while ((icom=get_command(io_ptr, "Used points after the holes have been cut>", 
			     break_com, 2, LEVEL, NULL, NULL)) != 0)
      {
	switch (icom)
	  {
	  case 1:		/* exit to calling routine */
	    return FALSE;
	    break;
	  case 2:		/* postscript copy */
	    if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
				       "hole.ps", &file_name, 'w', 0, 
				       FALSE) ) != NULL)
	      {
		col_flag = get_yes_no( io_ptr, "Do you want to save color "
				      "information>", LEVEL );
		PL_postscript( 0, fp, col_flag );
		fclose( fp );
	      }
	    break;
	  default:
	    ;
	  } /* end switch */
      } /* end while */
  } /* end if show_holes */

  if ( !(step_3_ok = step_3(io_ptr, comp_ptr)) ) 
    return step_3_ok;
  times(&time3);

#ifdef STOP_AFTER_STEP_3
  PL_erase(0);
  PL_start_plot(0);
  plot_overlapping( comp_ptr, 1+2 ); /* used grid points + interpolation points */
  PL_end_plot();
/* break or proceed? */    
  while ((icom=get_command(io_ptr, "After a successful main classification>", 
			   break_com, 2, LEVEL, NULL, NULL)) != 0)
    {
      switch (icom)
	{
	case 1:			/* exit to calling routine */
	  return FALSE;
	  break;
	case 2:			/* postscript copy */
	  if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
				     "classify.ps", &file_name, 'w', 0, 
				     FALSE) ) != NULL)
	    {
	      col_flag = get_yes_no( io_ptr, "Do you want to save color "
				    "information>", LEVEL );
	      PL_postscript( 0, fp, col_flag );
	      fclose( fp );
	    }
	  break;
	default:
	  ;
	} /* end switch */
    } /* end while */
#endif

  step_4(comp_ptr);
  times(&time4);

  step_5(comp_ptr);
  times(&time5);

  if (actual_iw){
    printf("Recomputing interpolation locations for interpolation width = 2\n"
	   "and implicit interpolation.\n");
    comp_ptr->interp_width = actual_iw;
    narrow_interp_width(comp_ptr);
  }

  step_6(comp_ptr);
  times(&time6);

  if (redo) step_7(comp_ptr);

  printf("\nTiming information for the overlap algorithm:\n");
  printf("Initialization:                          %.1f seconds\n", 
	 (time1.tms_utime-time0.tms_utime)/60.0);
  printf("Preparing physical and mixed boundaries: %.1f seconds\n", 
	 (time2.tms_utime-time1.tms_utime)/60.0);
  printf("Initial classification and explicit:     %.1f seconds\n", 
	 (time3.tms_utime-time2.tms_utime)/60.0);
  printf("Preparation for trimming:                %.1f seconds\n", 
	 (time4.tms_utime-time3.tms_utime)/60.0);
  printf("Trimming:                                %.1f seconds\n", 
	 (time5.tms_utime-time4.tms_utime)/60.0);
  printf("Changing sign on interpolation points:   %.1f seconds\n", 
	 (time6.tms_utime-time5.tms_utime)/60.0);
  printf("\nTotal time:                              %.1f seconds\n", 
	 (time6.tms_utime-time0.tms_utime)/60.0);
  
  return 1;
}

static void
narrow_interp_width(overlapping_grid *over_ptr){
  component_grid *grid1_, *grid_ptr;
  interp_point *interp_ptr;
  int i, j;
  
/* this function is only called with interpolation width = 2 */

  for (grid1_ = over_ptr->first; 
       grid1_ != NULL;
       grid1_ = grid1_->next){
    for (interp_ptr = grid1_->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev){

/* make the range macro work */
      grid_ptr = interp_ptr->grid_loc;

/* get the enclosing cell */
      i = range(1,1) + interp_ptr->r_loc * (range(2,1) - range(1,1));
      j = range(1,2) + interp_ptr->s_loc * (range(2,2) - range(1,2));

/* make sure the closest cell is inside the grid */
      interp_ptr->i_loc = int_min( int_max( i, range(1,1) ), range(2,1)-1 );
      interp_ptr->j_loc = int_min( int_max( j, range(1,2) ), range(2,2)-1 );

    } /* end for all interpolation points */
  } /* end for all grids */
}

static void
add_curve( linked_list *boundary_list, int curve_label, int n_points ){
  linked_list_member *this_boundary;
  boundary_part *boundary_ptr=NULL;

/* is there anything to add? */
  if (n_points <= 0) return;
/* go through the list and check if we already have a node for the current curve_label */
  for (this_boundary = boundary_list->first; this_boundary != NULL;
       this_boundary = this_boundary->next){
    boundary_ptr = this_boundary->data;
    if (boundary_ptr->curve_label == curve_label) break;
  } /* end for all present parts of the boundary */

/* did a link already exist? */
  if (this_boundary){
    boundary_ptr->n_points += n_points;
  }
  else{
    boundary_ptr = new_boundary_part( curve_label );
    boundary_ptr->n_points = n_points;
    new_link(boundary_list)->data = boundary_ptr;
  }
}

#define boundary_coord(i,j) compute_index_2d(boundary_ptr->coordinates, i, j)
static real
dist2(int i, int j, real t_x, real t_y, boundary_part *boundary_ptr ){
  real t_dist, n_dist, d;
#define sqr(x) ((x)*(x))
  t_dist = t_x * (boundary_coord(1,j) - boundary_coord(1,i)) +
    t_y * (boundary_coord(2,j) - boundary_coord(2,i));
  n_dist = -t_y * (boundary_coord(1,j) - boundary_coord(1,i)) +
    t_x * (boundary_coord(2,j) - boundary_coord(2,i));
/* Elliptic weighted distance. */
  if (t_dist >= 0.0){
/* Normal semi-axis 1, tangential semi-axis 5. */
    d = sqr( t_dist/5.0 ) + sqr( n_dist ) ;
  }
  else{
/* Normal semi-axis 1, tangential semi-axis 1/5.0. */
    d = sqr( t_dist*5.0 ) + sqr( n_dist ) ;
  }
  return d;
}

static void
sort_coordinates( boundary_part *boundary_ptr ){
  int i_p, j, i_swap;
  real x_tmp, y_tmp, min_dist, t_x, t_y;
/* compute the first local tangent and normal */
  i_p = 1;
  t_x = boundary_coord(1,i_p+1) - boundary_coord(1,i_p);
  t_y = boundary_coord(2,i_p+1) - boundary_coord(2,i_p);
  for (i_p=1; i_p<boundary_ptr->n_points; i_p++){
/* look for the closest point in the rest of the array */
    min_dist = dist2(i_p, i_p+1, t_x, t_y, boundary_ptr);
    i_swap = i_p+1;
    for (j=i_p+2; j<=boundary_ptr->n_points; j++){
      if (dist2(i_p, j, t_x, t_y, boundary_ptr) < min_dist){
	min_dist = dist2(i_p, j, t_x, t_y, boundary_ptr);
	i_swap = j;
      }
    } /* end for j */
/* swap i_p+1 and i_swap */
    if (i_swap != i_p+1){
      x_tmp = boundary_coord(1,i_p+1);
      y_tmp = boundary_coord(2,i_p+1);
      boundary_coord(1,i_p+1) = boundary_coord(1,i_swap);
      boundary_coord(2,i_p+1) = boundary_coord(2,i_swap);
      boundary_coord(1,i_swap) = x_tmp;
      boundary_coord(2,i_swap) = y_tmp;
    }
/* only compute a new tangent if the distance was finite */
    if (min_dist > 0.0){
/* compute next local tangent */
      t_x = boundary_coord(1,i_p+1) - boundary_coord(1,i_p);
      t_y = boundary_coord(2,i_p+1) - boundary_coord(2,i_p);
    }
  } /* end for i_p */
#undef sqr
}

static int 
inside_other_grid( real xp, real yp, component_grid *this_grid, 
		  overlapping_grid *over_ptr){
  component_grid *grid_ptr;
  grid_point gp;
  
  for (grid_ptr = over_ptr->first; 
       grid_ptr != NULL;
       grid_ptr = grid_ptr->next)
    if (grid_ptr != this_grid){
/* call the grid function to invert the mapping */
      gp.x = xp; gp.y = yp;
      if (invert_mapping(&gp, 0, 0, 0, grid_ptr, over_ptr)) return TRUE;
    }
/* no match */
  return FALSE;
}

/* STEP 1: Initialize flag and compute the global boundary */
static int 
step_1(overlapping_grid *over_ptr){
  component_grid *grid_ptr;
  linked_list_member *this_boundary;
  boundary_part *boundary_ptr;
  real r, s, x_test=0.0, y_test=0.0, n_x, n_y, length, max_eps, this_eps;
  int i, j, phys_points, i_p;

  printf("\n");
  for (grid_ptr = over_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){
    printf("Initializing the flag array and setting up the search tree in grid `%s'\n",
	   grid_ptr->grid_name);
    set_flag( grid_ptr, grid_ptr->priority, over_ptr->extra_period );
    subdivide_boundaries( grid_ptr );
  }

  printf("\n");
  printf("Scanning all components for physical boundaries\n");
  for (grid_ptr = over_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){

    if ( curve(1,1) > 0 ){
/* check for gaps */
      if (grid_ptr->gap_low_r[0] > 1.0 || grid_ptr->gap_low_r[1] < 0.0 ||
	  grid_ptr->gap_low_r[1] < grid_ptr->gap_low_r[0])
	phys_points = grid_ptr->s_points;
      else{
	phys_points = 0;
	for (j=range(1,2); j<=range(2,2); j++){
	  s = (j-range(1,2)) * grid_ptr->s_step;
	  if (!(grid_ptr->gap_low_r[0] <= s && s <= grid_ptr->gap_low_r[1])) 
	    phys_points++;
	} /* end for i */
      }
      add_curve( over_ptr->boundary_list, curve(1,1), phys_points );
    } /* end if curve(1,1) > 0 */

    if( curve(2,1) > 0 ){
/* check for gaps */
      if (grid_ptr->gap_high_r[0] > 1.0 || grid_ptr->gap_high_r[1] < 0.0 ||
	  grid_ptr->gap_high_r[1] < grid_ptr->gap_high_r[0])
	phys_points = grid_ptr->s_points;
      else{
	phys_points = 0;
	for (j=range(1,2); j<=range(2,2); j++){
	  s = (j-range(1,2)) * grid_ptr->s_step;
	  if (!(grid_ptr->gap_high_r[0] <= s && s <= grid_ptr->gap_high_r[1])) 
	    phys_points++;
	} /* end for i */
      }
      add_curve( over_ptr->boundary_list, curve(2,1), phys_points );
    } /* end if curve(2,1) > 0 */

    if( curve(1,2) > 0 ){
/* check for c-topology */
      if (grid_ptr->topology == 1){
	phys_points = grid_ptr->last_point - grid_ptr->first_point + 1;
	printf("Detected c-topology with %i physical points in grid %s\n", 
	       phys_points, grid_ptr->grid_name);
      }
/* check for gaps */
      else if (grid_ptr->gap_low_s[0] > 1.0 || grid_ptr->gap_low_s[1] < 0.0 ||
	  grid_ptr->gap_low_s[1] < grid_ptr->gap_low_s[0])
	phys_points = grid_ptr->r_points;
      else{
	phys_points = 0;
	for (i=range(1,1); i<=range(2,1); i++){
	  r = (i-range(1,1)) * grid_ptr->r_step;
	  if (!(grid_ptr->gap_low_s[0] <= r && r <= grid_ptr->gap_low_s[1])) 
	    phys_points++;
	} /* end for i */
      }
      add_curve( over_ptr->boundary_list, curve(1,2), phys_points );
    } /* end if curve(2,1) > 0 */

    if( curve(2,2) > 0 ){
/* check for gaps */
      if (grid_ptr->gap_high_s[0] > 1.0 || grid_ptr->gap_high_s[1] < 0.0 ||
	  grid_ptr->gap_high_s[1] < grid_ptr->gap_high_s[0])
	phys_points = grid_ptr->r_points;
      else{
	phys_points = 0;
	for (i=range(1,1); i<=range(2,1); i++){
	  r = (i-range(1,1)) * grid_ptr->r_step;
	  if (!(grid_ptr->gap_high_s[0] <= r && r <= grid_ptr->gap_high_s[1])) 
	    phys_points++;
	} /* end for i */
      }
      add_curve( over_ptr->boundary_list, curve(2,2), phys_points );
    }

/* mixed physical / interpolating sides */

    if ( curve(1,1) < 0 ){
/* Only the boundary points that are outside of all other components are */
/* considered to be physical boundary points */
      phys_points = 0;
      i = range(1,1);
      for (j=range(1,2); j<=range(2,2); j++){
	if (!inside_other_grid( x(i,j), y(i,j), grid_ptr, over_ptr )){
	  phys_points++;
	  left_interp(j) = 1;
	} /* end if */
      } /* end for all boundary points */
      add_curve( over_ptr->boundary_list, abs(curve(1,1)), phys_points );
    } /* end if curve(1,1) < 0 */

    if ( curve(2,1) < 0 ){
/* Only the boundary points that are outside of all other components are */
/* considered to be physical boundary points */
      phys_points = 0;
      i = range(2,1);
      for (j=range(1,2); j<=range(2,2); j++){
	if (!inside_other_grid( x(i,j), y(i,j), grid_ptr, over_ptr )){
	  phys_points++;
	  right_interp(j) = 1;
	} /* end if */
      } /* end for all boundary points */
      add_curve( over_ptr->boundary_list, abs(curve(2,1)), phys_points );
    } /* end if curve(2,1) < 0 */

    if ( curve(1,2) < 0 ){
/* Only the boundary points that are outside of all other components are */
/* considered to be physical boundary points */
      phys_points = 0;
      j = range(1,2);
      for (i=range(1,1); i<=range(2,1); i++){
	if (!inside_other_grid( x(i,j), y(i,j), grid_ptr, over_ptr )){
	  phys_points++;
	  lower_interp(i) = 1;
	} /* end if */
      } /* end for all boundary points */
      add_curve( over_ptr->boundary_list, abs(curve(1,2)), phys_points );
    } /* end if curve(1,2) < 0 */

    if ( curve(2,2) < 0 ){
/* Only the boundary points that are outside of all other components are */
/* considered to be physical boundary points */
      phys_points = 0;
      j = range(2,2);
      for (i=range(1,1); i<=range(2,1); i++){
	if (!inside_other_grid( x(i,j), y(i,j), grid_ptr, over_ptr )){
	  phys_points++;
	  upper_interp(i) = 1;
	} /* end if */
      } /* end for all boundary points */
      add_curve( over_ptr->boundary_list, abs(curve(2,2)), phys_points );
    } /* end if curve(2,2) < 0 */


  } /* end for all grids */

/* allocate space for the arrays */
  for (this_boundary = over_ptr->boundary_list->first; this_boundary != NULL;
       this_boundary = this_boundary->next){
    boundary_ptr = this_boundary->data;
    boundary_ptr->coordinates = create_real_array_2d( 2, boundary_ptr->n_points );
  }

/* insert the coordinates in the arrays */
  for (grid_ptr = over_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){

    if ( curve(1,1) > 0 ){
      if ((boundary_ptr = find_curve( over_ptr->boundary_list, curve(1,1) ))){
	i = range(1,1);
	for (j = range(1,2); j <= range(2,2); j++){
	  s = (j-range(1,2)) * grid_ptr->s_step;
	  if (!(grid_ptr->gap_low_r[0] <= s && s <= grid_ptr->gap_low_r[1])){
	    boundary_coord(1,++(boundary_ptr->top)) = x(i,j);
	    boundary_coord(2,boundary_ptr->top) = y(i,j);
	  }
	}
      }
      else
	printf("WARNING: Could not find the boundary part with curve_label = %i\n",
	       curve(1,1));
    } /* end if curve(1,1) > 0 */

    if( curve(2,1) > 0 ){
      if ((boundary_ptr = find_curve( over_ptr->boundary_list, curve(2,1) ))){
	i = range(2,1);
	for (j = range(1,2); j <= range(2,2); j++){
	  s = (j-range(1,2)) * grid_ptr->s_step;
	  if (!(grid_ptr->gap_high_r[0] <= s && s <= grid_ptr->gap_high_r[1])){
	    boundary_coord(1,++(boundary_ptr->top)) = x(i,j);
	    boundary_coord(2,boundary_ptr->top) = y(i,j);
	  }
	}
      }
      else
	printf("WARNING: Could not find the boundary part with curve_label = %i\n",
	       curve(2,1));
    } /* end if curve(2,1) > 0 */

    if( curve(1,2) > 0 ){
      if ((boundary_ptr = find_curve( over_ptr->boundary_list, curve(1,2) ))){
	j = range(1,2);
/* check for c-topology */
	if (grid_ptr->topology == 1){
/* 	  printf("Filling in boundary coordinates a c-grid\n"); */
	  for (i = grid_ptr->first_point; i <= grid_ptr->last_point; i++){
	    boundary_coord(1,++(boundary_ptr->top)) = x(i,j);
	    boundary_coord(2,boundary_ptr->top) = y(i,j);
	  }
/* 	  printf("First point is at (%e,%e)\n",  */
/* 		 x(grid_ptr->first_point,j), y(grid_ptr->first_point,j)); */
/* 	  printf("Last point is at (%e,%e)\n",  */
/* 		 x(grid_ptr->last_point,j), y(grid_ptr->last_point,j)); */
	}
	else{
	  for (i = range(1,1); i <= range(2,1); i++){
	    r = (i-range(1,1)) * grid_ptr->r_step;
	    if (!(grid_ptr->gap_low_s[0] <= r && r <= grid_ptr->gap_low_s[1])){
	      boundary_coord(1,++(boundary_ptr->top)) = x(i,j);
	      boundary_coord(2,boundary_ptr->top) = y(i,j);
	    }
	  }
	}
      }
      else
	printf("WARNING: Could not find the boundary part with curve_label = %i\n",
	       curve(1,2));
    } /* end if curve(1,2) > 0 */

    if( curve(2,2) > 0 ){
      if ((boundary_ptr = find_curve( over_ptr->boundary_list, curve(2,2) ))){
	j = range(2,2);
	for (i = range(1,1); i <= range(2,1); i++){
	  r = (i-range(1,1)) * grid_ptr->r_step;
	  if (!(grid_ptr->gap_high_s[0] <= r && r <= grid_ptr->gap_high_s[1])){
	    boundary_coord(1,++(boundary_ptr->top)) = x(i,j);
	    boundary_coord(2,boundary_ptr->top) = y(i,j);
	  }
	}
      }
      else
	printf("WARNING: Could not find the boundary part with curve_label = %i\n",
	       curve(2,2));
   } /* end if curve(2,2) > 0 */

/* mixed sides */
    if ( curve(1,1) < 0 ){
      if ((boundary_ptr = find_curve( over_ptr->boundary_list, abs(curve(1,1)) ))){
	i = range(1,1);
	for (j = range(1,2); j<= range(2,2); j++){
	  if (left_interp(j)){
	    boundary_coord(1,++(boundary_ptr->top)) = x(i,j);
	    boundary_coord(2,boundary_ptr->top) = y(i,j);
	  }
	}
      }
      else
	printf("WARNING: Could not find the boundary part with curve_label = %i\n",
	       curve(1,1));
    }
    if( curve(2,1) < 0 ){
      if ((boundary_ptr = find_curve( over_ptr->boundary_list, abs(curve(2,1)) ))){
	i = range(2,1);
	for (j = range(1,2); j<= range(2,2); j++){
	  if (right_interp(j)){
	    boundary_coord(1,++(boundary_ptr->top)) = x(i,j);
	    boundary_coord(2,boundary_ptr->top) = y(i,j);
	  }
	}
      }
      else
	printf("WARNING: Could not find the boundary part with curve_label = %i\n",
	       curve(2,1));
    }
    if( curve(1,2) < 0 ){
      if ((boundary_ptr = find_curve( over_ptr->boundary_list, abs(curve(1,2)) ))){
	j = range(1,2);
	for (i = range(1,1); i<= range(2,1); i++){
	  if (lower_interp(i)){
	    boundary_coord(1,++(boundary_ptr->top)) = x(i,j);
	    boundary_coord(2,boundary_ptr->top) = y(i,j);
	  }
	}
      }
      else
	printf("WARNING: Could not find the boundary part with curve_label = %i\n",
	       curve(1,2));
    }
    if( curve(2,2) < 0 ){
      if ((boundary_ptr = find_curve( over_ptr->boundary_list, abs(curve(2,2)) ))){
	j = range(2,2);
	for (i = range(1,1); i<= range(2,1); i++){
	  if (upper_interp(i)){
	    boundary_coord(1,++(boundary_ptr->top)) = x(i,j);
	    boundary_coord(2,boundary_ptr->top) = y(i,j);
	  }
	}
      }
      else
	printf("WARNING: Could not find the boundary part with curve_label = %i\n",
	       curve(2,2));
    }

  } /* end for all grids */

/* sort the coordinates */
  for (this_boundary = over_ptr->boundary_list->first; this_boundary != NULL;
       this_boundary = this_boundary->next){
    boundary_ptr = this_boundary->data;
    if (boundary_ptr->n_points != boundary_ptr->top){
      printf("ERROR: The number of allocated points does not equal the number of\n"
	     "actual point in the global boundary with curve_label = %i\n", 
	     boundary_ptr->curve_label);
    }
    sort_coordinates( boundary_ptr );
  }
  
/* check if the computational domain is inside or outside of the curve */
  for (this_boundary = over_ptr->boundary_list->first; this_boundary != NULL;
       this_boundary = this_boundary->next){
    boundary_ptr = this_boundary->data;
/* find a component grid that has a side with matching curve-label */
    for (grid_ptr = over_ptr->first; grid_ptr != NULL; 
	 grid_ptr = grid_ptr->next){
      if (abs(curve(1,1)) == boundary_ptr->curve_label ||
	  abs(curve(1,2)) == boundary_ptr->curve_label ||
	  abs(curve(2,1)) == boundary_ptr->curve_label ||
	  abs(curve(2,2)) == boundary_ptr->curve_label){
	x_test = x( (range(1,1)+range(2,1))/2 , (range(1,2)+range(2,2))/2 );
	y_test = y( (range(1,1)+range(2,1))/2 , (range(1,2)+range(2,2))/2 );
	break;
      }
    } /* end for all grids */
    if (!grid_ptr)
      printf("ERROR: Could not find a component with curve_label = %i\n", 
	     boundary_ptr->curve_label);
    else{
      boundary_ptr->inside = inside_cutter_curve( x_test, y_test, boundary_ptr );
      printf("The computational domain appears to be %s of curve_label %i\n",
	     (boundary_ptr->inside? "inside": "outside"), boundary_ptr->curve_label);
    }
  }

/* estimate the required mismatch tolerance */
  max_eps = 0.0;
  for (this_boundary = over_ptr->boundary_list->first; this_boundary != NULL;
       this_boundary = this_boundary->next){
    boundary_ptr = this_boundary->data;
    this_eps = 0.0;
    for (i_p=2; i_p<boundary_ptr->n_points; i_p++){
      n_x =  boundary_coord(2,i_p) - boundary_coord(2,i_p-1);
      n_y = -boundary_coord(1,i_p) + boundary_coord(1,i_p-1);
      length = sqrt( n_x*n_x + n_y*n_y );
      if (length > 0.0){
	n_x /= length; n_y /= length;
	this_eps = real_max(this_eps, 
			    fabs(n_x*(boundary_coord(1,i_p+1) - boundary_coord(1,i_p)) +
				 n_y*(boundary_coord(2,i_p+1) - boundary_coord(2,i_p)))
			    );
      } /* end if length > 0 */
    } /* end for all boundary points */
    printf("Boundary part %i. Estimated mismatch tolerance: %e\n", 
	   boundary_ptr->curve_label, this_eps);
    max_eps = real_max( max_eps, this_eps );
  } /* end for all parts of the boundary */
  printf("Estimated global mismatch tolerance: %e\n", max_eps);
  over_ptr->max_n_dist += max_eps;
  printf("Total global mismatch tolerance: %e\n", over_ptr->max_n_dist);
/* done */
  return OK;
}


#define x_h(i,j) compute_index_2d(hole_ptr->x_ptr, i, j)
#define y_h(i,j) compute_index_2d(hole_ptr->y_ptr, i, j)
#define h_flag(i,j) compute_index_2d(hole_ptr->flag_ptr, i, j)
#define h_range(i,j) compute_index_2d(hole_ptr->range_ptr, i, j)

static void
bridge_the_gap(int i0, int j0, real x0, real y0,
	       int i1, int j1, real x1, real y1, 
	       component_grid *hole_ptr, int cutter_curve, overlapping_grid *over_ptr){
  int i_cell, j_cell, i_vertex, j_vertex, i_gap, j_gap, hole_curve;
  real x_new, y_new;
  grid_point gp;

/* compute the true gap */
   if (hole_ptr->r_period) 
     i_gap = int_min( abs(i1-i0), hole_ptr->r_points - 1 - abs(i1-i0)); 
   else 
    i_gap = abs(i1-i0);
   if (hole_ptr->s_period) 
     j_gap = int_min( abs(j1-j0), hole_ptr->s_points - 1 - abs(j1-j0)); 
   else 
    j_gap = abs(j1-j0);

/* check if the gap needs to be filled in further */
  if (i_gap > 1 || j_gap > 1){
    gp.x = x_new = 0.5*(x0 + x1); 
    gp.y = y_new = 0.5*(y0 + y1);
/* check if the middle point (gp.x, gp.y) is enclosed by a cell in the hole grid */
    if (invert_mapping(&gp, 0, 0, cutter_curve, hole_ptr, over_ptr)&&
/* also check if the distance to all boundaries with matching curve_label */
/* exceeds over_ptr->max_n_dist */
      normal_distance(&gp, cutter_curve, hole_ptr, over_ptr)){
/* find the lower left corner of the enclosing cell in the hole grid */
      i_cell = floor( gp.r / hole_ptr->r_step ) + h_range(1,1);
      j_cell = floor( gp.s / hole_ptr->s_step ) + h_range(1,2);
/* avoid index violations */
      i_cell = int_min( i_cell, h_range(2,1)-1 ); 
      i_cell = int_max( i_cell, h_range(1,1) );
      j_cell = int_min( j_cell, h_range(2,2)-1 ); 
      j_cell = int_max( j_cell, h_range(1,2) );

      for (i_vertex = i_cell; i_vertex <= i_cell+1; i_vertex++)
	for (j_vertex = j_cell; j_vertex <= j_cell+1; j_vertex++){
/* set the curve label if the present point is on a real boundary */
	  hole_curve = get_curve(i_vertex, j_vertex, hole_ptr, over_ptr);
	  if ( inside_global_domain(x_h(i_vertex, j_vertex), y_h(i_vertex, j_vertex), 
				    hole_curve, over_ptr ) )
	    h_flag(i_vertex, j_vertex) = 'w';
	  else
	    h_flag(i_vertex, j_vertex) = 'h'; 
	} /* end for */

/* recurse both ways */
      bridge_the_gap(i0, j0, x0, y0, i_cell, j_cell, x_new, y_new, 
		     hole_ptr, cutter_curve, over_ptr);
      bridge_the_gap(i_cell, j_cell, x_new, y_new, i1, j1, x1, y1, 
		     hole_ptr, cutter_curve, over_ptr);
    } /* end if the middle point is inside the hole grid */
  } /* end if there is a gap */
} /* end bridge_the_gap */


static int
normal_distance(grid_point *gp_ptr, int cutter_curve, component_grid *grid_ptr, 
		overlapping_grid *comp_ptr){
  grid_point boundary_gp;
  real dist;

/* Always further than max_n_dist away from all boundaries when curve_label=0 */
  if (cutter_curve == 0) return TRUE;

/* s=0 */
  if ( abs(curve(1,2)) == cutter_curve ){
    boundary_gp.r = gp_ptr->r;
    boundary_gp.s = 0.0;
    forward_grid_mapping( &boundary_gp, grid_ptr );
/* normal distance */
    dist = ((gp_ptr->x - boundary_gp.x) * boundary_gp.yr - 
      (gp_ptr->y - boundary_gp.y) * boundary_gp.xr)/
	sqrt( boundary_gp.xr * boundary_gp.xr + boundary_gp.yr * boundary_gp.yr );
    if (fabs(dist) <= comp_ptr->max_n_dist) 
      return FALSE;
  }
/* s=1 */
  if ( abs(curve(2,2)) == cutter_curve ){
    boundary_gp.r = gp_ptr->r;
    boundary_gp.s = 1.0;
    forward_grid_mapping( &boundary_gp, grid_ptr );
/* normal distance */
    dist = -((gp_ptr->x - boundary_gp.x) * boundary_gp.yr +
      (gp_ptr->y - boundary_gp.y) * boundary_gp.xr)/
	sqrt( boundary_gp.xr * boundary_gp.xr + boundary_gp.yr * boundary_gp.yr );
    if (fabs(dist) <= comp_ptr->max_n_dist) 
      return FALSE;
  }
/* r=0 */
  if ( abs(curve(1,1)) == cutter_curve ){
    boundary_gp.r = 0.0;
    boundary_gp.s = gp_ptr->s;
    forward_grid_mapping( &boundary_gp, grid_ptr );
/* normal distance */
    dist = ((gp_ptr->x - boundary_gp.x) * boundary_gp.ys - 
      (gp_ptr->y - boundary_gp.y) * boundary_gp.xs)/
	sqrt( boundary_gp.xs * boundary_gp.xs + boundary_gp.ys * boundary_gp.ys );
    if (fabs(dist) <= comp_ptr->max_n_dist) 
      return FALSE;
  }
/* r=1 */
  if ( abs(curve(2,1)) == cutter_curve ){
    boundary_gp.r = 1.0;
    boundary_gp.s = gp_ptr->s;
    forward_grid_mapping( &boundary_gp, grid_ptr );
/* normal distance */
    dist = -((gp_ptr->x - boundary_gp.x) * boundary_gp.ys +
      (gp_ptr->y - boundary_gp.y) * boundary_gp.xs)/
	sqrt( boundary_gp.xs * boundary_gp.xs + boundary_gp.ys * boundary_gp.ys );
    if (fabs(dist) <= comp_ptr->max_n_dist) 
      return FALSE;
  }
  
/* further than max_n_dist away from all boundaries with matching curve_label */
  return TRUE;
}

static int
mark_one_cell(real x_c, real y_c, int cutter_curve, component_grid *hole_ptr, 
	      overlapping_grid *over_ptr, int *i_cell, int *j_cell){
  grid_point gp;
  int i_vertex, j_vertex, hole_curve;

  gp.x = x_c;
  gp.y = y_c;
/* invert the mapping to see if (xp,yp) is inside the hole grid */
  if (invert_mapping(&gp, 0, 0, cutter_curve, hole_ptr, over_ptr) &&
/* also check if the distance to all boundaries with matching curve_label */
/* exceeds over_ptr->max_n_dist */
      normal_distance(&gp, cutter_curve, hole_ptr, over_ptr)){
/* find the lower left corner of the cell in the hole grid */
    *i_cell = floor( gp.r / hole_ptr->r_step ) + h_range(1,1);
    *j_cell = floor( gp.s / hole_ptr->s_step ) + h_range(1,2);
/* avoid index violations */
    *i_cell = int_min( *i_cell, h_range(2,1)-1 ); 
    *i_cell = int_max( *i_cell, h_range(1,1) );
    *j_cell = int_min( *j_cell, h_range(2,2)-1 ); 
    *j_cell = int_max( *j_cell, h_range(1,2) );
/* check each vertex in the cell (i_cell, j_cell) in the hole-grid, */
/* to see if it is inside or outside of the computational domain */
    for (i_vertex = *i_cell; i_vertex <= *i_cell+1; i_vertex++)
      for (j_vertex = *j_cell; j_vertex <= *j_cell+1; j_vertex++){
/* set the curve label if the present point is on a real boundary */
	hole_curve = get_curve(i_vertex, j_vertex, hole_ptr, over_ptr);
	if ( inside_global_domain(x_h(i_vertex, j_vertex), y_h(i_vertex, j_vertex), 
				  hole_curve, over_ptr ) )
	  h_flag(i_vertex, j_vertex) = 'w';
	else
	  h_flag(i_vertex, j_vertex) = 'h'; 
      } /* end for */

    return TRUE;
  } /* the point (x_c, y_c) was sufficiently far inside of the grid hole_ptr */
  else
    return FALSE;
} /* end mark_one_cell */

static int
inside_global_domain( real xp, real yp, int curve_label, overlapping_grid *over_ptr ){
  linked_list_member *this_boundary;
  boundary_part *boundary_ptr;
  int inside_this_part;

  for (this_boundary = over_ptr->boundary_list->first;
       this_boundary != NULL; this_boundary = this_boundary->next){
    boundary_ptr = this_boundary->data;
/* a point with matching curve-label is by definition on the */
/* inside of the cutter curve */
    if (curve_label != boundary_ptr->curve_label){
      inside_this_part = inside_cutter_curve( xp, yp, boundary_ptr );
/* the point (xp, yp) is outside of the computational domain if it is */
/* not inside this part, but supposed to be, or */
/* inside this part, but not supposed to be */
      if ((!inside_this_part &&  boundary_ptr->inside) ||
	  (inside_this_part  && !boundary_ptr->inside) )
	return FALSE;
    } /* if NOT matching curve_label */
  }
/* the point is inside of the computational domain if it wasn't */
/* outside of any part of the boundary */
  return TRUE;
}

static int
inside_cutter_curve(real xp, real yp, boundary_part *boundary_ptr){
  int i, intersect=0, n;
  const real eps=1.0e-7;
  real x_star;

  for( i=1; i<boundary_ptr->n_points; i++){
/* is yp between boundary_coord(i) and boundary_coord(i+1) ? */
    if (yp <= real_max(boundary_coord(2,i),boundary_coord(2,i+1)) && 
	yp > real_min(boundary_coord(2,i),boundary_coord(2,i+1)) )
/* is xp to the right of any point? */
      if ( xp > real_min(boundary_coord(1,i),boundary_coord(1,i+1)) )
/* is xp to the right of both points? */
	if (xp >= real_max(boundary_coord(1,i),boundary_coord(1,i+1)) )
	  intersect++;
	else{
/* determine the x-coordinate of the straight line between 
   point (i) and (i+1) */
	  if (fabs(boundary_coord(2,i+1)-boundary_coord(2,i)) > eps){
	    x_star = (boundary_coord(1,i+1)*(yp-boundary_coord(2,i)) + 
		      boundary_coord(1,i)*(boundary_coord(2,i+1)-yp))/
			(boundary_coord(2,i+1)-boundary_coord(2,i));
	    if (x_star <= xp)
	      intersect++;
	  }
	  else
	    intersect++;
	}
  } /* end for all boundary points */

/* check the segment between the last and first point */
  n = boundary_ptr->n_points;
  i = 1;
/* is yp between boundary_coord(i) and boundary_coord(i+1) ? */
  if (yp <= real_max(boundary_coord(2,i),boundary_coord(2,n)) && 
      yp > real_min(boundary_coord(2,i),boundary_coord(2,n)) )
/* is xp to the right of any point? */
    if ( xp > real_min(boundary_coord(1,i),boundary_coord(1,n)) )
/* is xp to the right of both points? */
      if (xp >= real_max(boundary_coord(1,i),boundary_coord(1,n)) )
	intersect++;
      else{
/* determine the x-coordinate of the straight line between the
   points (i) and (n) */
	if (fabs(boundary_coord(2,n)-boundary_coord(2,i)) > eps){
	  x_star = (boundary_coord(1,n)*(yp-boundary_coord(2,i)) + 
		    boundary_coord(1,i)*(boundary_coord(2,n)-yp))/
		      (boundary_coord(2,n)-boundary_coord(2,i));
	  if (x_star <= xp)
	    intersect++;
	}
	else
	  intersect++;
      }

/* The point (xp,yp) was inside if intersect is odd */  
  return ((intersect % 2 == 0)? FALSE : TRUE);
}

/* STEP 2: cut holes in the grids */

static int 
step_2(overlapping_grid *over_ptr){
  component_grid *hole_ptr, *grid_ptr;
  int i, j, hole_curve, currently_inside, i_cell, j_cell, ih, jh
    , go_right, go_up, n, i_min, i_max, j_min, j_max, i_vertex, j_vertex;
  int i_prev=0, j_prev=0;
  boundary_part *boundary_ptr;
  linked_list_member *this_boundary;

  printf("\n");
  for (hole_ptr = over_ptr->first; hole_ptr != NULL; 
       hole_ptr = hole_ptr->next){
    printf("Making holes in grid `%s'\n", hole_ptr->grid_name);
/* go through all parts of the boundary */
    for (this_boundary = over_ptr->boundary_list->first;
	 this_boundary != NULL; this_boundary = this_boundary->next){
      boundary_ptr = this_boundary->data;
      currently_inside = FALSE;
      for (j = 1; j <= boundary_ptr->n_points; j++){

	if (mark_one_cell(boundary_coord(1,j), boundary_coord(2,j), 
			  boundary_ptr->curve_label, 
			  hole_ptr, over_ptr, &i_cell, &j_cell)){
/* check if the current cell is separated from the previous cell */
	  if (currently_inside){
/* subdivide recursively to bridge the gap */
	    bridge_the_gap(i_prev, j_prev, boundary_coord(1,j-1), boundary_coord(2,j-1),
			   i_cell, j_cell, boundary_coord(1,j),   boundary_coord(2,j), 
			   hole_ptr, boundary_ptr->curve_label, over_ptr);
	  }
	  else{
/* the previous point on the hole-cutting boundary was outside, but the present */
/* point was inside, so maybe it is necessary to check some more vertices */
	    if (i_cell == h_range(1,1)){
	      i_min = i_max = h_range(1,1);
	      j_min = int_max(j_cell-1, h_range(1,2));
	      j_max = int_min(j_cell+2, h_range(2,2));
	    }
	    else if (i_cell == h_range(2,1)-1){
	      i_min = i_max = h_range(2,1);
	      j_min = int_max(j_cell-1, h_range(1,2));
	      j_max = int_min(j_cell+2, h_range(2,2));
	    }
	    else if (j_cell == h_range(1,2)){
	      i_min = int_max(i_cell-1, h_range(1,1)); 
	      i_max = int_min(i_cell+2, h_range(2,1));
	      j_min = j_max = h_range(1,2);
	    }
	    else if (j_cell == h_range(2,2)-1){
	      i_min = int_max(i_cell-1, h_range(1,1)); 
	      i_max = int_min(i_cell+2, h_range(2,1));
	      j_min = j_max = h_range(2,2);
	    }
	    else{
/* empty domain */
	      i_min = 1; i_max = 0;
	      j_min = 1; j_max = 0;
	    }
/* check each vertex i_min <= i <= i_max, j_min <= j <= j_max in the hole-grid, */
/* to see if it is inside or outside of the computational domain */
	    for (i_vertex = i_min; i_vertex <= i_max; i_vertex++)
	      for (j_vertex = j_min; j_vertex <= j_max; j_vertex++){
/* set the curve label if the present point is on a real boundary */
		hole_curve = get_curve(i_vertex, j_vertex, hole_ptr, over_ptr);
		if ( inside_global_domain(x_h(i_vertex, j_vertex), 
					  y_h(i_vertex, j_vertex), hole_curve,
					  over_ptr ) )
		  h_flag(i_vertex, j_vertex) = 'w';
		else
		  h_flag(i_vertex, j_vertex) = 'h'; 
	      } /* end for */
	  }
/* remember that we are inside, and where we are */
	  currently_inside = TRUE;
	  i_prev = i_cell;
	  j_prev = j_cell;
	} /* end if mark_one_cell(boundary_coord(j), hole_ptr) */
	else{
/* if the present point on the hole-cutting boundary was outside, but the previous */
/* point was inside, it might be necessary to back up */
	  if (currently_inside){
	    if (i_prev == h_range(1,1)){
	      i_min = i_max = h_range(1,1);
	      j_min = int_max(j_prev-1, h_range(1,2));
	      j_max = int_min(j_prev+2, h_range(2,2));
	    }
	    else if (i_prev == h_range(2,1)-1){
	      i_min = i_max = h_range(2,1);
	      j_min = int_max(j_prev-1, h_range(1,2));
	      j_max = int_min(j_prev+2, h_range(2,2));
	    }
	    else if (j_prev == h_range(1,2)){
	      i_min = int_max(i_prev-1, h_range(1,1)); 
	      i_max = int_min(i_prev+2, h_range(2,1));
	      j_min = j_max = h_range(1,2);
	    }
	    else if (j_prev == h_range(2,2)-1){
	      i_min = int_max(i_prev-1, h_range(1,1)); 
	      i_max = int_min(i_prev+2, h_range(2,1));
	      j_min = j_max = h_range(2,2);
	    }
	    else{
/* empty domain */
	      i_min = 1; i_max = 0;
	      j_min = 1; j_max = 0;
	    }
/* check each vertex i_min <= i <= i_max, j_min <= j <= j_max in the hole-grid, */
/* to see if it is inside or outside of the computational domain */
	    for (i_vertex = i_min; i_vertex <= i_max; i_vertex++)
	      for (j_vertex = j_min; j_vertex <= j_max; j_vertex++){
/* set the curve label if the present point is on a real boundary */
		hole_curve = get_curve(i_vertex, j_vertex, hole_ptr, over_ptr);
		if ( inside_global_domain(x_h(i_vertex, j_vertex), 
					  y_h(i_vertex, j_vertex), hole_curve,
					  over_ptr ) )
		  h_flag(i_vertex, j_vertex) = 'w';
		else
		  h_flag(i_vertex, j_vertex) = 'h'; 
	      } /* end for */
	  }
	  currently_inside = FALSE;
	}
      } /* end for all points on the present hole-cutting boundary */

/* Should also bridge the gap for the segment between the last and first points */
      if (currently_inside){
	j = 1;
	n = boundary_ptr->n_points;
	if (mark_one_cell(boundary_coord(1,j), boundary_coord(2,j), 
			  boundary_ptr->curve_label, 
			  hole_ptr, over_ptr, &i_cell, &j_cell)){
/* subdivide recursively to bridge the gap */
	  bridge_the_gap(i_prev, j_prev, boundary_coord(1,n), boundary_coord(2,n),
			 i_cell, j_cell, boundary_coord(1,j), boundary_coord(2,j), 
			 hole_ptr, boundary_ptr->curve_label, over_ptr);
	} /* end if mark_one_cell(boundary_point(1), hole_ptr) */
	else{
/* the present point on the hole-cutting boundary is outside, but the previous */
/* point was inside, so it might be necessary to back up */
	  if (i_prev == h_range(1,1)){
	    i_min = i_max = h_range(1,1);
	    j_min = int_max(j_prev-1, h_range(1,2));
	    j_max = int_min(j_prev+2, h_range(2,2));
	  }
	  else if (i_prev == h_range(2,1)-1){
	    i_min = i_max = h_range(2,1);
	    j_min = int_max(j_prev-1, h_range(1,2));
	    j_max = int_min(j_prev+2, h_range(2,2));
	  }
	  else if (j_prev == h_range(1,2)){
	    i_min = int_max(i_prev-1, h_range(1,1)); 
	    i_max = int_min(i_prev+2, h_range(2,1));
	    j_min = j_max = h_range(1,2);
	  }
	  else if (j_prev == h_range(2,2)-1){
	    i_min = int_max(i_prev-1, h_range(1,1)); 
	    i_max = int_min(i_prev+2, h_range(2,1));
	    j_min = j_max = h_range(2,2);
	  }
	  else{
/* empty domain */
	    i_min = 1; i_max = 0;
	    j_min = 1; j_max = 0;
	  }
/* check each vertex i_min <= i <= i_max, j_min <= j <= j_max in the hole-grid, */
/* to see if it is inside or outside of the computational domain */
	  for (i_vertex = i_min; i_vertex <= i_max; i_vertex++)
	    for (j_vertex = j_min; j_vertex <= j_max; j_vertex++){
/* set the curve label if the present point is on a real boundary */
	      hole_curve = get_curve(i_vertex, j_vertex, hole_ptr, over_ptr);
	      if ( inside_global_domain(x_h(i_vertex, j_vertex), 
					y_h(i_vertex, j_vertex), hole_curve,
					over_ptr ) )
		h_flag(i_vertex, j_vertex) = 'w';
	      else
		h_flag(i_vertex, j_vertex) = 'h'; 
	    } /* end for */
	} /* end else... */
      } /* end if currently_inside */

    } /* end for all hole-cutting boundaries */
  } /* end for all grids */

#ifdef STOP_AFTER_HOLE_CUT_1
/* plot the overlapping grid */
  PL_erase(0);
  PL_start_plot(0);
  plot_overlapping( over_ptr, 1+8 ); /* the used grid points and the flag array */
  PL_end_plot();
  wait_for_key("You are looking at the flag array after the first stage in the hole "
	       "cut process.\n"
	       "Hit RETURN to proceed>");
#endif

/* STEP 2 */
/* Fill in the holes */
  printf("\n");
  for (grid_ptr = over_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){
    printf("Filling in the holes in grid `%s'\n", grid_ptr->grid_name);
    for (j=range(1,2); j<=range(2,2); j++)
      for (i=range(1,1); i<=range(2,1); i++){
/* check if it is a hole point */
	if (flag(i,j) == 'h'){
/* check the neighbors */
	  go_right = (i == range(1,1) || flag(i-1,j) == 'w' || flag(i-1,j) == 'h');
	  go_up = (j == range(1,2) || flag(i,j-1) == 'w' || flag(i,j-1) == 'h');
/* horizontal sweep */
	  if (go_right){
	    for (ih = i+1; ih <= range(2,1) && flag(ih,j) != 'h' && flag(ih,j) != 'w';
		 ih++)
	      flag(ih,j) = 0;
	  }
	  else{
	    for (ih = i-1; ih >= range(1,1) && flag(ih,j) != 'h' && flag(ih,j) != 'w';
		 ih--)
	      flag(ih,j) = 0;
	  }
/* vertical sweep */
	  if (go_up){
	    for (jh = j+1; jh <= range(2,2) && flag(i,jh) != 'h' && flag(i,jh) != 'w';
		 jh++)
	      flag(i,jh) = 0;
	  }
	  else{
	    for (jh = j-1; jh >= range(1,2) && flag(i,jh) != 'h' && flag(i,jh) != 'w';
		 jh--)
	      flag(i,jh) = 0;
	  }
	} /* end if 'h' point */
      } /* end for all grid points */

/* change 'w' to the  grid priority and 'h' to 0 */
    for (j=range(1,2); j<=range(2,2); j++)
      for (i=range(1,1); i<=range(2,1); i++){
	if (flag(i,j) == 'w')
	  flag(i,j) = grid_ptr->priority;
	else if (flag(i,j) == 'h')
	  flag(i,j) = 0;
      } /* end for all grid points */
  } /* end for all grids */


/* Mixed interpolation / physical boundaries */

  printf("\n");
  for (grid_ptr = over_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){
    printf("Looking for interpolation points on mixed "
	   "boundaries in grid `%s'\n", grid_ptr->grid_name);

/* left side */

    i = range(1,1);
    for (j = range(1,2); j<= range(2,2); j++){
      if (curve(1,1) < 0 && boundary_interp( x(i,j), y(i,j), grid_ptr, over_ptr)){
	left_interp(j) = curve(1,1);
      }
      else{
	left_interp(j) = abs(curve(1,1));
      }
    }

/* right side */

    i = range(2,1);
    for (j = range(1,2); j<= range(2,2); j++){
      if (curve(2,1) < 0 && boundary_interp( x(i,j), y(i,j), grid_ptr, over_ptr)){
	right_interp(j) = curve(2,1);
      }
      else{
	right_interp(j) = abs(curve(2,1));
      }
    }

/* lower side */

    j = range(1,2);
    for (i = range(1,1); i<= range(2,1); i++){
      if (curve(1,2) < 0 && boundary_interp( x(i,j), y(i,j), grid_ptr, over_ptr)){
	lower_interp(i) = curve(1,2);
      }
      else{
	lower_interp(i) = abs(curve(1,2));
      }
    }

/* upper side */

    j = range(2,2);
    for (i = range(1,1); i<= range(2,1); i++){
      if (curve(2,2) < 0 && boundary_interp( x(i,j), y(i,j), grid_ptr, over_ptr)){
	upper_interp(i) = curve(2,2);
      }
      else{
	upper_interp(i) = abs(curve(2,2));
      }
    }

  } /* end for all component grids */

/* done */
  return OK;
}


#define sq(x) ((x)*(x))
#define curve_d(i,j) compute_index_2d( donor_grid->curve_ptr, i, j)

static int
closest_boundary(int present_curve, component_grid *grid_ptr, 
		 interp_point *interp_ptr, real *x_dist, real *y_dist){
  real x_dist_tmp, y_dist_tmp;
  int found_closest_bndry = FALSE;
  grid_point gp;
  component_grid *donor_grid;

/* mixed physical / interpolating curves are not of interest here! */
  if (present_curve <= 0 ) return FALSE;

  donor_grid = interp_ptr->grid_loc;
  *x_dist = *y_dist = 1.e30;
/* low r in donor grid */
  if (abs(curve_d(1,1)) == abs(present_curve)){
    gp.r = 0.0;
    gp.s = interp_ptr->s_loc;
    forward_grid_mapping( &gp, donor_grid );
    *x_dist = gp.x - x(interp_ptr->i_point, interp_ptr->j_point);
    *y_dist = gp.y - y(interp_ptr->i_point, interp_ptr->j_point);
    found_closest_bndry = TRUE;
  }
/* high r in donor grid */
  if (abs(curve_d(2,1)) == abs(present_curve)){
    gp.r = 1.0;
    gp.s = interp_ptr->s_loc;
    forward_grid_mapping( &gp, donor_grid );
    x_dist_tmp = gp.x - x(interp_ptr->i_point, interp_ptr->j_point);
    y_dist_tmp = gp.y - y(interp_ptr->i_point, interp_ptr->j_point);
    if ((sq(x_dist_tmp) + sq(y_dist_tmp)) < (sq(*x_dist) + sq(*y_dist))){
      found_closest_bndry = TRUE;
      *x_dist = x_dist_tmp;
      *y_dist = y_dist_tmp;
    }	    
  }
/* low s in donor grid */
  if (abs(curve_d(1,2)) == abs(present_curve)){
    gp.r = interp_ptr->r_loc;
    gp.s = 0.0;
    forward_grid_mapping( &gp, donor_grid );
    x_dist_tmp = gp.x - x(interp_ptr->i_point, interp_ptr->j_point);
    y_dist_tmp = gp.y - y(interp_ptr->i_point, interp_ptr->j_point);
    if ((sq(x_dist_tmp) + sq(y_dist_tmp)) < (sq(*x_dist) + sq(*y_dist))){
      found_closest_bndry = TRUE;
      *x_dist = x_dist_tmp;
      *y_dist = y_dist_tmp;
    }	    
  }
/* high s in donor grid */
  if (abs(curve_d(2,2)) == abs(present_curve)){
    gp.r = interp_ptr->r_loc;
    gp.s = 1.0;
    forward_grid_mapping( &gp, donor_grid );
    x_dist_tmp = gp.x - x(interp_ptr->i_point, interp_ptr->j_point);
    y_dist_tmp = gp.y - y(interp_ptr->i_point, interp_ptr->j_point);
    if ((sq(x_dist_tmp) + sq(y_dist_tmp)) < (sq(*x_dist) + sq(*y_dist))){
      found_closest_bndry = TRUE;
      *x_dist = x_dist_tmp;
      *y_dist = y_dist_tmp;
    }	    
  }
  return found_closest_bndry;
}

static void
recompute_interpolation_location(interp_point *interp_ptr, component_grid *grid_ptr, 
				 overlapping_grid *comp_ptr, int job, real x_dist, 
				 real y_dist){
  grid_point gp;
  interp_point *interp2_ptr, new_interp;
  int present_curve;

/* tmp */
  printf("The mismatch at boundary ip (%i,%i) is (%e,%e)\n", 
	 interp_ptr->i_point, interp_ptr->j_point, x_dist, y_dist);

  for (interp2_ptr = grid_ptr->last_interp; interp2_ptr != NULL;
       interp2_ptr = interp2_ptr->prev){
    if ((job == 1 && interp2_ptr->i_point == interp_ptr->i_point) ||
	(job == 2 && interp2_ptr->j_point == interp_ptr->j_point)){
/* recompute */
      gp.x = x(interp2_ptr->i_point, interp2_ptr->j_point) + x_dist; 
      gp.y = y(interp2_ptr->i_point, interp2_ptr->j_point) + y_dist; 
/* set the curve label if the present point is on a real boundary */
      present_curve = get_curve(interp2_ptr->i_point, interp2_ptr->j_point, 
				grid_ptr, comp_ptr);
      if (invert_mapping(&gp, 0, 0, present_curve, interp2_ptr->grid_loc, 
			 comp_ptr)){
/* save the (r,s) coordinates */
	new_interp.r_loc = gp.r;
	new_interp.s_loc = gp.s;
/* check if (i,j) is a valid interpolation location*/
 	if (good_interp_loc(&new_interp, interp2_ptr->grid_loc, comp_ptr)){
/* update the coefficients */
	  interp2_ptr->r_loc = new_interp.r_loc;
	  interp2_ptr->s_loc = new_interp.s_loc;
	  interp2_ptr->i_loc = new_interp.i_loc;
	  interp2_ptr->j_loc = new_interp.j_loc;
	}
	else
	  printf("Warning, the interpolation location for ip (%i,%i) was not "
		 "updated because the new location was bad\n", 
		 interp2_ptr->i_point, interp2_ptr->j_point);
      }
      else
	printf("Warning, the interpolation location for ip (%i,%i) was not "
	       "updated because the new coordinate was outside of the donor "
	       "grid.\n", interp2_ptr->i_point, interp2_ptr->j_point);
    } /* end if along the same line */
  } /* end for all interpolation points */
}

/* STEP 3: classify all gridpoints as either interpolation points , 
   discretization points or exterior points */

static int 
step_3(input_output *io_ptr, overlapping_grid *comp_ptr){
  component_grid *grid_ptr, *other_ptr;
  int i, j, interpolee, iw1, ok, present_curve, bad_disc, bad_interp, n_interp
    , ok_interp, bad_explicit, ok_explicit;
  real x_dist, y_dist;
  interp_point new_interp, *interp_ptr;
#if (defined(STOP_AFTER_INTIAL_CLASS) || defined(STOP_AFTER_EXPLICIT))
  int icom, col_flag;
  FILE *fp;
  char * break_com[3], *file_name;
  const int LEVEL=1;
  break_com[0] = "proceed";
  break_com[1] = "break";
  break_com[2] = "postscript-copy";
#endif

  printf("\n");
  for (grid_ptr = comp_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){
    printf("Classifying all points in grid `%s'\n",grid_ptr->grid_name);
    for (i = range(1,1); i <= range(2,1); i++){

/* the first point in each column has no initial guess */
      new_interp.r_loc = 0.0; new_interp.s_loc = 0.0;
      new_interp.grid_loc = NULL;
      for (j = range(1,2); j <= range(2,2); j++){
/* don't mess with the points that were flagged exterior in step 2 */
	if (flag(i,j) != 0){
/* save the index of the present point */
	  new_interp.i_point = i;
	  new_interp.j_point = j;
/* set the curve label if the present point is on a real boundary */
	  present_curve = get_curve(i, j, grid_ptr, comp_ptr);

/* test if (i,j) can interpolalate from a higher grid */
	  interpolee = interp_from_higher( x(i,j), y(i,j), present_curve,
					  &new_interp, grid_ptr, comp_ptr);
	  if (interpolee > 0)
	    flag(i,j) = interpolee;
	  else if ( disc_point( i, j, grid_ptr, comp_ptr) )
	    flag(i,j) = grid_ptr->priority;
	  else
	    flag(i,j) = interp_from_lower( x(i,j), y(i,j), present_curve,
					  &new_interp, grid_ptr, grid_ptr, comp_ptr);
	  if (flag(i,j) != grid_ptr->priority && flag(i,j) != 0){
	    grid_ptr->last_interp = 
	      create_interp_point(&new_interp, grid_ptr->last_interp, grid_ptr);
	  }
	  else{
/* reset the initial guess */
	    new_interp.r_loc = 0.0;
	    new_interp.s_loc = 0.0;
	    new_interp.grid_loc = NULL;
	  }
	}
      }
    }
  } /* end for all grids */
  
#ifdef STOP_AFTER_INITIAL_CLASS
  PL_erase(0);
  PL_start_plot(0);
  plot_overlapping( comp_ptr, 1+2 ); /* used grid points + interpolation points */
  PL_end_plot();
/* break or proceed? */    
  while ((icom=get_command(io_ptr, "After the initial classification>", 
			   break_com, 2, LEVEL, NULL, NULL)) != 0)
    {
      switch (icom)
	{
	case 1:			/* exit to calling routine */
	  return FALSE;
	  break;
	case 2:			/* postscript copy */
	  if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
				     "classify.ps", &file_name, 'w', 0, 
				     FALSE) ) != NULL)
	    {
	      col_flag = get_yes_no( io_ptr, "Do you want to save color "
				    "information>", LEVEL );
	      PL_postscript( 0, fp, col_flag );
	      fclose( fp );
	    }
	  break;
	default:
	  ;
	} /* end switch */
    } /* end while */
#endif

  if (comp_ptr->correct_mismatch){
    printf("\n");
    for (grid_ptr = comp_ptr->first; grid_ptr != NULL; 
	 grid_ptr = grid_ptr->next){
      printf("Mismatch correction in grid `%s'\n",grid_ptr->grid_name);
      for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	   interp_ptr = interp_ptr->prev){

/* check if we are at a boundary interpolation point */
/* low r */
	if (interp_ptr->i_point == range(1,1) && 
	    left_interp(interp_ptr->j_point) > 0 &&
	    closest_boundary(curve(1,1), grid_ptr, interp_ptr, &x_dist, &y_dist)){
	  recompute_interpolation_location(interp_ptr, grid_ptr, comp_ptr, 2, 
					   x_dist, y_dist);
	}
/* high r */
	else if (interp_ptr->i_point == range(2,1) && 
		 right_interp(interp_ptr->j_point) > 0 &&
		 closest_boundary(curve(2,1), grid_ptr, interp_ptr, &x_dist, &y_dist)){
	  recompute_interpolation_location(interp_ptr, grid_ptr, comp_ptr, 2, 
					   x_dist, y_dist);
	}
/* low s */
	if (interp_ptr->j_point == range(1,2) && 
	    lower_interp(interp_ptr->i_point) > 0 &&
	    closest_boundary(curve(1,2), grid_ptr, interp_ptr, &x_dist, &y_dist)){
	  recompute_interpolation_location(interp_ptr, grid_ptr, comp_ptr, 1, 
					   x_dist, y_dist);
	} /* end low s */
/* high s */
	else if (interp_ptr->j_point == range(2,2) && 
		 upper_interp(interp_ptr->i_point) > 0 &&
		 closest_boundary(curve(2,2), grid_ptr, interp_ptr, &x_dist, &y_dist)){
	  recompute_interpolation_location(interp_ptr, grid_ptr, comp_ptr, 1, 
					   x_dist, y_dist);
	}

      } /* end for all interpolation points */
    } /* end for all component grids */
  } /* end if correct_mismatch */

  /* ensure explicit interpolation, if wanted */
  if (comp_ptr->interp_type == 'e'){
    explicit(io_ptr, comp_ptr);
#ifdef STOP_AFTER_EXPLICIT
  PL_erase(0);
  PL_start_plot(0);
  plot_overlapping( comp_ptr, 1+2 ); /* used grid points + interpolation points */
  PL_end_plot();
/* break or proceed? */    
  while ((icom=get_command(io_ptr, "After explicit interpolation>", 
			   break_com, 2, LEVEL, NULL, NULL)) != 0)
    {
      switch (icom)
	{
	case 1:			/* exit to calling routine */
	  return FALSE;
	  break;
	case 2:			/* postscript copy */
	  if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
				     "classify.ps", &file_name, 'w', 0, 
				     FALSE) ) != NULL)
	    {
	      col_flag = get_yes_no( io_ptr, "Do you want to save color "
				    "information>", LEVEL );
	      PL_postscript( 0, fp, col_flag );
	      fclose( fp );
	    }
	  break;
	default:
	  ;
	} /* end switch */
    } /* end while */
#endif
  }
  
  /* now check the consistency of the flag array */
  iw1 = comp_ptr->interp_width - 1;
  ok = 1;
  printf("\nChecking the consistency of all grid points...\n");
  for (grid_ptr = comp_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){
    /* interior points */
    bad_disc = 0;
    for (i = range(1,1); i <= range(2,1); i++)
      for (j = range(1,2); j <= range(2,2); j++){
	/* set the curve label if the present point is on a real boundary */
	present_curve = get_curve( i, j, grid_ptr, comp_ptr);

	if (flag(i,j) == grid_ptr->priority && 
	    !disc_point( i, j, grid_ptr, comp_ptr)){
	  insert_bad_point(i, j, BOX, grid_ptr);
	  bad_disc++;
	}
      }
    /* interpolation points */
    bad_interp = 0;
    bad_explicit = 0;
    n_interp = 0;
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL; 
	 interp_ptr = interp_ptr->prev){
      if (interp_ptr->active){
	n_interp++;
	other_ptr = interp_ptr->grid_loc;
	ok_interp = 1;
	ok_explicit = 1;
	for (i = interp_ptr->i_loc; i <= interp_ptr->i_loc+iw1; i++)
	  for (j = interp_ptr->j_loc; j <= interp_ptr->j_loc+iw1; j++){
	    if (other_flag(i,j) == 0){
	      insert_bad_point(i, j, CROSS, other_ptr);
	      ok_interp = 0;
	    }
	    /* check for non-explicit interpolation */
	    else if (comp_ptr->interp_type == 'e' &&
		     abs(other_flag(i,j)) != other_ptr->priority){
	      insert_bad_point(i, j, FILLED_CIRCLE, other_ptr);
	      ok_explicit = 0;
	    }
	  }	    
	if (!ok_interp){
	  insert_bad_point(interp_ptr->i_point,interp_ptr->j_point,
			   TRIANGLE,grid_ptr);
	  bad_interp++;
	}
	else if (!ok_explicit){
	  insert_bad_point(interp_ptr->i_point,interp_ptr->j_point,
			   INV_TRIANGLE, grid_ptr);
	  bad_explicit++;
	}
      }
      else{
	printf("A deactivated interpolation point found while checking ");
	printf("the constistency\n");
      }
    }
    /* report if there are any errors */
    if (bad_disc > 0 || bad_interp > 0 || bad_explicit > 0){
      ok = 0;
      printf("\nInconsistencies found in grid `%s'.\n", grid_ptr->grid_name);
      if (bad_disc>0)
	printf("It has %i bad discretization points\n", bad_disc);
      if (bad_interp>0)
	printf("It has %i bad interpolation points because of dead donors\n", 
	       bad_interp);
      if (comp_ptr->interp_type == 'e' && bad_explicit>0)
	printf("It has %i non-explicit interpolation points\n", bad_explicit);
    }
    
    if (n_interp != grid_ptr->n_interp)
      printf("n_interp missmatch in grid `%s'\n",grid_ptr->grid_name);
  }
  
  
  if (ok) printf("All discretization and interpolation points OK.\n");
  return ok;
}



/* mark interpolation points which are needed for interpolation by higher 
   priority grids */

static int 
step_4(overlapping_grid *comp_ptr){
  int iw1,i,j;
  component_grid *other_ptr,*grid_ptr;
  interp_point *this_interp_ptr;

  /* only necessary to mark points when the interpolation type is implicit */
  if (comp_ptr->interp_type == 'e')
    return 1;

  iw1 = comp_ptr->interp_width - 1;
  for (other_ptr = comp_ptr->first; other_ptr != NULL; 
       other_ptr = other_ptr->next){
    for(this_interp_ptr = other_ptr->last_interp; this_interp_ptr != NULL;
	this_interp_ptr = this_interp_ptr->prev)
      /* check if this interpolation point interpolates from a grid with lower
	 priority than other_ptr->priority */
      if (this_interp_ptr->grid_loc->priority < other_ptr->priority){
	/* set grid_ptr to point to the interpolation location, so that the
	   inline macros for flag will work */
	grid_ptr = this_interp_ptr->grid_loc;
	/* change sign of all the points in that grid that are needed in that 
	   interpolation formula */
	for (i=this_interp_ptr->i_loc; i<=this_interp_ptr->i_loc + iw1; i++)
	  for (j=this_interp_ptr->j_loc; j<=this_interp_ptr->j_loc + iw1; j++)
	    flag(i,j) = -abs(flag(i,j));
      }
  }
  return 1;
}


/* delete interpolation points which are not needed and change interpolation 
   points into discretization points */

static int 
step_5(overlapping_grid *comp_ptr){
  component_grid *grid_ptr, *other_grid;
  interp_point *interp_ptr;
  int iw1, points_removed, i, j;
  
  iw1 = comp_ptr->interp_width - 1;
  
  for (grid_ptr = comp_ptr->last; grid_ptr != NULL; 
       grid_ptr = grid_ptr->prev){
    printf("\nTrimming grid `%s'.\n", grid_ptr->grid_name);
    
/* maximize the width of the overlap */
/* ***************************************** */
/*   trim-style == 1 starts here               */
/* ***************************************** */
    if (comp_ptr->trim_style == 1){
/* if this interpolation point can be an interior point instead, 
   set flag(i_point,j_point) = grid_ptr->priority and delete this 
   interpolation point from the stack */
      points_removed = 0;

      for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL; 
	   interp_ptr = interp_ptr->prev){
      
	if (interp_ptr->active &&
	    disc_point(interp_ptr->i_point, interp_ptr->j_point, grid_ptr, 
		       comp_ptr)){
	
/* change the flag and remove this interpolation point from the stack */
	
	  flag(interp_ptr->i_point,interp_ptr->j_point) = grid_ptr->priority;
	  interp_ptr->active = 0;
	  points_removed++;
	}
      }

/* update the number of interpolation points */
      grid_ptr->n_interp += -points_removed;
      if (points_removed > 0)
	printf("Changed %i interpolation points to discretization points\n", 
	       points_removed);
    } /* end if trim-style == 1 */
/* ***************************************** */
/*   trim-style == 1 ends here               */
/* ***************************************** */


/* search all interpolation points in all other grids for donors in the */
/* present grid */
    if (comp_ptr->interp_type == 'i'){
      for (other_grid = comp_ptr->first; other_grid != NULL;
	   other_grid = other_grid->next){
	if (other_grid == grid_ptr) continue;
	for (interp_ptr = other_grid->last_interp; interp_ptr != NULL;
	     interp_ptr = interp_ptr->prev){
	  if (interp_ptr->grid_loc == grid_ptr){
/* mark the donor points */
	    for (i=interp_ptr->i_loc; i<=interp_ptr->i_loc+iw1; i++)
	      for (j=interp_ptr->j_loc; j<=interp_ptr->j_loc+iw1; j++)
		flag(i,j) = -abs(flag(i,j));
	  } /* end if donor is grid_ptr */
	} /* end for all interpolation points */
      } /* end for al grids except grid_ptr */
    } /* end if implicit interpolation */

    points_removed = 0;
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL; 
	 interp_ptr = interp_ptr->prev){
      
/* don't mess with unused or marked points, i.e. points where flag <= 0 */
      if (flag(interp_ptr->i_point,interp_ptr->j_point) > 0){

/* if this interpolation point is not needed, set flag(i_point,j_point) = 0
   and delete this interpolation point from the stack */

	if (!needed_by_disc(interp_ptr->i_point, interp_ptr->j_point, 
			    grid_ptr, comp_ptr) 
/* 	    &&  */
/* 	    !needed_by_interp(interp_ptr->i_point, interp_ptr->j_point,  */
/* 			      grid_ptr->priority, comp_ptr)  */
	    ){
/* remove this interpolation point */
	  flag(interp_ptr->i_point,interp_ptr->j_point) = 0;
	  interp_ptr->active = 0;
	  points_removed++;
	}
      }
    } /* end for all interpolation points */

/* update the number of interpolation points */
    grid_ptr->n_interp += -points_removed;
    if (points_removed > 0) 
      printf("%i interpolation points not needed\n", points_removed);

/* minimize the width of the overlap */
    if (comp_ptr->trim_style == 0){
/* if this interpolation point can be an interior point instead, 
   set flag(i_point,j_point) = grid_ptr->priority and delete this 
   interpolation point from the stack */
      points_removed = 0;

      for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL; 
	   interp_ptr = interp_ptr->prev){
      
	if (interp_ptr->active &&
	    disc_point(interp_ptr->i_point, interp_ptr->j_point, grid_ptr, 
		       comp_ptr)){
	
/* change the flag and remove this interpolation point from the stack */
	
	  flag(interp_ptr->i_point,interp_ptr->j_point) = grid_ptr->priority;
	  interp_ptr->active = 0;
	  points_removed++;
	}
      } /* end for all interpolation points */

/* update the number of interpolation points */
      grid_ptr->n_interp += -points_removed;
      if (points_removed > 0)
	printf("Changed %i interpolation points to discretization points\n", 
	       points_removed);
    }

    /* remove the de-activated interpolation points */
    remove_deactive_interp( grid_ptr );

    /* the marking is only necessary when the interpolation type is implicit */
    if (comp_ptr->interp_type == 'i'){
      for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	   interp_ptr = interp_ptr->prev)
	if (interp_ptr->grid_loc->priority > grid_ptr->priority)
	  /* set flag = - |flag| for all points in the interpolation formula */
	  change_sign(interp_ptr->i_loc, interp_ptr->j_loc, iw1,
		      interp_ptr->grid_loc);
    }

  } /* end for all grids */
  return 0;
}


/* change sign of the entries in the flag array so that discretization points 
   are positive and interpolation points are negative */

static int 
step_6(overlapping_grid *comp_ptr){
  component_grid *grid_ptr;
  interp_point *interp_ptr;
  int i, j, n_interp;

  for (grid_ptr = comp_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){

/* update the sign of the flag array */
    for (i = range(1,1); i <= range(2,1); i++)
      for (j = range(1,2); j <= range(2,2); j++)
	if (abs(flag(i,j)) == grid_ptr->priority)
	  flag(i,j) = abs(flag(i,j));
	else if (flag(i,j) != 0)
	  flag(i,j) = -abs(flag(i,j));

/* take care of the periodic case */
    if (grid_ptr->r_period){
/* remove any spurious signs from the flag array for the periodic points */
      for (j = 1; j <= grid_ptr->s_dim; j++){
	for (i = 1; i <= range(1,1)-1; i++)
	  if (abs(flag(i,j)) == grid_ptr->priority)
	    flag(i,j) = abs(flag(i,j));
	  else if (flag(i,j) != 0)
	    flag(i,j) = -abs(flag(i,j));
	for (i = grid_ptr->r_dim; i >= range(2,1); i--)
	  if (abs(flag(i,j)) == grid_ptr->priority)
	    flag(i,j) = abs(flag(i,j));
	  else if (flag(i,j) != 0)
	    flag(i,j) = -abs(flag(i,j));
      }

/* remove any interpolation points that have i_point=range(2,1) */
      for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	   interp_ptr = interp_ptr->prev){

	if (interp_ptr->i_point == range(2,1)){
/* change the flag and remove this interpolation point from the stack */
	  flag(interp_ptr->i_point,interp_ptr->j_point) = grid_ptr->priority;
	  interp_ptr->active = 0;  
/*	  printf("\nRemoved an interpolation point with i_point == " */
/*		 "range(2,1) in grid %s.\n", grid_ptr->grid_name );  */
	  grid_ptr->n_interp--;
	}
      }
    }
    if (grid_ptr->s_period){
/* remove any spurious signs from the flag array for the periodic points */
      for (i = 1; i <= grid_ptr->r_dim; i++){
	for (j = 1; j <= range(1,2)-1; j++)
	  if (abs(flag(i,j)) == grid_ptr->priority)
	    flag(i,j) = abs(flag(i,j));
	  else if (flag(i,j) != 0)
	    flag(i,j) = -abs(flag(i,j));
	for (j = grid_ptr->s_dim; j >= range(2,2); j--)
	  if (abs(flag(i,j)) == grid_ptr->priority)
	    flag(i,j) = abs(flag(i,j));
	  else if (flag(i,j) != 0)
	    flag(i,j) = -abs(flag(i,j));
      }

/* remove any interpolation points that have j_point=range(2,2) */
      for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	   interp_ptr = interp_ptr->prev){

	if (interp_ptr->j_point == range(2,2)){
/* change the flag and remove this interpolation point from the stack */
	  flag(interp_ptr->i_point,interp_ptr->j_point) = grid_ptr->priority;
	  interp_ptr->active = 0;
/*	  printf("\nRemoved an interpolation point with j_point == " */
/*		 "range(2,2) in grid %s.\n", grid_ptr->grid_name );  */
	  grid_ptr->n_interp--;
	}
      }
    }

/* remove any de-activated interpolation points */
    remove_deactive_interp( grid_ptr );

/* check the number of interpolation points */
    n_interp = 0;
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev)
      n_interp++;
    if (n_interp != grid_ptr->n_interp){
      printf("Number of interpolation points missmatch in grid `%s'\n"
	     ,grid_ptr->grid_name);
      printf("grid_ptr->n_interp = %i, while the actual numer was %i\n"
	     ,grid_ptr->n_interp, n_interp);
    }

  } /* end for all grids */

  return 0;
}

/* remove any interpolation points with zero flag */
static int 
step_7(overlapping_grid *comp_ptr){
  component_grid *grid_ptr;
  interp_point *interp_ptr;

  for (grid_ptr = comp_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL;
	 interp_ptr = interp_ptr->prev){
      if (flag(interp_ptr->i_point, interp_ptr->j_point) == 0){
	interp_ptr->active = 0;
	grid_ptr->n_interp--;
      }
    }
/* remove any de-activated interpolation points */
    remove_deactive_interp( grid_ptr );
  }
  return OK;
}

static int 
needed_by_disc(int i_point, int j_point, component_grid *grid_ptr,
			  overlapping_grid *comp_ptr){
/* check if it is needed by a discretization point in the same grid */
  int i_min, j_min, i_max, j_max, i, j, dw2, nw1, tw2, cw1;

  dw2 = (comp_ptr->disc_width - 1)/2;

/* needed by an interior point? */
  i_min = int_max(i_point - dw2, range(1,1));
  i_max = int_min(i_point + dw2, range(2,1));
  j_min = int_max(j_point - dw2, range(1,2));
  j_max = int_min(j_point + dw2, range(2,2));

  for (i=i_min; i<= i_max; i++)
    for (j=j_min; j<= j_max; j++)
      if (abs(flag(i,j)) == grid_ptr->priority)
	return 1;

/* check if the point is needed by a boundary formula */
  nw1 = comp_ptr->normal_width-1;
  tw2 = (comp_ptr->tangent_width-1)/2;

/* needed by a boundary point on the left side? */
  if (curve(1,1) != 0 && i_point <= range(1,1) + nw1){
    i_min = range(1,1);
    i_max = range(1,1)+nw1;
    j_min = int_max( j_point - tw2, range(1,2) );
    j_max = int_min( j_point + tw2, range(2,2) );
    for (i=i_min; i<= i_max; i++)
      for (j=j_min; j<= j_max; j++)
	if (abs(flag(i,j)) == grid_ptr->priority)
	  return 1;
  }

/* needed by a boundary point on the right side? */
  if (curve(2,1) != 0 && i_point >= range(2,1) - nw1){
    i_min = range(2,1)-nw1;
    i_max = range(2,1);
    j_min = int_max( j_point - tw2, range(1,2) );
    j_max = int_min( j_point + tw2, range(2,2) );
    for (i=i_min; i<= i_max; i++)
      for (j=j_min; j<= j_max; j++)
	if (abs(flag(i,j)) == grid_ptr->priority)
	  return 1;
  }

/* needed by a boundary point on the lower side? */
  if (curve(1,2) != 0 && j_point <= range(1,2) + nw1){
    i_min = int_max( i_point - tw2, range(1,1) );
    i_max = int_min( i_point + tw2, range(2,1) );
    j_min = range(1,2);
    j_max = range(1,2) + nw1;
    for (i=i_min; i<= i_max; i++)
      for (j=j_min; j<= j_max; j++)
	if (abs(flag(i,j)) == grid_ptr->priority)
	  return 1;
  }

/* needed by a boundary point on the upper side? */
  if (curve(2,2) != 0 && j_point >= range(2,2) - nw1){
    i_min = int_max( i_point - tw2, range(1,1) );
    i_max = int_min( i_point + tw2, range(2,1) );
    j_min = range(2,2) - nw1;
    j_max = range(2,2);
    for (i=i_min; i<= i_max; i++)
      for (j=j_min; j<= j_max; j++)
	if (abs(flag(i,j)) == grid_ptr->priority)
	  return 1;
  }

/* check if the point is needed by a corner formula */
  cw1 = comp_ptr->corner_width - 1;

/* lower left corner */
  if (curve(1,1) != 0 && curve(1,2) != 0 && 
      i_point <= range(1,1)+cw1 && j_point <= range(1,2)+cw1){
    i_min = range(1,1);
    i_max = range(1,1)+cw1;
    j_min = range(1,2);
    j_max = range(1,2)+cw1;
    for (i=i_min; i<= i_max; i++)
      for (j=j_min; j<= j_max; j++)
	if (abs(flag(i,j)) == grid_ptr->priority)
	  return 1;
  }

/* upper left corner */
  if (curve(1,1) != 0 && curve(2,2) != 0 && 
      i_point <= range(1,1)+cw1 && j_point >= range(2,2)-cw1){
    i_min = range(1,1);
    i_max = range(1,1)+cw1;
    j_min = range(2,2)-cw1;
    j_max = range(2,2);
    for (i=i_min; i<= i_max; i++)
      for (j=j_min; j<= j_max; j++)
	if (abs(flag(i,j)) == grid_ptr->priority)
	  return 1;
  }

/* lower right corner */
  if (curve(2,1) != 0 && curve(1,2) != 0 && 
      i_point >= range(2,1)-cw1 && j_point <= range(1,2)+cw1){
    i_min = range(2,1)-cw1;
    i_max = range(2,1);
    j_min = range(1,2);
    j_max = range(1,2)+cw1;
    for (i=i_min; i<= i_max; i++)
      for (j=j_min; j<= j_max; j++)
	if (abs(flag(i,j)) == grid_ptr->priority)
	  return 1;
  }

/* upper right corner */
  if (curve(2,1) != 0 && curve(2,2) != 0 && 
      i_point >= range(2,1)-cw1 && j_point >= range(2,2)-cw1){
    i_min = range(2,1)-cw1;
    i_max = range(2,1);
    j_min = range(2,2)-cw1;
    j_max = range(2,2);
    for (i=i_min; i<= i_max; i++)
      for (j=j_min; j<= j_max; j++)
	if (abs(flag(i,j)) == grid_ptr->priority)
	  return 1;
  }

/* if you get this far, it means that the point (i_point, j_point) is not
    needed by any discretization point */

  return 0;
}

static int 
needed_by_interp(int i_point, int j_point, int priority,
		 overlapping_grid *comp_ptr){
/* check if it is needed by another interpolation point in another grid */
/* this is only possible if the interpolation type is implicit */ 
  component_grid *grid_ptr;
  interp_point *interp_ptr;
  int iw1;

  iw1 = comp_ptr->interp_width - 1;

  if (comp_ptr->interp_type == 'i'){
    for (grid_ptr = comp_ptr->last; grid_ptr != NULL; 
	 grid_ptr = grid_ptr->prev)
      for (interp_ptr = grid_ptr->last_interp; 
	   interp_ptr != NULL;
	   interp_ptr = interp_ptr->prev)
	if (interp_ptr->grid_loc->priority == priority &&
	    i_point >= interp_ptr->i_loc &&
	    i_point <= interp_ptr->i_loc + iw1 &&
	    j_point >= interp_ptr->j_loc &&
	    j_point <= interp_ptr->j_loc + iw1)
	  return 1;
  }
  return 0;
}


static void 
set_flag(component_grid *grid_ptr, int flag_value, int extra_period){
  int i, i_max, i_min, j, j_min, j_max;

  if (grid_ptr->r_period){
    i_min = range(1,1) - extra_period;
    i_max = range(2,1) + extra_period - 1;
  }
  else{
    i_min = range(1,1);
    i_max = range(2,1);
  }

  if (grid_ptr->s_period){
    j_min = range(1,2) - extra_period;
    j_max = range(2,2) + extra_period - 1;
  }
  else{
    j_min = range(1,2);
    j_max = range(2,2);
  }

/* first initialize the whole array including ghostpoints to zero */
  for (i = 1; i<= grid_ptr->flag_ptr->n1; i++)
    for (j = 1; j<= grid_ptr->flag_ptr->n2; j++)
      flag(i,j) = 0;
  
/* now initialize the real points to flag_value */
  for (i = i_min; i <= i_max; i++)
    for (j = j_min; j <= j_max; j++)
      flag(i,j) = flag_value;
    
  
} /* end set_flag */


static int
close_to_boundary(real xp, real yp, int curve_label, overlapping_grid *over_ptr){
  linked_list_member *this_boundary;
  boundary_part *boundary_ptr=NULL;
  real dist, t_x, t_y, min_dist;
  int i_p, i_p_min;

/* Always further than max_n_dist away from all boundaries when curve_label=0 */
  if (curve_label == 0) return FALSE;

  for (this_boundary = over_ptr->boundary_list->first; this_boundary != NULL;
       this_boundary = this_boundary->next){
    boundary_ptr = this_boundary->data;
    if (boundary_ptr->curve_label == curve_label){
      break;
    }
  }

  if (boundary_ptr == NULL){
    printf("ERROR: No part of the global boundary has curve_label=%i\n", curve_label);
    return FALSE;
  }

/* find the closest boundary point */
  min_dist = 1.e10;
  i_p_min = 0;
  for (i_p=1; i_p<boundary_ptr->n_points; i_p++){
/* compute the normal distance */
    dist = sq(xp - boundary_coord(1,i_p)) + sq(yp - boundary_coord(2,i_p));
    if (dist < min_dist){
      i_p_min = i_p;
      min_dist = dist;
    }
  }

/* estimate the tangent */
  i_p = i_p_min;
  if (i_p < boundary_ptr->n_points && 
      (fabs(boundary_coord(1,i_p+1) - boundary_coord(1,i_p)) + 
       fabs(boundary_coord(2,i_p+1) - boundary_coord(2,i_p)) > 0.0)){
    t_x = boundary_coord(1,i_p+1) - boundary_coord(1,i_p);
    t_y = boundary_coord(2,i_p+1) - boundary_coord(2,i_p);
  }
  else if (i_p > 1 && 
	   (fabs(boundary_coord(1,i_p) - boundary_coord(1,i_p-1)) + 
	    fabs(boundary_coord(2,i_p) - boundary_coord(2,i_p-1)) > 0.0)){
    t_x = boundary_coord(1,i_p) - boundary_coord(1,i_p-1);
    t_y = boundary_coord(2,i_p) - boundary_coord(2,i_p-1);
  }
  else{
    printf("Error: close_to_boundary: Failt to estimate the tangent\n");
    return FALSE;
  }

/* compute the normal component of the distance */
  dist = ((xp - boundary_coord(1,i_p)) * t_y - (yp - boundary_coord(2,i_p)) * t_x) / 
    sqrt( t_x*t_x + t_y*t_y );
  if (fabs(dist) <= over_ptr->max_n_dist) 
    return TRUE;
  else
    return FALSE;
}

static int
converted_to_disc_point(int i, int j, component_grid *grid_ptr, 
			overlapping_grid *comp_ptr ){
/* for mixed boundaries, it might be necessary to first change the boundary */
/* point to physical */
  if (i==range(1,1) && left_interp(j)<0 && 
      close_to_boundary(x(i,j), y(i,j), abs(curve(1,1)), comp_ptr)){
    left_interp(j) = abs(left_interp(j));
    if (disc_point(i, j, grid_ptr, comp_ptr)){
      return TRUE;
    }
  }
  if (i==range(2,1) && right_interp(j)<0 && 
      close_to_boundary(x(i,j), y(i,j), abs(curve(2,1)), comp_ptr)){
    right_interp(j) = abs(right_interp(j));
    if (disc_point(i, j, grid_ptr, comp_ptr)){
      return TRUE;
    }
  }
  if (j==range(1,2) && lower_interp(i)<0 && 
      close_to_boundary(x(i,j), y(i,j), abs(curve(1,2)), comp_ptr)){
    lower_interp(i) = abs(lower_interp(i));
    if (disc_point(i, j, grid_ptr, comp_ptr)){
      return TRUE;
    }
  }
  if (j==range(2,2) && upper_interp(i)<0 && 
      close_to_boundary(x(i,j), y(i,j), abs(curve(2,2)), comp_ptr)){
    upper_interp(i) = abs(upper_interp(i));
    if (disc_point(i, j, grid_ptr, comp_ptr)){
      return TRUE;
    }
  }
  return FALSE;
}

static int 
explicit(input_output *io_ptr, overlapping_grid *comp_ptr){
  component_grid *grid_ptr, *donor_ptr;
  int i, j, iw1, implicit, points_removed, convert, interp_lowered, present_curve
    , did_lower, new_flag, convert_ip, converted_location;
  interp_point *interp_ptr, *donor_interp;
  interp_point new_interp;
#ifdef STOP_DURING_EXPLICIT
  int icom, col_flag;
  FILE *fp;
  char * break_com[3], *file_name;
  const int LEVEL=1;
  break_com[0] = "proceed";
  break_com[1] = "break";
  break_com[2] = "postscript-copy";
#endif

#define donor_flag(i,j)  compute_index_2d(donor_ptr->flag_ptr,i,j)

  iw1 = comp_ptr->interp_width - 1;

  for (grid_ptr = comp_ptr->first; grid_ptr != NULL; 
       grid_ptr = grid_ptr->next){ 
    printf("\nExplicit> Grid `%s'.\n",grid_ptr->grid_name);
    points_removed = 0;
    convert = 0;
    convert_ip = 0;
    interp_lowered = 0;
    for (interp_ptr = grid_ptr->last_interp; interp_ptr != NULL; 
	 interp_ptr = interp_ptr->prev){
/* only consider active interpolation points */
      if (interp_ptr->active){
/* make sure all the points in the interpolation formula are 
	   either discretization points or interpolation points that don't
	   interpolate from the present grid */
	do {
	  donor_ptr = interp_ptr->grid_loc;
	  did_lower = 0;
	  implicit = 0;
	  converted_location = FALSE;
	  for (i=interp_ptr->i_loc; i<=interp_ptr->i_loc + iw1; i++)
	    for (j=interp_ptr->j_loc; j<=interp_ptr->j_loc + iw1; j++)
	      if (abs(donor_flag(i,j)) != donor_ptr->priority){
/* first try to convert (i,j) in the donor grid to an interior point */
/* set the curve label if the present point is on a real boundary */
		if (!implicit &&
		    (disc_point( i, j, donor_ptr, comp_ptr) ||
		     converted_to_disc_point( i, j, donor_ptr, comp_ptr)) ){
/* find the interpolation point in the donor grids interpolation list and
		     deactivate it */
		  convert++; converted_location = TRUE;
		  for (donor_interp = donor_ptr->last_interp;
		       donor_interp != NULL; donor_interp = donor_interp->prev){
		    if (donor_interp->i_point == i && donor_interp->j_point == j){
		      donor_interp->active = 0;
		      donor_ptr->n_interp--;
/* also change the point's flag. */
		      donor_flag(i,j) = donor_ptr->priority;
		    }
		  }
		}
		else
/* this implicit interpolation was NOT possible to repair by changing the interpolation
   location to a discretization point */
		  implicit = 1;
	      }
/* change the point to discretization, lower interpolate or unused 
   if implicit interpolation was encountered */
	  if (implicit){
	    i = interp_ptr->i_point; j = interp_ptr->j_point;

/* no initial guess */
	    new_interp.r_loc = 0.0; new_interp.s_loc = 0.0;
	    new_interp.grid_loc = NULL;

/* set the curve label if the present point is on a real boundary */
	    present_curve = get_curve( i, j, grid_ptr, comp_ptr);

/* it should first be checked if the point can interpolate from a grid which is lower
	       than the present donor but higher than itself */

	    new_flag = interp_from_lower( x(i,j), y(i,j), present_curve, &new_interp,
					 grid_ptr, interp_ptr->grid_loc, comp_ptr);

	    if (new_flag > grid_ptr->priority){
	      flag(i, j) = new_flag;
	    }
/* if the new donor priority is lower than the present grid, try to make it 
	       into a discretization point instead */
	    else if (disc_point(i, j, grid_ptr, comp_ptr) ||
		     converted_to_disc_point( i, j, donor_ptr, comp_ptr)){
	      convert_ip++;
	      flag(i, j) = grid_ptr->priority;
	    }
/* interpolate from a grid with lower priority */
	    else{
	      flag(i, j) = new_flag;
	    }


	    if (flag(i, j) == 0 || abs(flag(i, j)) == grid_ptr->priority){
/* deactivate this interpolation point */
	      interp_ptr->active = 0;
	      grid_ptr->n_interp--;
	      points_removed++;
	    }
	    else{
/* update the interpolation */
	      interp_ptr->grid_loc = new_interp.grid_loc;
	      interp_ptr->r_loc = new_interp.r_loc;
	      interp_ptr->s_loc = new_interp.s_loc;
	      interp_ptr->i_loc = new_interp.i_loc;
	      interp_ptr->j_loc = new_interp.j_loc;
	      interp_lowered++;
	      did_lower = 1;
/* 	      printf("The implicit interpolation point (%i, %i) ",  */
/* 		     interp_ptr->i_point, interp_ptr->j_point); */
/* 	      printf("now interpolates from grid `%s'\n" */
/* 		     , interp_ptr->grid_loc->grid_name); */
	    }
	  }
#ifdef STOP_DURING_EXPLICIT
	  if (implicit || converted_location){
	    PL_erase(0);
	    PL_start_plot(0);
	    plot_overlapping( comp_ptr, 1+2 ); 
	    PL_end_plot();
	    printf("Doing interpolation point (%i, %i)\n", interp_ptr->i_point, 
		   interp_ptr->j_point);
	    printf("%s\n", (implicit? "It was implicit": 
			    "The donors were reclassified"));
/* break or proceed? */    
	    while ((icom=get_command(io_ptr, "During explicit>", 
				     break_com, 2, LEVEL, NULL, NULL)) != 0)
	      {
		switch (icom)
		  {
		  case 1:			/* exit to calling routine */
		    return FALSE;
		    break;
		  case 2:			/* postscript copy */
		    if (( fp = open_ascii_file( io_ptr, "postscript file name: ",
					       "classify.ps", &file_name, 'w', 0, 
					       FALSE) ) != NULL)
		      {
			col_flag = get_yes_no( io_ptr, "Do you want to save color "
					      "information>", LEVEL );
			PL_postscript( 0, fp, col_flag );
			fclose( fp );
		      }
		    break;
		  default:
		    ;
		  } /* end switch */
	      } /* end while */
          }
#endif

	} while (did_lower);
      }
    }
/* print a report of the changes on the present grid */
    if (convert > 0)
      printf("Converted %i implicit donor points to discretization points\n", 
	     convert);
    if (interp_lowered > 0)
      printf("Lowered the interpolation location for %i implicit interpolation "
	     "points\n", interp_lowered);
    if (convert_ip > 0)
      printf("Converted %i implicit interpolation points to discretization points\n", 
	     convert_ip);
    if (points_removed-convert_ip > 0)
      printf("Converted %i implicit interpolation points to exterior points\n", 
	     points_removed-convert_ip);
  }

/* remove all deactivated interpolation points */
  for (grid_ptr = comp_ptr->last; grid_ptr != NULL; 
       grid_ptr = grid_ptr->prev){
    remove_deactive_interp( grid_ptr );
  } 

  return 0;

#undef donor_flag
}

static void 
remove_deactive_interp(component_grid *grid_ptr){
/* are there any interpolation points? */
  interp_point *interp_ptr, *prev_victim, *next_victim;

  if (grid_ptr->last_interp != NULL){
    interp_ptr = grid_ptr->last_interp;
    prev_victim = NULL;
    do{
      next_victim = interp_ptr->prev;
      if (!interp_ptr->active){
/* remove this interpolation point */
	if (prev_victim == NULL)
	  grid_ptr->last_interp = interp_ptr->prev;
	else
	  prev_victim->prev = interp_ptr->prev;
	free(interp_ptr);
      }
      else
	prev_victim = interp_ptr;

      interp_ptr = next_victim;
    }
    while (interp_ptr != NULL);
  } 
}

static interp_point *
create_interp_point(interp_point *template_ptr, 
		    interp_point *prev_ptr, 
		    component_grid *grid_ptr){
  interp_point *new_ptr;
  
  if ((new_ptr = (interp_point *)malloc(sizeof(interp_point)) ) == NULL)
    printf("memory error in create_interp_point\n");
  
  new_ptr->i_point = template_ptr->i_point;
  new_ptr->j_point = template_ptr->j_point;
  new_ptr->grid_loc = template_ptr->grid_loc;
  new_ptr->r_loc = template_ptr->r_loc;
  new_ptr->s_loc = template_ptr->s_loc;
  new_ptr->i_loc = template_ptr->i_loc;
  new_ptr->j_loc = template_ptr->j_loc;
  new_ptr->active = 1;
  
/* setup the pointer to the previos interpolation point */
  new_ptr->prev = prev_ptr;
  
/* count the new interpolation point */
  grid_ptr->n_interp++;
  
  return new_ptr;
}


static void 
insert_bad_point(int i, int j, int type, component_grid *grid_ptr){
  bad_point *new_bad_point;

  if ((new_bad_point = (bad_point *)malloc( sizeof(bad_point)) ) == NULL)
    printf("memory error in insert_bad_point\n");
  new_bad_point->i = i;
  new_bad_point->j = j;
  new_bad_point->type = type;

  if (grid_ptr->last_bad_point == NULL)
    new_bad_point->prev = NULL;
  else
    new_bad_point->prev = grid_ptr->last_bad_point;
    
  grid_ptr->last_bad_point = new_bad_point;
}


static int 
boundary_interp( real xp, real yp, component_grid *this_grid, 
		overlapping_grid *comp_ptr){
  component_grid *grid_ptr;
  grid_point gp;
  int width;
  
/*   if (comp_ptr->interp_type == 'e') */
/*     width = comp_ptr->interp_width/2 + (comp_ptr->disc_width - 1)/2; */
/*   else */
  width = comp_ptr->interp_width/2;
  
  for (grid_ptr = comp_ptr->first; 
       grid_ptr != NULL;
       grid_ptr = grid_ptr->next)
    if (grid_ptr != this_grid){
/* call the external grid function to invert the mapping */
      gp.x = xp; gp.y = yp;
      if (invert_mapping(&gp, 0, 0, 0, grid_ptr, comp_ptr)){

/* we first test the interior case: (xp, yp) must be sufficiently far inside 
   the other grid to be an interpolation point. If the point is close to a 
   boundary, it must also be sufficiently far away from the endpoints of the 
   boundary. */
 	if ( (curve(1,1) > 0 || gp.r >= width * grid_ptr->r_step) &&  
 	     (curve(2,1) > 0 || gp.r <= 1.0 - width * grid_ptr->r_step) && 
 	     (curve(1,2) > 0 || gp.s >= width * grid_ptr->s_step) &&  
 	     (curve(2,2) > 0 || gp.s <= 1.0 - width * grid_ptr->s_step) ) 
 	  return 1; 
      } /* end if inside grid_ptr */
    } /* end grid_ptr != this_grid */
/* no match, return zero */
  return 0;
}


static int 
invert_mapping(grid_point *query_ptr, int i_cell, int j_cell,
	       int query_curve, component_grid *grid_ptr, 
	       overlapping_grid *comp_ptr){
/* 
     Invert_mapping investigates if a point with Cartesian coordinates
     (query_ptr->x, query_ptr->y) is inside the grid pointed to by
     grid_ptr.  If the point is found to be outside and query_curve !=
     0, every side of the grid with curve = query_curve is a candidate
     for the following further investigations: The closest point on
     the boundary in question is first located. The point
     (query_ptr->x, query_ptr->y) is considered to be inside of the
     grid if the normal distance from the closest boundary point is
     less than comp_ptr->max_n_dist and the tangential distance is
     less than the local grid size along that side of the grid.
     
     If the point is found to be inside, the (r,s) coordinates in
     grid_ptr corresponding to (query_ptr->x, query_ptr->y) is found
     by a Newton iteration. These quantities are returned in the
     structure pointed to by query_ptr.  
*/
  int i_coord, j_coord, ip, jp, curve_match, inside;
  real missfit, n_dist, parameter;
  
  inside = 0;
  curve_match = 0;
  
/* First test if the point (query_ptr->x, query_ptr->y) is outside the grid's
   bounding box extended by the missfit allowance */
  missfit = (query_curve == 0)? 0.0 : comp_ptr->max_n_dist;
  if (grid_ptr->x_min - missfit > query_ptr->x ||
      grid_ptr->x_max + missfit < query_ptr->x ||
      grid_ptr->y_min - missfit > query_ptr->y ||
      grid_ptr->y_max + missfit < query_ptr->y)
    return 0;

/* Secondly, test if the point (query_ptr->x, query_ptr->y) is inside the grid 
   according to the mapping */
  
/* check if there is an inverse mapping function */
  if( grid_ptr->inverse_known == 1 && grid_ptr->inverse_grid_mapping != NULL){

/* call the inverse mapping function */
    backward_grid_mapping(query_ptr, grid_ptr);

/* was the point inside? */
    if (query_ptr->r >= 0.0 && query_ptr->r <= 1.0 && 
	query_ptr->s >= 0.0 && query_ptr->s <= 1.0){

      inside = 1;
    }
  }
/* When there is no inverse mapping function, we use the ray method to 
     find out if (query_ptr->x, query_ptr->y) is inside */
  else if (locate_cell(query_ptr->x, query_ptr->y, &i_coord, &j_coord, 
		       i_cell, j_cell, grid_ptr)){

    inside = 1;

/* find the (r,s) corresponding to (x,y). Take the (r,s) coordinates in the 
       middle of the cell as initial guess */

    if (grid_ptr->grid_mapping != NULL)
      newton( grid_ptr, query_ptr, 
	     (i_coord - range(1,1) + 0.5) * grid_ptr->r_step,
	     (j_coord - range(1,2) + 0.5) * grid_ptr->s_step );
    else
      local_newton(grid_ptr, query_ptr, 
		   comp_ptr->interp_width, i_coord, j_coord,
		   (i_coord - range(1,1) + 0.5) * grid_ptr->r_step,
		   (j_coord - range(1,2) + 0.5) * grid_ptr->s_step);
  }
  
  
/* Check for curve matches if the point wasn't inside according to the inverse 
   mapping or by the polygonal approximation of the grid boundary. */
  
  if (!inside && query_curve != 0){

/* left curve */
    if (query_curve == abs(curve(1,1))){
      ip = range(1,1);
/* find the closest point */
      closest_boundary_point( grid_ptr, query_ptr->x, query_ptr->y, 1, 
			     &ip, &jp, &n_dist, &parameter );
/* check if it is close enough the be considered as a boundary interpolation point */
      if (fabs(n_dist) <= comp_ptr->max_n_dist && 
	  parameter >= 0.0 && parameter <= 1.0 )
	curve_match = 1;
    }

/* right curve */
    if (!curve_match && query_curve == abs(curve(2,1))){
      ip = range(2,1);
      closest_boundary_point( grid_ptr, query_ptr->x, query_ptr->y, 1, 
			     &ip, &jp, &n_dist, &parameter );
      if (fabs(n_dist) <= comp_ptr->max_n_dist &&
	  parameter >= 0.0 && parameter <= 1.0 )
	curve_match = 1;
    }

/* lower curve */
    if (!curve_match && query_curve == abs(curve(1,2))){
      jp = range(1,2);
      closest_boundary_point( grid_ptr, query_ptr->x, query_ptr->y, 2, 
			     &ip, &jp, &n_dist, &parameter );
      if (fabs(n_dist) <= comp_ptr->max_n_dist &&
	  parameter >= 0.0 && parameter <= 1.0 )
	curve_match = 1;
    }

/* upper curve */
    if (!curve_match && query_curve == abs(curve(2,2))){
      jp = range(2,2);
      closest_boundary_point( grid_ptr, query_ptr->x, query_ptr->y, 2, 
			     &ip, &jp, &n_dist, &parameter );
      if (fabs(n_dist) <= comp_ptr->max_n_dist &&
	  parameter >= 0.0 && parameter <= 1.0 )
	curve_match = 1;
    }

/* find the (r,s) coordinates */
    if (curve_match){
      inside = 1;

/* use the inverse mapping function if there is one */
      if( grid_ptr->inverse_grid_mapping != NULL){
	backward_grid_mapping(query_ptr, grid_ptr);
      }
      else{
/* Invert the forward mapping funtion with Newton's method. 
   Use the (r,s) coordinates of the closest boundary point as
   initial guess */
	if (grid_ptr->grid_mapping != NULL)
	  newton( grid_ptr, query_ptr, 
		 (ip-range(1,1)) * grid_ptr->r_step,
		 (jp-range(1,2)) * grid_ptr->s_step );
	else
	  local_newton(grid_ptr, query_ptr, 
		       comp_ptr->interp_width, ip, jp,
		       (ip - range(1,1)) * grid_ptr->r_step,
		       (jp - range(1,2)) * grid_ptr->s_step);
      }
    }
  }
  
  return inside;
  
}


static void 
closest_boundary_point( component_grid *grid_ptr, real x0, real y0, 
		       int job, int *ip, int *jp, 
		       real *n_dist, real *parameter ){
  real min_dist, dist, x_dist, y_dist, x_step, y_step, local_step_2, increment, d_par;
  int i,j;
  
/* the i-index is fixed, i.e. this affects the left or right boundary */
  if (job == 1){
    j = range(1,2);

    min_dist = (x0 - x(*ip,j)) * (x0 - x(*ip,j)) + 
      (y0 - y(*ip,j)) * (y0 - y(*ip,j));
    *jp = j;
    x_dist = x0 - x(*ip,j);
    y_dist = y0 - y(*ip,j);
    x_step = x(*ip,j+1) - x(*ip,j);
    y_step = y(*ip,j+1) - y(*ip,j);
    d_par = 1.0;

    for (j=range(1,2)+1; j<=range(2,2)-1; j++){
      dist = (x0 - x(*ip,j)) * (x0 - x(*ip,j)) + 
	(y0 - y(*ip,j)) * (y0 - y(*ip,j));
      if (dist < min_dist){
	min_dist = dist;
	*jp = j;
	x_dist = x0 - x(*ip,j);
	y_dist = y0 - y(*ip,j);
	x_step = x(*ip,j+1) - x(*ip,j-1);
	y_step = y(*ip,j+1) - y(*ip,j-1);
	d_par = 2.0;
      }
    }

    j = range(2,2);
    dist = (x0 - x(*ip,j)) * (x0 - x(*ip,j)) + 
      (y0 - y(*ip,j)) * (y0 - y(*ip,j));
    if (dist < min_dist){
      min_dist = dist;
      *jp = j;
      x_dist = x0 - x(*ip,j);
      y_dist = y0 - y(*ip,j);
      x_step = x(*ip,j) - x(*ip,j-1);
      y_step = y(*ip,j) - y(*ip,j-1);
      d_par = 1.0;
    }
  }

/* the j-index is fixed, i.e. this affects the lower or upper boundary */
  else if (job == 2){
    i = range(1,1);

    min_dist = (x0 - x(i,*jp))*(x0 - x(i,*jp)) + 
      (y0 - y(i,*jp))*(y0 - y(i,*jp));
    *ip = i;
    x_dist = x0 - x(i,*jp);
    y_dist = y0 - y(i,*jp);
    x_step = x(i+1,*jp) - x(i,*jp);
    y_step = y(i+1,*jp) - y(i,*jp);
    d_par = 1.0;

    for (i=range(1,1)+1; i<=range(2,1)-1; i++){
      dist = (x0 - x(i,*jp))*(x0 - x(i,*jp)) + 
	(y0 - y(i,*jp))*(y0 - y(i,*jp));
      if (dist < min_dist){
	min_dist = dist;
	*ip = i;
	x_dist = x0 - x(i,*jp);
	y_dist = y0 - y(i,*jp);
	x_step = x(i+1,*jp) - x(i-1,*jp);
	y_step = y(i+1,*jp) - y(i-1,*jp);
	d_par = 2.0;
      }
    }

    i = range(2,1);
    dist = (x0 - x(i,*jp))*(x0 - x(i,*jp)) + 
      (y0 - y(i,*jp))*(y0 - y(i,*jp));
    if (dist < min_dist){
      min_dist = dist;
      *ip = i;
      x_dist = x0 - x(i,*jp);
      y_dist = y0 - y(i,*jp);
      x_step = x(i,*jp) - x(i-1,*jp);
      y_step = y(i,*jp) - y(i-1,*jp);
      d_par = 1.0;
    }
  }
  else{
    printf("ERROR in boundary_normal: Unknown job = %i\n", job);
    exit(1);
  }
  
/* compute the normal component of the distance, and the parameter value of the */
/* projected coordinate */
  local_step_2 = x_step*x_step + y_step*y_step;
  *n_dist = (x_dist * y_step - y_dist * x_step) / sqrt( local_step_2 );

  increment = ( (x0 - x(*ip,*jp)) * x_step + (y0 - y(*ip,*jp)) * y_step )/
    local_step_2 * d_par;

  if (job == 1){
    *parameter = (*jp - range(1,2) + increment) * grid_ptr->s_step;
  }
  else if (job == 2){
    *parameter = (*ip - range(1,1) + increment) * grid_ptr->r_step;
  }

/* done */
}



static int 
newton(component_grid *grid_ptr, grid_point *goal_ptr, real r0, real s0){
  int itr = 0;
  real res, dx, dy, jac;
  grid_point guess;

/* copy the initial guess */
  guess.r = r0; guess.s = s0;

/* newton iteration */

  forward_grid_mapping( &guess, grid_ptr );
  res = fabs(guess.x - goal_ptr->x) + fabs(guess.y - goal_ptr->y);

  while( res > NEWTON_EPS && itr < 5){
    dx = goal_ptr->x - guess.x;
    dy = goal_ptr->y - guess.y;
    if ((jac = guess.xr*guess.ys - guess.xs*guess.yr) == 0.0){
      printf(
"Warning, Newton didn't converge in invert_mapping because the Jacobian\n"
"of the forward mapping became singular at a poor intermediate guess.\n"
"Grid `%s', (x, y) = (%e, %e), initial guess (r, s) = (%e, %e)\n", 
	     grid_ptr->grid_name, goal_ptr->x, goal_ptr->y, r0, s0);
      goal_ptr->r = 0.0; goal_ptr->s = 0.0;
      return 0;
    }
    guess.r += ( guess.ys*dx - guess.xs*dy)/jac;
    guess.s += (-guess.yr*dx + guess.xr*dy)/jac;
    itr++;
    if (forward_grid_mapping( &guess, grid_ptr ) != 0){
      printf(
"Warning, Newton didn't converge in invert_mapping because the\n"
"forward mapping could not be evaluated at a poor intermediate\n"
"guess. Grid `%s', (x, y) = (%e, %e), initial guess (r, s) = (%e, %e)\n", 
	     grid_ptr->grid_name, goal_ptr->x, goal_ptr->y, r0, s0);
      goal_ptr->r = 0.0; goal_ptr->s = 0.0;
      return 0;
    }      
    res = fabs(guess.x-goal_ptr->x) + fabs(guess.y-goal_ptr->y);
/*	printf("iteration %i. residual = %e\n",itr,res); */
  }

/* copy the solution */
  goal_ptr->r = guess.r; goal_ptr->s = guess.s;

/* check for convergence */
  if (res > NEWTON_EPS){
    printf("Warning, Newton didn't converge in invert_mapping. ");
    printf("Grid `%s', (x, y) = (%e, %e), initial guess (r, s) = (%e, %e)\n", 
	   grid_ptr->grid_name, goal_ptr->x, goal_ptr->y, r0, s0);
    printf("Residual = %e\n",res);
    return 0;
  }
/* convergent iteration */
  return 1;
}


static int 
local_newton(component_grid *grid_ptr, grid_point *goal_ptr,
	     int interp_width, int i_close, int j_close,
	     real r0, real s0){
  int itr = 0, iw_low, iw_high, iw1, i_loc, j_loc;
  real res, dx, dy;
  grid_point guess;

/* make sure the interpolation width is even */
  if (interp_width % 2 != 0) interp_width++;

  iw1 = interp_width - 1;
  iw_low = (interp_width-1)/2;
  iw_high = interp_width/2;

/* find the interpolation location */
  if (i_close-iw_low < range(1,1))
    i_loc = range(1,1);
  else if (i_close+iw_high > range(2,1))
    i_loc = range(2,1)-iw1;
  else
    i_loc = i_close-iw_low;

  if (j_close-iw_low < range(1,2))
    j_loc = range(1,2);
  else if (j_close+iw_high > range(2,2))
    j_loc = range(2,2)-iw1;
  else
    j_loc = j_close-iw_low;

/* copy the initial guess */
  guess.r = r0; guess.s = s0;

/* newton iteration */

  if (interp_width == 2)
    bi_linear(i_loc, j_loc, &guess, grid_ptr);
  else if (interp_width == 4)
    bi_cubic(i_loc, j_loc, &guess, grid_ptr);
  else{
    printf("Error in local_newton: only interp_width == 2 and 4 are implemented\n");
    exit(1);
  }
  res = fabs(guess.x - goal_ptr->x) + fabs(guess.y - goal_ptr->y);

  while( res > NEWTON_EPS && itr < 5){
    dx = goal_ptr->x - guess.x;
    dy = goal_ptr->y - guess.y;
    guess.r += ( guess.ys*dx - guess.xs*dy)/(guess.xr*guess.ys - guess.xs*guess.yr);
    guess.s += (-guess.yr*dx + guess.xr*dy)/(guess.xr*guess.ys - guess.xs*guess.yr);
    itr++;
    if (interp_width == 2)
      bi_linear(i_loc, j_loc, &guess, grid_ptr);
    else if (interp_width == 4)
      bi_cubic(i_loc, j_loc, &guess, grid_ptr);
    else{
      printf("Error in local_newton: only interp_width == 2 and 4 are implemented\n");
      exit(1);
    }
    res = fabs(guess.x-goal_ptr->x) + fabs(guess.y-goal_ptr->y);
/*	printf("iteration %i. residual = %e\n",itr,res); */
  }

/* copy the solution */
  goal_ptr->r = guess.r; goal_ptr->s = guess.s;

/* check for convergence */
  if (res > NEWTON_EPS){
    printf("Warning, Newton didn't converge in local_newton. ");
    printf("Grid `%s', (x, y) = (%e, %e), initial guess (r, s) = (%e, %e)\n", 
	   grid_ptr->grid_name, goal_ptr->x, goal_ptr->y, r0, s0);
    printf("Residual = %e\n",res);
    return 0;
  }
/* convergent iteration */
  return 1;
}


static int 
disc_point(int i, int j, component_grid *grid_ptr, 
	   overlapping_grid *comp_ptr){
  int disc_pnt, ip, jp, dw2, cw1, tw2, nw1;

/* The point (i,j) is a valid discretization point if sufficiently many of its
   neighbors are either discretization points or interpolation points */

  dw2 = (comp_ptr->disc_width -1) / 2;
  nw1 = comp_ptr->normal_width - 1;
  cw1 = comp_ptr->corner_width - 1;
  tw2 = (comp_ptr->tangent_width -1) / 2;
  if (tw2 > dw2){
    printf("Disc_point> Warning: tw2>dw2 not implemented\n");
  }

/* discretization point until otherwise proven */
  disc_pnt = 1;

/* 
   PERIODIC IN THE R-DIRECTION 
*/
  if (grid_ptr->r_period){
/* shift `i' if necessary */
    if (i < range(1,1))
      i += range(2,1) - range(1,1);
    else if (i > range(2,1)-1)
      i += -(range(2,1) - range(1,1));
/* interior in the s-direction */
    if (j >= range(1,2)+dw2 && j <= range(2,2)-dw2){
      for (ip=i - dw2; ip<= i + dw2; ip++)
	for (jp=j - dw2; jp<= j + dw2; jp++)
	  if (flag(ip,jp) == 0) disc_pnt = 0;
    }
/* lower side */
    else if (j < range(1,2)+dw2){
      if (curve(1,2) == 0) disc_pnt = 0;
      else{
	for (ip=i-tw2; ip<=i+tw2; ip++)
	  for (jp=range(1,2); jp<=range(1,2)+nw1; jp++)
	    if (flag(ip,jp) == 0) disc_pnt = 0;
      }
    }
/* upper side */
    else if (j > range(2,2)-dw2){
      if (curve(2,2) == 0) disc_pnt = 0;
      else{
	for (ip=i-tw2; ip<=i+tw2; ip++)
	  for (jp=range(2,2)-nw1; jp<=range(2,2); jp++)
	    if (flag(ip,jp) == 0) disc_pnt = 0;
      }
    }
/* this should never happen */
    else
      disc_pnt = 0;
  }
/* 
   PERIODIC IN THE S-DIRECTION 
*/
  else if (grid_ptr->s_period){
/* shift `j' if necessary */
    if (j < range(1,2))
      j += range(2,2) - range(1,2);
    else if (j > range(2,2)-1)
      j += -(range(2,2) - range(1,2));
/* interior in the r-direction */
    if (i >= range(1,1)+dw2 && i <= range(2,1)-dw2){
      for (ip=i - dw2; ip<= i + dw2; ip++)
	for (jp=j - dw2; jp<= j + dw2; jp++)
	  if (flag(ip,jp) == 0) disc_pnt = 0;
    }
/* left side */
    else if (i < range(1,1)+dw2){
      if (curve(1,1) == 0) disc_pnt = 0;
      else{
	for (ip=range(1,1); ip<= range(1,1)+nw1; ip++)
	  for (jp=j - tw2; jp<= j + tw2; jp++)
	    if (flag(ip,jp) == 0) disc_pnt = 0;
      }
    }
/* right side */
    else if (i > range(2,1)-dw2){
      if (curve(2,1) == 0) disc_pnt = 0;
      else{
	for (ip=range(2,1)-nw1; ip<= range(2,1); ip++)
	  for (jp=j - tw2; jp<= j + tw2; jp++)
	    if (flag(ip,jp) == 0) disc_pnt = 0;
      }
    }
/* this should never happen */
    else
      disc_pnt = 0;
  }
/*
   NON-PERIODIC CASE 
*/
  else{
/* interior in i-direction points */
    if (i >= range(1,1)+dw2 && i <= range(2,1)-dw2){
/* interior in both directions */
      if (j >= range(1,2)+dw2 && j <= range(2,2)-dw2){
	for (ip=i - dw2; ip<= i + dw2; ip++)
	  for (jp=j - dw2; jp<= j + dw2; jp++)
	    if (flag(ip,jp) == 0) disc_pnt = 0;
      }
/* lower side away from corners */
      else if (j < range(1,2)+dw2){
	if ( curve(1,2) == 0 || 
	    lower_interp(i) < 0
	    )
	  disc_pnt = 0;
	else{
	  for (ip=i-tw2; ip<=i+tw2; ip++)
	    for (jp=range(1,2); jp<=range(1,2)+nw1; jp++)
	      if (flag(ip,jp) == 0) disc_pnt = 0;
	}
      }
/* upper side away from corners */
      else if (j > range(2,2)-dw2){
	if (curve(2,2) == 0 || 
	    upper_interp(i) < 0
	    )
	  disc_pnt = 0;
	else{
	  for (ip=i-tw2; ip<=i+tw2; ip++)
	    for (jp=range(2,2)-nw1; jp<=range(2,2); jp++)
	      if (flag(ip,jp) == 0) disc_pnt = 0;
	}
      }
/* this should never happen */
      else
	disc_pnt = 0;
    }
/* left side */
    else if (i < range(1,1)+dw2){
/* away from corners */
      if (j>= range(1,2)+dw2 && j<= range(2,2)-dw2){
	if (curve(1,1) == 0 || 
	    left_interp(j) < 0
	    )
	  disc_pnt = 0;
	else{
	  for (ip=range(1,1); ip<= range(1,1)+nw1; ip++)
	    for (jp=j - tw2; jp<= j + tw2; jp++)
	      if (flag(ip,jp) == 0) disc_pnt = 0;
	}
      }
/* lower left corner */
      else if (j < range(1,2)+dw2){
	if (curve(1,1) == 0 || 
	    left_interp(j) < 0
	    ||
	    curve(1,2) == 0 || 
	    lower_interp(i) < 0
	    )
	  disc_pnt = 0;
	else{
	  for (ip=range(1,1); ip <= range(1,1)+cw1; ip++)
	    for (jp=range(1,2); jp<= range(1,2)+cw1; jp++)
	      if (flag(ip,jp) == 0) disc_pnt = 0;
	}
      }
/* upper left corner */
      else if (j > range(2,2)-dw2){
	if (curve(1,1) == 0 || 
	    left_interp(j) < 0
	    ||
	    curve(2,2) == 0 || 
	    upper_interp(i) < 0
	    )
	  disc_pnt = 0;
	else{
	  for (ip=range(1,1); ip <= range(1,1)+cw1; ip++)
	    for (jp=range(2,2)-cw1; jp<= range(2,2); jp++)
	      if (flag(ip,jp) == 0) disc_pnt = 0;
	}
      }
/* this should never happen */
      else
	disc_pnt = 0;
    }
/* right side */
    else if (i > range(2,1)-dw2){
/* away from corners */
      if (j>= range(1,2)+dw2 && j<= range(2,2)-dw2){
	if (curve(2,1) == 0 || 
	    right_interp(j) < 0
	    )
	  disc_pnt = 0;
	else{
	  for (ip=range(2,1)-nw1; ip<= range(2,1); ip++)
	    for (jp=j - tw2; jp<= j + tw2; jp++)
	      if (flag(ip,jp) == 0) disc_pnt = 0;
	}
      }
/* lower right corner */
      else if (j < range(1,2)+dw2){
	if (curve(2,1) == 0 || 
	    right_interp(j) < 0
	    ||
	    curve(1,2) == 0 || 
	    lower_interp(i) < 0
	    )
	  disc_pnt = 0;
	else{
	  for (ip=range(2,1)-cw1; ip <= range(2,1); ip++)
	    for (jp=range(1,2); jp<= range(1,2)+cw1; jp++)
	      if (flag(ip,jp) == 0) disc_pnt = 0;
	}
      }
/* upper right corner */
      else if (j > range(2,2)-dw2){
	if (curve(2,1) == 0 || 
	    right_interp(j) < 0
	    ||
	    curve(2,2) == 0 || 
	    upper_interp(i) < 0
	    )
	  disc_pnt = 0;
	else{
	  for (ip=range(2,1)-cw1; ip <= range(2,1); ip++)
	    for (jp=range(2,2)-cw1; jp<= range(2,2); jp++)
	      if (flag(ip,jp) == 0) disc_pnt = 0;
	}
      }
/* this should never happen */
      else
	disc_pnt = 0;
    }
  }/* end non-periodic case */

  return disc_pnt;
}


static int 
get_curve(int i, int j, component_grid *grid_ptr, 
	  overlapping_grid *comp_ptr){
  int present_curve;
  real left_dist, right_dist, lower_dist, upper_dist, dist_2, r, s;

  dist_2 = comp_ptr->max_n_dist * comp_ptr->max_n_dist;

  left_dist = right_dist = lower_dist = upper_dist = 10.0 * dist_2;

  if (curve(1,1) != 0){
/* corresponding parameter value along the s=0 boundary */
    s = (j - range(1,2)) * grid_ptr->s_step;
/* Make sure that this parameter is NOT in the external gap */
    if (!(grid_ptr->gap_low_r[0] <= s && s <= grid_ptr->gap_low_r[1])){
      if (i == range(1,1))
	left_dist = 0.0;
      else
	left_dist = (x(i,j)-x(range(1,1),j)) * (x(i,j)-x(range(1,1),j)) +
	  (y(i,j)-y(range(1,1),j)) * (y(i,j)-y(range(1,1),j));
    }
  } /* end if curve(1,1) != 0 */

  if (curve(2,1) != 0){
/* corresponding parameter value along the s=0 boundary */
    s = (j - range(1,2)) * grid_ptr->s_step;
/* Make sure that this parameter is NOT in the external gap */
    if (!(grid_ptr->gap_high_r[0] <= s && s <= grid_ptr->gap_high_r[1])){
      if (i == range(2,1))
	right_dist = 0.0;
      else
	right_dist = (x(i,j)-x(range(2,1),j)) * (x(i,j)-x(range(2,1),j)) +
	  (y(i,j)-y(range(2,1),j)) * (y(i,j)-y(range(2,1),j));
    }
  } /* end if curve(2,1) != 0 */

  if (curve(1,2) != 0){
/* corresponding parameter value along the s=0 boundary */
    r = (i - range(1,1)) * grid_ptr->r_step;
/* Make sure that this parameter is NOT in the external gap */
    if (!(grid_ptr->gap_low_s[0] <= r && r <= grid_ptr->gap_low_s[1])){
      if (j == range(1,2))
	lower_dist = 0.0;
      else
	lower_dist = (x(i,j)-x(i,range(1,2))) * (x(i,j)-x(i,range(1,2))) +
	  (y(i,j)-y(i,range(1,2))) * (y(i,j)-y(i,range(1,2)));
    }
  } /* end if curve(1,2) != 0 */

  if (curve(2,2) != 0){
/* corresponding parameter value along the s=0 boundary */
    r = (i - range(1,1)) * grid_ptr->r_step;
/* Make sure that this parameter is NOT in the external gap */
    if (!(grid_ptr->gap_high_s[0] <= r && r <= grid_ptr->gap_high_s[1])){
      if (j == range(2,2))
	upper_dist = 0.0;
      else
	upper_dist = (x(i,j)-x(i,range(2,2))) * (x(i,j)-x(i,range(2,2))) +
	  (y(i,j)-y(i,range(2,2))) * (y(i,j)-y(i,range(2,2)));
    }
  } /* end if curve(2,2) != 0 */
    

  if (left_dist <= dist_2)
    if (lower_dist <= dist_2){
      if (curve(1,1) != 0)
	present_curve = curve(1,1);
      else
	present_curve = curve(1,2);
    }
    else if (upper_dist <= dist_2) {
      if (curve(1,1) != 0)
	present_curve = curve(1,1);
      else
	present_curve = curve(2,2);
    }
    else
      present_curve = curve(1,1);

  else if (right_dist <= dist_2)
    if (lower_dist <= dist_2){
      if (curve(2,1) != 0)
	present_curve = curve(2,1);
      else
	present_curve = curve(1,2);
    }
    else if (upper_dist <= dist_2) {
      if (curve(2,1) != 0)
	present_curve = curve(2,1);
      else
	present_curve = curve(2,2);
    }
    else
      present_curve = curve(2,1);

  else if (lower_dist <= dist_2)
    present_curve = curve(1,2);

  else if (upper_dist <= dist_2)
    present_curve = curve(2,2);

  else
    present_curve = 0;

  return abs(present_curve);
}


static int 
interp_from_higher(real x, real y, int curve, interp_point *new_interp_ptr, 
		   component_grid *other_grid_ptr, overlapping_grid *comp_ptr){
  component_grid *grid_ptr;
  grid_point gp;
  int i_cell, j_cell;

  for (grid_ptr = comp_ptr->first; 
       grid_ptr != NULL && grid_ptr != other_grid_ptr; 
       grid_ptr = grid_ptr->next){

/* use the initial guess if there is one */
    if (grid_ptr == new_interp_ptr->grid_loc){
      i_cell = range(1,1) + new_interp_ptr->r_loc * (range(2,1) - range(1,1));
      j_cell = range(1,2) + new_interp_ptr->s_loc * (range(2,2) - range(1,2));
/* make sure the closest cell is inside the grid */
      i_cell = int_min( int_max( i_cell, range(1,1) ), range(2,1)-1 );
      j_cell = int_min( int_max( j_cell, range(1,2) ), range(2,2)-1 );
    } 
    else{
      i_cell = 0; j_cell = 0;
    }

/* call the external grid function to invert the mapping */
    gp.x = x; gp.y = y;
    if (invert_mapping(&gp, i_cell, j_cell, curve, grid_ptr, comp_ptr)){
/* save the priority and the (r,s) coordinates */
      new_interp_ptr->grid_loc = grid_ptr;
      new_interp_ptr->r_loc = gp.r;
      new_interp_ptr->s_loc = gp.s;
/* check if (i,j) is a valid interpolation location*/
      if (good_interp_loc(new_interp_ptr, grid_ptr, comp_ptr))
	return grid_ptr->priority;
    }
  }
/* no match, return zero */
  return 0;
}


static int 
interp_from_lower(real x, real y, int curve, interp_point *new_interp_ptr,
		  component_grid *own_grid_ptr, 
		  component_grid *other_grid_ptr, overlapping_grid *comp_ptr){
  component_grid *grid_ptr;
  grid_point gp;
  int i_cell, j_cell;
  
  for (grid_ptr = other_grid_ptr->next; 
       grid_ptr != NULL;
       grid_ptr = grid_ptr->next)
    if (grid_ptr != own_grid_ptr){

/* use the initial guess if there is one */
      if (grid_ptr == new_interp_ptr->grid_loc){
	i_cell = range(1,1) + new_interp_ptr->r_loc * (range(2,1) - range(1,1));
	j_cell = range(1,2) + new_interp_ptr->s_loc * (range(2,2) - range(1,2));
/* make sure the closest cell is inside the grid */
	i_cell = int_min( int_max( i_cell, range(1,1) ), range(2,1)-1 );
	j_cell = int_min( int_max( j_cell, range(1,2) ), range(2,2)-1 );
      } 
      else{
	i_cell = 0; j_cell = 0;
      }

/* call the external grid function to invert the mapping */
      gp.x = x; gp.y = y;
      if (invert_mapping(&gp, i_cell, j_cell, curve, grid_ptr, comp_ptr)){
/* save the priority and the (r,s) coordinates */
	new_interp_ptr->grid_loc = grid_ptr;
	new_interp_ptr->r_loc = gp.r;
	new_interp_ptr->s_loc = gp.s;
/* check if (i,j) is a valid interpolation location*/
	if (good_interp_loc(new_interp_ptr, grid_ptr, comp_ptr))
	  return grid_ptr->priority;
      }
    }
/* no match, return zero */
  return 0;
}


static int 
good_interp_loc(interp_point *new_interp_ptr, component_grid *grid_ptr,
		overlapping_grid *comp_ptr){
  int iw1, iw_low, iw_high, i_min, i_max, j_min, j_max, ip, jp, i, j;

/* first check that the point is at least 0.5 grid cells away from all interpolation
   boundaries. */
  if (!grid_ptr->r_period){
    if ((new_interp_ptr->r_loc < 0.5*grid_ptr->r_step && curve(1,1) == 0 ) ||
	(new_interp_ptr->r_loc > 1.0 - 0.5*grid_ptr->r_step && curve(2,1) == 0))
      return 0;
  }
  if (!grid_ptr->s_period){
    if ((new_interp_ptr->s_loc < 0.5*grid_ptr->s_step && curve(1,2) == 0) ||
	(new_interp_ptr->s_loc > 1.0 - 0.5*grid_ptr->s_step && curve(2,2) == 0))
      return 0;
  }

  iw1 =  comp_ptr->interp_width - 1;
  iw_low = (comp_ptr->interp_width - 1)/2;
  iw_high = comp_ptr->interp_width/2;

/* even interpolation width */
  if (comp_ptr->interp_width % 2 == 0){

/* get the enclosing cell */
    i = range(1,1) + new_interp_ptr->r_loc * (range(2,1) - range(1,1));
    j = range(1,2) + new_interp_ptr->s_loc * (range(2,2) - range(1,2));

/* make sure the closest cell is inside the grid */
    i = int_min( int_max( i, range(1,1) ), range(2,1)-1 );
    j = int_min( int_max( j, range(1,2) ), range(2,2)-1 );
  }

/* odd interpolation width */
  else{

/* get the closest point from (r, s) */
    i = range(1,1) + new_interp_ptr->r_loc * (range(2,1) - range(1,1)) + 0.5;
    j = range(1,2) + new_interp_ptr->s_loc * (range(2,2) - range(1,2)) + 0.5;

/* make sure the closest point is inside the grid */
    i = int_min( int_max( i, range(1,1) ), range(2,1) );
    j = int_min( int_max( j, range(1,2) ), range(2,2) );
  }

/* periodic in r */
  if (grid_ptr->r_period){
    if (i < range(1,1)){
      i += range(2,1) - range(1,1);
      new_interp_ptr->r_loc += 1.0;
    }
    else if (i > range(2,1)-1){
      i += -(range(2,1) - range(1,1));
      new_interp_ptr->r_loc += -1.0;
    }
    i_min = i - iw_low;
  }

/* non-periodic in r */
  else {

/* check if (i,j) is close to a boundary */
/* i-direction */
    if (i >= range(1,1)+iw_low && i <= range(2,1)-iw_high){
      i_min = i - iw_low; 
    }
    else if (i < range(1,1)+iw_low){
      i_min = range(1,1); 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (curve(1,1) == 0)	return 0;
    }
    else if (i > range(2,1)-iw_high){
      i_min = range(2,1)-iw1; 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (curve(2,1) == 0)	return 0;
    }
/* this can't happen */
    else
      return 0;
  }

/* periodic in s */
  if (grid_ptr->s_period){
    if (j < range(1,2)){
      j += range(2,2) - range(1,2);
      new_interp_ptr->s_loc += 1.0;
    }
    else if (j > range(2,2)-1){
      j += -(range(2,2) - range(1,2));
      new_interp_ptr->s_loc += -1.0;
    }
    j_min = j - iw_low;
  }

/* non-periodic in s */
  else {
    if (j >= range(1,2)+iw_low && j <= range(2,2)-iw_low){
      j_min = j - iw_low; 
    }
    else if (j < range(1,2)+iw_low){
      j_min = range(1,2); 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (curve(1,2) == 0)	return 0;
    }
    else if (j > range(2,2)-iw_high){
      j_min = range(2,2)-iw1; 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (curve(2,2) == 0)	return 0;
    }
/* this can't happen */
    else
      return 0;
  }

/* save the interpolation location */
  new_interp_ptr->i_loc = i_min;
  new_interp_ptr->j_loc = j_min;

  i_max = i_min + iw1;
  j_max = j_min + iw1;

/* determine if it is a good interpolation location. At this point it is
   not possible to assure explicit interpolation. */

  for (ip = i_min; ip <= i_max; ip++)
    for (jp = j_min; jp <= j_max; jp++)
      if (flag(ip,jp) == 0)
	return 0;

/* the interpolation location was good if you got this far */
  return 1;
}


static void 
change_sign(int i_loc, int j_loc, int iw1, component_grid *grid_ptr){
  int i,j;
  for (i=i_loc; i<=i_loc+iw1; i++)
    for (j=j_loc; j<=j_loc+iw1; j++)
      flag(i,j) = - abs(flag(i,j));
}


static int 
locate_cell(real xp, real yp, int *i_coord, int *j_coord, 
	    int i_cell, int j_cell, component_grid *grid_ptr){
  int i_min=0, i_max=0, i_half, j_min=0, j_max=0, j_half, bad_guess, level
    , half_width, inside00;
  const int max_level=3;

  bad_guess = 1;

/* is there an initial guess? */
  if (i_cell > 0 && j_cell > 0){
/* expand the domain around the initial guess until the point is inside, or */
/* the domain gets too large, or it is not possible to expand the domain */
/* further */
    level = 0;
    half_width = 1;
    do{
/* update the coarsening level; max_level=3 gives maximum half_width = 4. */
      level++;
/* make sure the new domain is within the boundaries */
      i_min = int_max( range(1,1), i_cell - half_width );
      i_max = int_min( range(2,1), i_cell + half_width );
      j_min = int_max( range(1,2), j_cell - half_width );
      j_max = int_min( range(2,2), j_cell + half_width );
/* is it inside now? */      
      bad_guess = !inside_sub_region( xp, yp, i_min, i_max, j_min, j_max, grid_ptr);
/* next size is twice as large */
      half_width *= 2;
    } while( bad_guess && level < max_level );
/* tmp */
/*     if (!bad_guess) */
/*       printf("The initial guess helped at level = %i\n", level); */
/*     else */
/*       printf("The initial guess did not help\n"); */
  }


/* is (xp,yp) inside? */
  if( bad_guess && inside_region(xp, yp, grid_ptr) == 0)
    return 0;

/* start with the whole domain */
  if ( bad_guess ){
    i_min = range(1,1); i_max = range(2,1);
    j_min = range(1,2); j_max = range(2,2);
  }

/* divide the domain into 4 pieces and check in which piece (xp,yp) is */
  do{
    i_half = (i_min + i_max)/2;
    j_half = (j_min + j_max)/2;

/*     inside00 = inside_sub_region(xp,yp,i_min,i_half,j_min,j_half,grid_ptr); */
/*     inside10 = inside_sub_region(xp,yp,i_half,i_max,j_min,j_half,grid_ptr); */
/*     inside01 = inside_sub_region(xp,yp,i_min,i_half,j_half,j_max,grid_ptr); */
/*     inside11 = inside_sub_region(xp,yp,i_half,i_max,j_half,j_max,grid_ptr); */

    if ((inside00 = inside_sub_region(xp,yp,i_min,i_half,j_min,j_half,grid_ptr))
	|| inside_sub_region(xp,yp,i_half,i_max,j_min,j_half,grid_ptr))
      j_max = j_half;
    else
      j_min = j_half;

    if (inside00 || inside_sub_region(xp,yp,i_min,i_half,j_half,j_max,grid_ptr))
      i_max = i_half;
    else
      i_min = i_half;

  } while( i_max-i_min > 1 || j_max-j_min > 1);

  *i_coord = i_min;
  *j_coord = j_min;

  return 1;
}


static int 
inside_region(real xp, real yp, component_grid *grid_ptr){
  int intersect=0;

/* a point is outside if it is outside of the boudning box */
  if ( xp < grid_ptr->x_min || xp > grid_ptr->x_max || 
      yp < grid_ptr->y_min || yp > grid_ptr->y_max)
    return 0;

  if ( !grid_ptr->s_period ){
    intersect += ray_intersect( xp, yp, grid_ptr->lower, grid_ptr ) +
      ray_intersect( xp, yp, grid_ptr->upper, grid_ptr );
  }

  if ( !grid_ptr->r_period ){
    intersect += ray_intersect( xp, yp, grid_ptr->left, grid_ptr ) +
      ray_intersect( xp, yp, grid_ptr->right, grid_ptr );
  }  

/* The point (xp,yp) was outside if intersect is even */  
  return ((intersect % 2 == 0)? 0 : 1); 

}


static int 
inside_sub_region(real xp, real yp, int i_min, int i_max, int j_min, 
		  int j_max, component_grid *grid_ptr){
  int i,j,intersect;
  const real eps=1.0e-7;
  real x_star;

  intersect = 0;

/* a point is always outside a flat or empty domain */
  if (i_min == i_max || j_min == j_max)
    return 0;

  j = j_min;
  for( i=i_min; i<i_max; i++)
/* is yp between y(i,j) and y(i+1,j) ? */
    if ( yp <= real_max(y(i,j),y(i+1,j)) && yp > real_min(y(i,j),y(i+1,j)) )
/* is xp to the right of any point? */
      if ( xp > real_min(x(i,j),x(i+1,j)) )
/* is xp to the right of both points? */
	if (xp >= real_max(x(i,j),x(i+1,j)) )
	  intersect++;
	else{
/* determine the x-coordinate of the straight line between 
   point (i,j) and (i+1,j) */
	  if (fabs(y(i+1,j)-y(i,j)) > eps){
	    x_star = (x(i+1,j)*(yp-y(i,j)) + x(i,j)*(y(i+1,j)-yp))/
	      (y(i+1,j)-y(i,j));
	    if (x_star <= xp)
	      intersect++;
	  }
	  else
	    intersect++;
	}
/*    printf("Inside_region after j_min: intersect = %d\n",intersect);*/
  
  j = j_max;
  for( i=i_min; i<i_max; i++)
/* is yp between y(i,j) and y(i+1,j) ? */
    if ( yp <= real_max(y(i,j),y(i+1,j)) && yp > real_min(y(i,j),y(i+1,j)) )
/* is xp to the right of any point? */
      if ( xp > real_min(x(i,j),x(i+1,j)) )
/* is xp to the right of both points? */
	if (xp >= real_max(x(i,j),x(i+1,j)) )
	  intersect++;
	else{
/* determine the x-coordinate of the straight line between 
   point (i,j) and (i+1,j) */
	  if (fabs(y(i+1,j)-y(i,j)) > eps){
	    x_star = (x(i+1,j)*(yp-y(i,j)) + x(i,j)*(y(i+1,j)-yp))/
	      (y(i+1,j)-y(i,j));
	    if (x_star <= xp)
	      intersect++;
	  }
	  else
	    intersect++;
	}
/*    printf("Inside_region after j_max: intersect = %d\n",intersect);*/


  i = i_min;
  for( j=j_min; j<j_max; j++)
/* is yp between y(i,j) and y(i,j+1) ? */
    if ( yp <= real_max(y(i,j),y(i,j+1)) && yp > real_min(y(i,j),y(i,j+1)) )
/* is xp to the right of any point? */
      if ( xp > real_min(x(i,j),x(i,j+1)) )
/* is xp to the right of both points? */
	if (xp >= real_max(x(i,j),x(i,j+1)) )
	  intersect++;
	else{
/* determine the x-coordinate of the straight line between 
   point (i,j) and (i,j+1) */
	  if (fabs(y(i,j+1)-y(i,j)) > eps){
	    x_star = (x(i,j+1)*(yp-y(i,j)) + x(i,j)*(y(i,j+1)-yp))/
	      (y(i,j+1)-y(i,j));
	    if (x_star <= xp)
	      intersect++;
	  }
	  else
	    intersect++;
	}

/*    printf("Inside_region after i_min: intersect = %d\n",intersect);*/
    
  i = i_max;
  for( j=j_min; j<j_max; j++)
/* is yp between y(i,j) and y(i,j+1) ? */
    if ( yp <= real_max(y(i,j),y(i,j+1)) && yp > real_min(y(i,j),y(i,j+1)) )
/* is xp to the right of any point? */
      if ( xp > real_min(x(i,j),x(i,j+1)) )
/* is xp to the right of both points? */
	if (xp >= real_max(x(i,j),x(i,j+1)) )
	  intersect++;
	else{
/* determine the x-coordinate of the straight line between 
   point (i,j) and (i,j+1) */
	  if (fabs(y(i,j+1)-y(i,j)) > eps){
	    x_star = (x(i,j+1)*(yp-y(i,j)) + x(i,j)*(y(i,j+1)-yp))/
	      (y(i,j+1)-y(i,j));
	    if (x_star <= xp)
	      intersect++;
	  }
	  else
	    intersect++;
	}
/*  printf("Inside_region after i_max: intersect = %d\n",intersect);*/

/* The point (xp,yp) was outside if intersect is even */  
  return ((intersect % 2 == 0)? 0 : 1);
}


static int 
ray_intersect( real xp, real yp, boundary_info *info, component_grid *grid_ptr ){
  int intersect=0;
  int i, j;
  const real eps=1.0e-7;
  real x_star=0.0;

/* there can not be an intersection if the ray is under, over, or starts to the left */
/* of the bounding box */
  if (yp < info->ymin || yp > info->ymax || xp < info->xmin){
    return 0;
  }

/* use the sub-tree if there is one */
   if (info->part1 != NULL && info->part2 != NULL){ 
     intersect = ray_intersect( xp, yp, info->part1, grid_ptr ) +  
       ray_intersect( xp, yp, info->part2, grid_ptr ); 
   } 

/* Count the number of intersections if there is no subtree */
  else{

    if( info->i_const == 0 ){
      j = info->j_const;
      for( i=info->start; i<info->stop; i++)
/* is yp between y(i,j) and y(i+1,j) ? */
	if ( yp <= real_max(y(i,j),y(i+1,j)) && yp > real_min(y(i,j),y(i+1,j)) )
/* is xp to the right of any point? */
	  if ( xp > real_min(x(i,j),x(i+1,j)) )
/* is xp to the right of both points? */
	    if (xp >= real_max(x(i,j),x(i+1,j)) )
	      intersect++;
	    else{
/* determine the x-coordinate of the straight line between 
   point (i,j) and (i+1,j) */
	      if (fabs(y(i+1,j)-y(i,j)) > eps){
		x_star = (x(i+1,j)*(yp-y(i,j)) + x(i,j)*(y(i+1,j)-yp))/
		  (y(i+1,j)-y(i,j));
		if (x_star <= xp)
		  intersect++;
	      }
	      else
		intersect++;
	    }
    }
    else if( info->j_const == 0 ){
      i = info->i_const;
      for( j=info->start; j<info->stop; j++)
/* is yp between y(i,j) and y(i,j+1) ? */
	if ( yp <= real_max(y(i,j),y(i,j+1)) && yp > real_min(y(i,j),y(i,j+1)) )
/* is xp to the right of any point? */
	  if ( xp > real_min(x(i,j),x(i,j+1)) )
/* is xp to the right of both points? */
	    if (xp >= real_max(x(i,j),x(i,j+1)) )
	      intersect++;
	    else{
/* determine the x-coordinate of the straight line between 
   point (i,j) and (i,j+1) */
	      if (fabs(y(i,j+1)-y(i,j)) > eps){
		x_star = (x(i,j+1)*(yp-y(i,j)) + x(i,j)*(y(i,j+1)-yp))/
		  (y(i,j+1)-y(i,j));
		if (x_star <= xp)
		  intersect++;
	      }
	      else
		intersect++;
	    }
    }
    else{
      printf("ERROR in ray_intersect: i_const = %i, j_const = %i\n", 
	     info->i_const, info->j_const );
      exit(1);
    }
  }

/* done */
  return intersect;
}

   
