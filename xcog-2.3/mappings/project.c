
  char *prompt;
  int icom, quit=0, side = 0, dir = 0, n1, n2, n3, i, j, k;
  discrete_surface *sgrid_ptr;
  real_array_2d *x_gap_, *y_gap_, *z_gap_, *x_project, *y_project, *z_project;
  hyp_grid_data *info;
  inverse_point *interp_ptr;
  grid_point gp;
  real s, w, r_loc, s_loc, n_loc, l_scale;
  const real alpha=18.42; /* a pretty universial number */

#define x_gap(i,j) compute_index_2d(x_gap_,i,j)
#define y_gap(i,j) compute_index_2d(y_gap_,i,j)
#define z_gap(i,j) compute_index_2d(z_gap_,i,j)

#include "compute_projection_com.h"

  prompt = "Project face>";

  do{

    switch( get_command( io_ptr, prompt, COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL, SAVE_ON_COPY, ARGUMENT)) 
	     == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case LEFT:
      side = 1; dir = 1;
      quit = 1;
      break;

    case RIGHT:
      side = 2; dir = 1;
      quit = 1;
      break;

    case LOWER:
      side = 1; dir = 2;
      quit = 1;
      break;

    case UPPER:
      side = 2; dir = 2;
      quit = 1;
      break;

/*     case NEAR: */
/*       printf("Sorry, projection of the near (r3=0) face is not implemented yet\n"); */
/*       side = 1; dir = 3; */
/*       break; */

/*     case MY_FAR: */
/*       printf("Sorry, projection of the upper (r3=1) face is not implemented yet\n"); */
/*       side = 2; dir = 3; */
/*       break; */

    case CANCEL:
      return;
      break;

    default:
      break;

    }
  } while (!quit);

/* do the projection for the selected face */

/* construct a discrete surface to approximately invert the surface */
/* mapping `project' */
  surface_coordinates(project, &x_project, &y_project, &z_project);
  sgrid_ptr = new_discrete_surface(project->name, x_project, y_project, z_project, 0, NULL);

/* make the tolerance independent of the mesh size in the `project' patch */
  l_scale = sqrt(sqr(sgrid_ptr->bb->x_max - sgrid_ptr->bb->x_min) +
		 sqr(sgrid_ptr->bb->y_max - sgrid_ptr->bb->y_min) +
		 sqr(sgrid_ptr->bb->z_max - sgrid_ptr->bb->z_min));
  sgrid_ptr->max_n_dist = 0.1 * l_scale;
/* tmp */
/*  printf("Mismatch tolerance for projection surface: %e\n", sgrid_ptr->max_n_dist);*/
/* end tmp */

/* assumes that the mapping is of type hyperbolic */
  info = (hyp_grid_data *) volume->volume_data_ptr;

/* the correspondence between the directions in the mother grid and the component */
/* surface grid are as follows */
/* dir  r_dim   s_dim  */
/*  1   r2_dim  r3_dim */
/*  2   r3_dim  r1_dim */
/*  3   r1_dim  r2_dim */

/* left or right face */
  if (dir == 1){
    n2 = info->n2;
    n3 = info->n3;
/* pick side */
    if (side == 1)
      i = 1;
    else
      i = info->n1;

/* allocate space for the projected surface */
    x_gap_ = create_real_array_2d( n2, n3 );
    y_gap_ = create_real_array_2d( n2, n3 );
    z_gap_ = create_real_array_2d( n2, n3 );

/* do the projection */
    printf("Projecting...\n");
    for (j=1; j<=n2; j++){
      k = 1;
      if ((interp_ptr = search_quad_tree(x3_vol(i,j,k), y3_vol(i,j,k), z3_vol(i,j,k), 
					 sgrid_ptr->quad_tree, sgrid_ptr, 0.0, 0, 0))){
/* save the solution to use as initial guess for next iteration */
	r_loc = interp_ptr->r_loc;
	s_loc = interp_ptr->s_loc;
	n_loc = interp_ptr->n_loc;
/* free the space taken by the interpolation point */
	free(interp_ptr);

/* evaluate the discrete surface */
	gp.r = r_loc; gp.s = s_loc;
	bi_linear(&gp, sgrid_ptr);
/* save the gap */
	x_gap(j,k) = gp.x - x3_vol(i,j,k);
	y_gap(j,k) = gp.y - y3_vol(i,j,k);
	z_gap(j,k) = gp.z - z3_vol(i,j,k);

	for (k=2; k<=n3; k++){
	  newton_search(x3_vol(i,j,k), y3_vol(i,j,k), z3_vol(i,j,k), 
			&r_loc, &s_loc, &n_loc, sgrid_ptr);
/* evaluate the discrete surface */
	  gp.r = r_loc; gp.s = s_loc;
	  bi_linear(&gp, sgrid_ptr);
/* save the gap */
	  x_gap(j,k) = gp.x - x3_vol(i,j,k);
	  y_gap(j,k) = gp.y - y3_vol(i,j,k);
	  z_gap(j,k) = gp.z - z3_vol(i,j,k);
	}
      }
      else{
	printf("compute_projection: Warning: unable to project the point (%e, %e, %e) "
	       "onto the surface mapping `%s'\n", x3_vol(i,j,k), y3_vol(i,j,k), 
	       z3_vol(i,j,k), project->name);
	x_gap(j,k) = 0.0;
	y_gap(j,k) = 0.0;
	z_gap(j,k) = 0.0;
      }
    } /* end for i */

/* update the grid coordinates */
    for (i=1; i<=info->n1; i++){
/* scaled parameter */
      s = (i-1)/((real) info->n1 - 1);
/* use an exponential weighting function */
      if (side == 1){
/* the weighting function is 1 at s=0 and 0 at s=1 */
	w = (exp(-alpha*s) - exp(-alpha))/(1.0-exp(-alpha));
      }
      else{
/* the weighting function is 0 at s=0 and 1 at s=1 */
	w = (exp(-alpha*(1.0-s)) - exp(-alpha))/(1.0-exp(-alpha));
      }

/* Move the k=1 surface eventhough it defines the shape of the surface. In this way */
/* the grid can stay smooth all the way up to the boundary when it is stretched in  */
/* the k-direction */

      for (k=1; k<=info->n3; k++){
	for (j=1; j<=info->n2; j++){
	  x3_vol(i,j,k) = x3_vol(i,j,k) + w*x_gap(j,k);
	  y3_vol(i,j,k) = y3_vol(i,j,k) + w*y_gap(j,k);
	  z3_vol(i,j,k) = z3_vol(i,j,k) + w*z_gap(j,k);
	}
      }
    } /* end for i */

/* cleanup */
    delete_real_array_2d( x_gap_ );
    delete_real_array_2d( y_gap_ );
    delete_real_array_2d( z_gap_ );

  } /* end left or right face */
/* lower or upper face */
  else if (dir == 2){
    n1 = info->n1;
    n3 = info->n3;
/* pick side */
    if (side == 1)
      j = 1;
    else
      j = info->n2;

/* allocate space for the projected surface */
    x_gap_ = create_real_array_2d( n1, n3 );
    y_gap_ = create_real_array_2d( n1, n3 );
    z_gap_ = create_real_array_2d( n1, n3 );

/* do the projection */
    printf("Projecting...\n");
    for (i=1; i<=n1; i++){
      k = 1;
      if ((interp_ptr = search_quad_tree(x3_vol(i,j,k), y3_vol(i,j,k), z3_vol(i,j,k), 
					 sgrid_ptr->quad_tree, sgrid_ptr, 0.0, 0, 0))){
/* save the solution to use as initial guess for next iteration */
	r_loc = interp_ptr->r_loc;
	s_loc = interp_ptr->s_loc;
	n_loc = interp_ptr->n_loc;
/* free the space taken by the interpolation point */
	free(interp_ptr);

/* evaluate the discrete surface */
	gp.r = r_loc; gp.s = s_loc;
	bi_linear(&gp, sgrid_ptr);
/* save the gap */
	x_gap(i,k) = gp.x - x3_vol(i,j,k);
	y_gap(i,k) = gp.y - y3_vol(i,j,k);
	z_gap(i,k) = gp.z - z3_vol(i,j,k);

	for (k=1; k<=n3; k++){
	  newton_search(x3_vol(i,j,k), y3_vol(i,j,k), z3_vol(i,j,k), 
			&r_loc, &s_loc, &n_loc, sgrid_ptr);
/* evaluate the discrete surface */
	  gp.r = r_loc; gp.s = s_loc;
	  bi_linear(&gp, sgrid_ptr);
/* save the gap */
	  x_gap(i,k) = gp.x - x3_vol(i,j,k);
	  y_gap(i,k) = gp.y - y3_vol(i,j,k);
	  z_gap(i,k) = gp.z - z3_vol(i,j,k);
	}
      }
      else{
	printf("compute_projection: Warning: unable to project the point (%e, %e, %e) "
	       "onto the surface mapping `%s'\n", x3_vol(i,j,k), y3_vol(i,j,k), 
	       z3_vol(i,j,k), project->name);
	x_gap(i,k) = 0.0;
	y_gap(i,k) = 0.0;
	z_gap(i,k) = 0.0;
      }
    } /* end for i */

/* update the grid coordinates */
    for (j=1; j<=info->n2; j++){
/* scaled parameter */
      s = (j-1)/((real) info->n2 - 1);
/* use an exponential weighting function */
      if (side == 1){
/* the weighting function is 1 at s=0 and 0 at s=1 */
	w = (exp(-alpha*s) - exp(-alpha))/(1.0-exp(-alpha));
      }
      else{
/* the weighting function is 0 at s=0 and 1 at s=1 */
	w = (exp(-alpha*(1.0-s)) - exp(-alpha))/(1.0-exp(-alpha));
      }

/* Move the k=1 surface eventhough it defines the shape of the surface. In this way */
/* the grid can stay smooth all the way up to the boundary when it is stretched in  */
/* the k-direction */

      for (k=1; k<=info->n3; k++){
	for (i=1; i<=info->n1; i++){
	  x3_vol(i,j,k) = x3_vol(i,j,k) + w*x_gap(i,k);
	  y3_vol(i,j,k) = y3_vol(i,j,k) + w*y_gap(i,k);
	  z3_vol(i,j,k) = z3_vol(i,j,k) + w*z_gap(i,k);
	}
      }
    }

/* cleanup */
    delete_real_array_2d( x_gap_ );
    delete_real_array_2d( y_gap_ );
    delete_real_array_2d( z_gap_ );

  } /* end lower or upper face */

/* delete the component grid used to invert the project_onto mapping */
  delete_discrete_surface( sgrid_ptr );
