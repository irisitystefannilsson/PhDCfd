#include "mappings.h"

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

typedef struct{
  int n1;
  real_array_1d *x_, *y_;
} grid_function;

typedef struct{
  real lambda_1r, lambda_1i, lambda_2r, lambda_2i;
  grid_function *x1, *x2, *Ax1, *Ax2, *u_eps, *u_eps_t;
} eigen_val;

/* private member functions */
static void
compute_projection( input_output *io_ptr, generic_mapping *volume, generic_curve *project );
static int
grid_normal(real *n_vec, grid_function *u, int i, hyp_grid_data *info);
static void
set_hyp_bc( input_output *io_ptr, hyp_grid_data *info );
static void
power_method_1( grid_function *u, grid_function *u_t, real t, hyp_grid_data *info,
	     eigen_val *eval, real n_sign );
static real
mean_curvature( grid_function *u, real_array_1d *kappa_, int periodic );
static int 
hyp_grid_mapping( grid_point *gp, void *data_ptr);
static grid_function *
new_grid_function(int n1);
static grid_function *
delete_grid_function(grid_function *u);
static void
time_derivative(grid_function *u, real t, grid_function *u_t, hyp_grid_data *info, real n_sign);
static void 
compute_hyp_grid(input_output *io_ptr, generic_mapping *grid_ptr);
static void 
set_hyp_grid(input_output *io_ptr, generic_mapping *vgrid, generic_curve *first_c);
static void *
delete_hyp_grid( void *vgrid_data_ptr );
static void 
write_hyp_map( int32 dir, void *data_ptr );
static void *
copy_hyp_map( void *mapping_data_ptr );
/* mappings for smoothed polygon grids */

static void 
set_hyp_grid(input_output *io_ptr, generic_mapping *grid_ptr, generic_curve *first_c){
  hyp_grid_data *info;
  int quit=0, replot=TRUE, icom;
  const real eps=1.e-4;
  char *bc_names[3];
  generic_curve *project_curve;
#include "set_hyp_com.h"

  bc_names[0] = "free";
  bc_names[1] = "x-constant";
  bc_names[2] = "y-constant";

/* cast the data pointer to the right type */
  info = (hyp_grid_data *) grid_ptr->mapping_data_ptr;
  PL_window(1, grid_ptr->xyab);

  do{

    if (replot){
/* plot the mapping */
      PL_erase(1);
      PL_start_plot(1);
/* if info->status != 1, only draw the curve */
      if (info->status == 1)
	plot_generic_mapping( grid_ptr, grid_ptr->grid_plot_mode );
      else
	plot_curve( info->curve_ptr, 1|2 );
      PL_end_plot();
      replot = FALSE;
    }

    switch( get_command( io_ptr, "hyperbolic mapping>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case BOTH_SIDES:
      info->both_sides = !info->both_sides;
      if (!info->both_sides)
	printf("The grid will only be grown out from ONE side of the curve\n");
      else
	{
	  printf("The grid will be grown out from BOTH sides of the curve\n");
	}
      info->status = 0;
      break;

    case PROJECT_SIDE:
      if (info->status == 1){
	printf("Identify the curve to project onto:\n");
	project_curve = get_curve_ptr( io_ptr, first_c );
/* compute the projection */  
	if (project_curve){
	  compute_projection(io_ptr, grid_ptr, project_curve);
	  replot = TRUE;
	}
      }
      else{
	printf("You must update the hyperbolic grid before you can project "
	       "the grid points. This is done with the command `compute-mapping'.\n");
      }
      break;

    case C_GRID:
/* check if the curve is periodic */
      if (!info->curve_ptr->periodic){
	printf("Sorry, the C-grid mapping only works for periodic curves\n");
	info->c_grid = FALSE;
	break;
      }
      info->c_grid = ~info->c_grid;
      if (info->c_grid){
	grid_ptr->r_period = FALSE;
	info->periodic = FALSE;
	info->n_cut = 6;
	info->x_branch = 1.0;
	info->y_branch = 0.0;
      }
      else{
	info->periodic = grid_ptr->r_period = info->curve_ptr->periodic;
      }
      printf("The C-grid mode is now %s\n", (info->c_grid? "ON":"OFF"));
/* recompute the mapping */
      info->status = 0;
      break;

    case N_CUT:
      if (info->c_grid){
	info->n_cut = int_max(1, 
			      get_int(io_ptr, 
				      "Number of grid points on the branch-cut (>0):", 
				      info->n_cut, NO_INDENT));
/* recompute the mapping */
	info->status = 0;
      }
      else{
	printf("This command is only applicable to c-grids.\n");
      }
      break;

    case BRANCH_LOC:
      if (info->c_grid){
	info->x_branch = get_real(io_ptr, "Relative x-coordinate of start point "
				  "for branch-cut:", info->x_branch, LEVEL+1);
	info->y_branch = get_real(io_ptr, "Relative y-coordinate of start point "
				  "for branch-cut:", info->y_branch, LEVEL+1);
/* recompute the mapping */
	info->status = 0;
      }
      else{
	printf("This command is only applicable to c-grids.\n");
      }
      break;

    case QUALITY:
      if (grid_ptr->grid_plot_mode & 256)
	{
	  grid_ptr->grid_plot_mode &= ~256;
	  printf("Grid quality check toggled OFF\n");
	}
      else
	{
	  grid_ptr->grid_plot_mode |= 256;
	  printf("Grid quality check toggled ON\n");
	}
      replot = TRUE;
      break;

    case COMPUTE:
      if (info->status == 1){
	printf("The mapping is up to date with the parameters, \n"
	       "so there is no need to recompute it!\n");
      }
      else{
	compute_hyp_grid( io_ptr, grid_ptr );
	PL_window(1, grid_ptr->xyab);
	replot = TRUE;
      }
      break;

    case CURV_FACTOR:
      info->curv_factor = real_max( 0.0, get_real(io_ptr, "curvature factor:", 
						  info->curv_factor, NO_INDENT));
      info->status = 0;
      break;

    case SMOOTH_FACTOR:
      info->smooth_factor = 
	real_min( 1.0, real_max( 0.0, get_real(io_ptr, "averaging factor (>=0, <=1):", 
					       info->smooth_factor, NO_INDENT)) );
      info->status = 0;
      break;

    case V_MIN:
      info->v_min = real_max( 0.0, get_real(io_ptr, "smallest normal velocity:", 
					    info->v_min, NO_INDENT));
      info->status = 0;
      break;

    case SHOW:
      printf("This mapping is based on the curve `%s'\n", 
	     info->curve_ptr->curve_name);
      printf("averaging coefficient: %e\n", info->smooth_factor);
      printf("curvature coefficient: %e\n", info->curv_factor);
      printf("velocity threshold: %e\n", info->v_min);
      printf("thickness: %e\n", info->width);
      printf("Boundary condition at r1=0: %s, r1=1: %s.\n",
	     bc_names[hyp_bc(1)], bc_names[hyp_bc(2)]);
      printf("r-points: %i\ns-points: %i\n",
	     grid_ptr->r_points, grid_ptr->s_points);
      printf("\n");
      break;

    case BOUNDARY_CONDITION:
      set_hyp_bc( io_ptr, info );
/* recompute the mapping */
      info->status = 0;
      break;

    case WIDTH:
      info->width = real_max( eps, get_real( io_ptr, "thickness:", 
					    info->width, NO_INDENT));
/* recompute recompute the mapping */
      info->status = 0;
      break;

    case CHANGE_CURVE:
/* change the curve */
      set_curve( io_ptr, info->curve_ptr );
/* recompute the mapping */
      info->status = 0;
      break;

    case FLIP:
/* flip curve parametrization */
      info->flip = !info->flip;
/* recompute the mapping */
      info->status = 0;
      break;

    case R1_LINES:
      grid_ptr->r_points = int_max(2, get_int(io_ptr, "Number of r-lines:", 
					     grid_ptr->r_points, NO_INDENT));
      info->status = 0;
      break;

    case R2_LINES:
      grid_ptr->s_points = int_max(2, get_int(io_ptr, "Number of s-lines:", 
					     grid_ptr->s_points, NO_INDENT));
      info->status = 0;
      break;

    case S_STRETCHING:
/* normal stretching */
      if (info->normal_stretch == NULL){
	info->normal_stretch = choose_stretching(io_ptr, NULL, &grid_ptr->s_points );
	grid_ptr->grid_plot_mode |= 8;
	info->status = 0;
      }
      else{
	if (set_stretching( io_ptr, info->normal_stretch, NULL, &grid_ptr->s_points ) ){
	  info->status = 0;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->normal_stretch ) );
	}
      }
      break;

    case NO_S_STRETCHING:
/* no normal stretching */
      if (info->normal_stretch != NULL){
	info->normal_stretch = delete_generic_stretching( info->normal_stretch );
	grid_ptr->grid_plot_mode |= 8;
	info->status = 0;
      }
      else{
	printf("There is no stretching in the r3-direction.\n");
      }
      break;

    case R_STRETCHING:
/* tangential stretching */
      if (info->tangential_stretch == NULL){
	info->tangential_stretch = choose_stretching(io_ptr, info->curve_ptr, 
						     &grid_ptr->r_points );
	grid_ptr->grid_plot_mode |= 8;
	info->status = 0;
      }
      else{
	if (set_stretching(io_ptr, info->tangential_stretch, info->curve_ptr, 
			   &grid_ptr->r_points ) ){
	  info->status = 0;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->tangential_stretch ) );
	}
      }
      break;

    case NO_R_STRETCHING:
/* no tangential stretching */
      if (info->tangential_stretch != NULL){
	info->tangential_stretch = delete_generic_stretching( info->tangential_stretch );
	grid_ptr->grid_plot_mode |= 8;
	info->status = 0;
      }
      else{
	printf("There is no stretching in the r3-direction.\n");
      }
      break;

/* change the plot mode */
    case PLOT_MODE:
      if (info->status == 1){
	set_grid_plot_mode(io_ptr, grid_ptr); 
/* don't need to redraw the surface grids! */
      }
      else{
	printf("Sorry, you must compute the hyperbolic grid before you can look at it!\n"
	       "This is done with the command `compute-mapping'.\n");
      }
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL, NULL, NULL)) 
	     == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
/* make sure that we exit with a valid grid */
      if (info->status == 1){
	quit = 1;
      }
      else{
	printf(
"Sorry, you can only exit after the mapping has been successfully computed.\n"
"If you have trouble making the grid, I suggest that you make the thickness\n"
"very small, increase the averaging coefficient and decrease the curvature\n"
"coefficient.\n");
      }
      break;

    default:
      break;

    }
  } while (!quit);
}

static void
set_hyp_bc( input_output *io_ptr, hyp_grid_data *info ){
  char *bc_names[3];
  int icom, quit=0;

#include "set_hyp_bc_com.h"

  bc_names[0] = "free";
  bc_names[1] = "x-constant";
  bc_names[2] = "y-constant";

  do{

    switch (get_command( io_ptr, "boundary condition>", COMMAND, LAST_COM, LEVEL+1, 
			SAVE_ON_COPY, ARGUMENT)){

    case LOW_R1:
      hyp_bc(1) = get_command( io_ptr, "Low r1 bc>", bc_names, 2, NO_INDENT, 
			      NULL, NULL);
      break;

    case HIGH_R1: 
      hyp_bc(2) = get_command( io_ptr, "High r1 bc>", bc_names, 2, NO_INDENT, 
			      NULL, NULL);
      break;

    case SHOW:
      printf("Boundary condition at r1=0: %s, r1=1: %s.\n",
	     bc_names[hyp_bc(1)], bc_names[hyp_bc(2)]);
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
      quit = 1;
      break;

    default:
      break;
    }
  }
  while(!quit);
}


static grid_function *
new_grid_function(int n1){
  grid_function *u;
  u = (grid_function *) malloc( sizeof(grid_function) );
  u->n1 = n1;
  u->x_ = create_real_array_1d(n1);
  u->y_ = create_real_array_1d(n1);
  return u;
}

static grid_function *
delete_grid_function(grid_function *u){
  if( u ){
    delete_real_array_1d(u->x_);
    delete_real_array_1d(u->y_);
    free(u);
  }
  return NULL;
}

static void 
compute_hyp_grid(input_output *io_ptr, generic_mapping *grid_ptr){
  hyp_grid_data *info;
  grid_function *u, *uh, *u_t;
  curve_point cp;
  int i, j, k, itr, n2_2;
  real dr, dt, t, t_min, t_max, first_lambda=-1.0; /* just to keep the gcc compiler quiet */
  real r_stretch, dstr_duni, d2str_duni2, n_sign=1.0;
  real x_0, y_0, x_1, y_1, alpha, dr_c;
  const real eps=1.e-10;
  real_array_1d *kappa_;
  eigen_val *eval;

#define u1(i) compute_index_1d(u->x_, i)
#define u2(i) compute_index_1d(u->y_, i)

#define uh1(i) compute_index_1d(uh->x_, i)
#define uh2(i) compute_index_1d(uh->y_, i)

#define u_t1(i) compute_index_1d(u_t->x_, i)
#define u_t2(i) compute_index_1d(u_t->y_, i)

#define c_1(u,i) compute_index_1d(u->x_, i)
#define c_2(u,i) compute_index_1d(u->y_, i)

/* cast the data pointer to the right type */
  info = (hyp_grid_data *) grid_ptr->mapping_data_ptr;

/*   0. free old x, y, z arrays */
  info->x_ = delete_real_array_2d( info->x_ );
  info->y_ = delete_real_array_2d( info->y_ );

/*   1. allocate new x, y, z arrays */
  info->n1 = grid_ptr->r_points;
  info->n2 = grid_ptr->s_points;

  if (info->c_grid){
    info->n_curve = info->n1 - 2*info->n_cut;
  }
  else{
    info->n_curve = info->n1;
    info->n_cut = 0;
  }

  dr = info->d1 = 1.0/(info->n1-1);
  info->d2 = 1.0/(info->n2-1);

  if (!info->both_sides)
    n2_2 = 1;
  else
    {
      n2_2 = (info->n2 - 1)/2 + 1;
    }

  info->x_ = create_real_array_2d( info->n1, info->n2 );
  info->y_ = create_real_array_2d( info->n1, info->n2 );

/*   2. allocate 3 gridfunctions, 1 for the present solution, 1 for the middle step, */
/*      and 1 for the time-derivative. */
  u   = new_grid_function(info->n1);
  uh  = new_grid_function(info->n1);
  u_t = new_grid_function(info->n1);

/*   2.5 assign initial data */
  if (info->c_grid)
    {
/* get the coordinates of the periodic point */
      cp.r = 0.0;
      curve_function( &cp, info->curve_ptr );
      x_0 = cp.x;
      y_0 = cp.y;
/* the starting point of the branch-cut */
      x_1 = x_0 + info->x_branch;
      y_1 = y_0 + info->y_branch;
/* first side of the branch */
      for (i=1; i<=info->n_cut; i++)
	{
	  alpha = (i-1)/((real)info->n_cut);
	  u1(i) = x_hyp(i,n2_2) = x_1 + alpha * (x_0 - x_1);
	  u2(i) = y_hyp(i,n2_2) = y_1 + alpha * (y_0 - y_1);
	}
/* along the curve */
      dr_c = 1.0/(info->n_curve-1);
      for (i=1; i<=info->n_curve; i++)
	{
/* assign the parameters in the curve structure */
	  cp.r = dr_c*(i-1);
/* flip the parameter? */
	  if (info->flip)
	    cp.r = 1.0 - cp.r;

/* stretching in the tangential direction ? */
	  uniform_to_stretch(cp.r, info->tangential_stretch, &r_stretch, 
			     &dstr_duni, &d2str_duni2 );
/* Evaluate the curve */
	  cp.r = r_stretch;
	  curve_function( &cp, info->curve_ptr );
/* copy values */
	  u1(info->n_cut+i) = x_hyp(info->n_cut+i,n2_2) = cp.x;
	  u2(info->n_cut+i) = y_hyp(info->n_cut+i,n2_2) = cp.y;
	}
/* copy the second side of the branch from the first */
      for (i=1; i<=info->n_cut; i++)
	{
	  u1(info->n_cut+info->n_curve+i) = x_hyp(info->n_cut+info->n_curve+i,n2_2) = 
	    u1(info->n_cut-i+1);
	  u2(info->n_cut+info->n_curve+i) = y_hyp(info->n_cut+info->n_curve+i,n2_2) = 
	    u2(info->n_cut-i+1);
	}
    }
  else
    {
      for (i=1; i<=info->n1; i++)
	{
/* assign the parameters in the curve structure */
	  cp.r = dr*(i-1);
/* flip the parameter? */
	  if (info->flip)
	    cp.r = 1.0 - cp.r;

/* stretching in the tangential direction ? */
	  uniform_to_stretch(cp.r, info->tangential_stretch, &r_stretch, 
			     &dstr_duni, &d2str_duni2 );
/* Evaluate the curve */
	  cp.r = r_stretch;
	  curve_function( &cp, info->curve_ptr );
/* copy values */
	  u1(i) = x_hyp(i,n2_2) = cp.x;
	  u2(i) = y_hyp(i,n2_2) = cp.y;
	}
    }

  kappa_ = create_real_array_1d( u->n1 );
  printf("Max kappa from initial data: %e\n", 
	 mean_curvature( u, kappa_, info->periodic ));
  delete_real_array_1d( kappa_ );

/* initialize the eigenvalue structure*/

  dr = 2.0*((double)M_PI)/(u->n1-1);

/* allocate space */
  eval = (eigen_val *) malloc( sizeof(eigen_val) );
  eval->lambda_1r = eval->lambda_1i = eval->lambda_2r = eval->lambda_2i = 1.e10;
  eval->u_eps   = new_grid_function( u->n1 );
  eval->u_eps_t = new_grid_function( u->n1 );
  eval->x1    = new_grid_function( u->n1 );
  eval->x2    = new_grid_function( u->n1 );
  eval->Ax1   = new_grid_function( u->n1 );
  eval->Ax2   = new_grid_function( u->n1 );

/* initial guesses for Ax1 and Ax2 */
  for (i=1; i<=u->n1; i++)
    {
      c_1(eval->Ax1,i) = (1.0 - 2.0*((i)%2)) * cos( (i-1)*dr );
      c_2(eval->Ax1,i) = (1.0 - 2.0*((i)%2)) * cos( (i-1)*dr );
      c_1(eval->Ax2,i) = (1.0 - 2.0*((i)%2)) * sin( (i-1)*dr );
      c_2(eval->Ax2,i) = (1.0 - 2.0*((i)%2)) * sin( (i-1)*dr );
    }

/* done initializing the eigenvalue structure */

  printf("Layer %i is initial data\n", n2_2);
/*   3. for layer = n2_2+1 to layer = info->n2, step 1: */
  n_sign = 1.0;
  for (k=n2_2 + 1; k<=info->n2; k++){
    printf("Working on layer %i\n", k);
/* copy previous layer to the time-integration variables */
    for (i=1; i<= info->n1; i++)
      {
	u1(i) = x_hyp(i,k-1);
	u2(i) = y_hyp(i,k-1);
      }
/*   4. Integrate u from time t_{k-1} to t_{k} with RK-2. */
    t     = (k-2)/((real) info->n2-1);
    t_max = (k-1)/((real) info->n2-1);
/*    printf("Working on layer %i. t0 = %e and tmax = %e.\n", k, t, t_max);*/
    itr = 0;
    do{
/* first stage */
      time_derivative( u, t, u_t, info, n_sign );

/* compute the eigenvalues. We could if necessary do this in every time-step, if */
/* it turns out to be hard to estimate them analytically. */
      power_method_1( u, u_t, t, info, eval, n_sign );
      if (eval->lambda_1r > 0.0)
	printf("Warning: Positive dominating eigenvalue = %e\n", eval->lambda_1r);
      if (itr == 0){
	first_lambda = fabs(eval->lambda_1r);
      }
      else if (fabs(eval->lambda_1r) > 3.0*first_lambda){
	printf(
"Warning: The grid generation was prematurely halted because the time-integration \n"
"seemed to go unstable. Please try again with different parameters. It usually helps \n"
"to decrease the curvature coefficient and increase the averaging coefficient.\n"
"Furthermore, it also helps to make the grid thinner.\n");
/* unsuccessful computation */
	info->status = -1;

/*  delete the gridfunctions the present solution, the middle step, */
/*  and the time-derivative. */
	delete_grid_function(u);
	delete_grid_function(uh);
	delete_grid_function(u_t);

/* free the eigenvalue structure */
	delete_grid_function( eval->u_eps );
	delete_grid_function( eval->u_eps_t );
	delete_grid_function( eval->x1 );
	delete_grid_function( eval->x2 );
	delete_grid_function( eval->Ax1 );
	delete_grid_function( eval->Ax2 );
	free(eval);
	return;
      }
/* estimate largest stable dt for R-K-2 */
      dt = 2.0*info->time_factor / fabs(eval->lambda_1r);

/* do not integrate past t_max */
/*      printf("t_max-t = %e, dt = %e\n", t_max-t, dt);*/
      dt = real_min( t_max-t, dt );
/*      printf("Actual dt = %e\n", dt);*/

/* update */
      for (i=1; i<=info->n1; i++)
	{
	  uh1(i) = u1(i) + 0.5 * dt * u_t1(i);
	  uh2(i) = u2(i) + 0.5 * dt * u_t2(i);
	}
/* apply bc */

/* second stage */
      time_derivative( uh, t+0.5*dt, u_t, info, n_sign );
/* update */
      for (i=1; i<=info->n1; i++)
	{
	  u1(i) = u1(i) + dt * u_t1(i);
	  u2(i) = u2(i) + dt * u_t2(i);
	}
/* apply bc */

/* update time */
      t += dt;
/* increment the iteration counter */
      itr++;
/*      printf("Present time = %e\n", t);*/
    } while (t < t_max-eps);
    

/* assign the new layer from the time-integration variables */
    for (i=1; i<= info->n1; i++)
      {
	x_hyp(i,k) = u1(i);
	y_hyp(i,k) = u2(i);
      }
  } /* end for all k */

/* backwards if both_sided */
  if (info->both_sides)
    {
      n_sign=-1.0;
      printf("Doing the second side...\n");
      printf("Layer %i is initial data\n", n2_2);
/* initial guesses for Ax1 and Ax2 */
      for (i=1; i<=u->n1; i++)
	{
	  c_1(eval->Ax1,i) = (1.0 - 2.0*((i)%2)) * cos( (i-1)*dr );
	  c_2(eval->Ax1,i) = (1.0 - 2.0*((i)%2)) * cos( (i-1)*dr );
	  c_1(eval->Ax2,i) = (1.0 - 2.0*((i)%2)) * sin( (i-1)*dr );
	  c_2(eval->Ax2,i) = (1.0 - 2.0*((i)%2)) * sin( (i-1)*dr );
	}
      eval->lambda_1r = eval->lambda_1i = eval->lambda_2r = eval->lambda_2i = 1.e10;
/* done initializing the eigenvalue structure */

  for (k=n2_2 - 1; k>=1; k--){
    printf("Working on layer %i\n", k);
/* copy previous layer to the time-integration variables */
    for (i=1; i<= info->n1; i++)
      {
	u1(i) = x_hyp(i,k+1);
	u2(i) = y_hyp(i,k+1);
      }
/*   4. Integrate u from time t_{k-1} to t_{k} with RK-2. */
    t     = k/((real) info->n2-1);
    t_min = (k-1)/((real) info->n2-1);
/*    printf("Working on layer %i. t0 = %e and tmax = %e.\n", k, t, t_max);*/
    itr = 0;
    do{
/* first stage */
      time_derivative( u, t, u_t, info, n_sign );

/* compute the eigenvalues. We could if necessary do this in every time-step, if */
/* it turns out to be hard to estimate them analytically. */
      power_method_1( u, u_t, t, info, eval, n_sign );
      if (eval->lambda_1r > 0.0)
	printf("Warning: Positive dominating eigenvalue = %e\n", eval->lambda_1r);
      if (itr == 0)
	{
	first_lambda = fabs(eval->lambda_1r);
      }
      else if (fabs(eval->lambda_1r) > 3.0*first_lambda)
	{
	printf(
"Warning: The grid generation was prematurely halted because the time-integration \n"
"seemed to go unstable. Please try again with different parameters. It usually helps \n"
"to decrease the curvature coefficient and increase the averaging coefficient.\n"
"Furthermore, it also helps to make the grid thinner.\n");
/* unsuccessful computation */
	info->status = -1;

/*  delete the gridfunctions the present solution, the middle step, */
/*  and the time-derivative. */
	delete_grid_function(u);
	delete_grid_function(uh);
	delete_grid_function(u_t);

/* free the eigenvalue structure */
	delete_grid_function( eval->u_eps );
	delete_grid_function( eval->u_eps_t );
	delete_grid_function( eval->x1 );
	delete_grid_function( eval->x2 );
	delete_grid_function( eval->Ax1 );
	delete_grid_function( eval->Ax2 );
	free(eval);
	return;
      }
/* estimate largest stable dt for R-K-2 */
      dt = 2.0*info->time_factor / fabs(eval->lambda_1r);

/* do not integrate past t_max */
/*      printf("t-t_min = %e, dt = %e\n", t-t_min, dt);*/
      dt = real_min( t-t_min, dt );
/*      printf("Actual dt = %e\n", dt);*/

/* update */
      for (i=1; i<=info->n1; i++)
	{
	  uh1(i) = u1(i) + 0.5 * dt * u_t1(i);
	  uh2(i) = u2(i) + 0.5 * dt * u_t2(i);
	}
/* apply bc */

/* second stage */
      time_derivative( uh, t+0.5*dt, u_t, info, n_sign );
/* update */
      for (i=1; i<=info->n1; i++)
	{
	  u1(i) = u1(i) + dt * u_t1(i);
	  u2(i) = u2(i) + dt * u_t2(i);
	}
/* apply bc */

/* update time */
      t -= dt;
/* increment the iteration counter */
      itr++;
/*      printf("Present time = %e\n", t);*/
    } while (t > t_min+eps);
    

/* assign the new layer from the time-integration variables */
    for (i=1; i<= info->n1; i++)
      {
	x_hyp(i,k) = u1(i);
	y_hyp(i,k) = u2(i);
      }
  } /* end for all k */
    } /* end if both-sided */

/* update the bounding box */
  grid_ptr->x_min_0 = 1.e10;
  grid_ptr->x_max_0 = -1.e10;
  grid_ptr->y_min_0 = 1.e10;
  grid_ptr->y_max_0 = -1.e10;
  for (j=1; j<=info->n2; j++)
    for (i=1; i<=info->n1; i++)
      {
	grid_ptr->x_min_0 = real_min( x_hyp(i,j), grid_ptr->x_min_0 );
	grid_ptr->x_max_0 = real_max( x_hyp(i,j), grid_ptr->x_max_0 );
	grid_ptr->y_min_0 = real_min( y_hyp(i,j), grid_ptr->y_min_0 );
	grid_ptr->y_max_0 = real_max( y_hyp(i,j), grid_ptr->y_max_0 );
      }

  reset_view_port( grid_ptr );

/* successful computation */
  info->status = 1;

/*  delete the gridfunctions the present solution, the middle step, */
/*  and the time-derivative. */
  delete_grid_function(u);
  delete_grid_function(uh);
  delete_grid_function(u_t);

/* free the eigenvalue structure */
  delete_grid_function( eval->u_eps );
  delete_grid_function( eval->u_eps_t );
  delete_grid_function( eval->x1 );
  delete_grid_function( eval->x2 );
  delete_grid_function( eval->Ax1 );
  delete_grid_function( eval->Ax2 );
  free(eval);

} /* end compute_hyp_grid */


static real
scalar_prod( grid_function *u1, grid_function *u2 ){
  real sp;
  int i;

  sp = 0.0;
  for (i=1; i<= u1->n1; i++)
    {
      sp += c_1(u1,i) * c_1(u2,i) + c_2(u1,i) * c_2(u2,i);
    }

  return sp;
}

static void
power_method_1( grid_function *u, grid_function *u_t, real t, hyp_grid_data *info,
	     eigen_val *eval, real n_sign ){
  const real eps=1.e-4;
  real Ax1_norm, lambda_1r, lambda_1i, rel_diff;
  int i, n1, itr=0;
  grid_function *x1, *Ax1, *u_eps, *u_eps_t;

  Ax1 = eval->Ax1;
  x1 = eval->x1;
  u_eps = eval->u_eps;
  u_eps_t = eval->u_eps_t;

  n1 = u->n1;

/* copy the eigenvalues to get the iteration going */
  lambda_1r = eval->lambda_1r; lambda_1i = eval->lambda_1i;

/* Iterate: */
  do{
/* save the eigenvalues from the previous iteration */
    eval->lambda_1r = lambda_1r; eval->lambda_1i = lambda_1i;

/* Normalize Ax1, put the result in x1 */
    Ax1_norm = sqrt( scalar_prod( Ax1, Ax1 ) );
    for (i=1; i<=n1; i++)
      {
	c_1(x1,i) = c_1(Ax1,i)/Ax1_norm;
	c_2(x1,i) = c_2(Ax1,i)/Ax1_norm;
      }

/* Assign u_eps = u + eps*x1, compute the time_derivative corresponding to u_eps */
/* => u_eps_t */
    for (i=1; i<=n1; i++)
      {
	c_1(u_eps,i) = c_1(u,i) + eps * c_1(x1,i);
	c_2(u_eps,i) = c_2(u,i) + eps * c_2(x1,i);
      }
    time_derivative( u_eps, t, u_eps_t, info, n_sign );
/* Ax1 = (u_eps_t - u_t)/eps. */
    for (i=1; i<=n1; i++)
      {
	c_1(Ax1,i) = (c_1(u_eps_t,i) - c_1(u_t,i))/eps;
	c_2(Ax1,i) = (c_2(u_eps_t,i) - c_2(u_t,i))/eps;
      }
/* Compute the eigenvalue */
    lambda_1r = scalar_prod( x1, Ax1 );
    
/*     printf("lambda_1 = %e + i %e\n", lambda_1r, lambda_1i); */

    rel_diff = fabs(lambda_1r - eval->lambda_1r) / fabs(lambda_1r);
/* iterate until the eigenvalues have converged, but not for more than 10 iterations */
  } while( ++itr <= 10 && rel_diff > 0.05 );

/* display the eigenvalues */
/*  printf("Dominating eigenvalue converged after %i iterations:\n", itr); */
/*  printf("lambda_1 = %e\n", lambda_1r); */

}

static int
grid_normal(real *n_vec, grid_function *u, int i, hyp_grid_data *info){
  real length, xr, yr;
  int q;
  const real eps = 1.e-10;

/* r-differences */
  if (i == 1){
    if (info->periodic)
      {
	xr = u1(2) - u1(u->n1-1);
	yr = u2(2) - u2(u->n1-1);
      }
    else
      {
	xr = -u1(i) + u1(i+1);
	yr = -u2(i) + u2(i+1);
      }
  }
  else if (i == u->n1){
    if (info->periodic)
      {
	xr = u1(2) - u1(u->n1-1);
	yr = u2(2) - u2(u->n1-1);
      }
    else
      {
	xr = u1(i) - u1(i-1);
	yr = u2(i) - u2(i-1);
      }
  }
  else{
    xr = u1(i+1) - u1(i-1);
    yr = u2(i+1) - u2(i-1);
  }

/* compute normal */
  n_vec[0] =  yr;
  n_vec[1] = -xr;

/* apply boundary conditions */
  if (i==1 && hyp_bc(1) != 0){
    if (hyp_bc(1) == 1)
      n_vec[0] = 0.0;
    else if (hyp_bc(1) == 2)
      n_vec[1] = 0.0;
  }
  else if (i==info->n1 && hyp_bc(2) != 0){
    if (hyp_bc(2) == 1)
      n_vec[0] = 0.0;
    else if (hyp_bc(2) == 2)
      n_vec[1] = 0.0;
  }

/* normalize */
  length = sqrt( n_vec[0] * n_vec[0] + n_vec[1] * n_vec[1] );
  if (length > eps){
    for (q=0; q<2; q++) n_vec[q] = n_vec[q]/length;
    return OK;
  }
  else
    return ERROR;
}

static real
mean_curvature( grid_function *u, real_array_1d *kappa_, int periodic ){
  real length;
  real x1[2], x11[2];
  real dr, dr2, max_value;
  int i0, j0, i, j;

#define kappa(i) compute_index_1d(kappa_, i)

  dr = 1.0/(u->n1-1); dr2 = dr*dr;

  max_value = -1.e10;
  for (i0=1; i0<= u->n1; i0++)
    {
      i = i0; j = j0;

/* r-differences */
      if (i == 1)
	{
	  if (periodic)
	    {
	      x1[0] = 0.5*(u1(2) - u1(u->n1-1))/dr;
	      x1[1] = 0.5*(u2(2) - u2(u->n1-1))/dr;
	    }
	  else
	    {
	      x1[0] = (u1(i+1) - u1(i))/dr;
	      x1[1] = (u2(i+1) - u2(i))/dr;
	    }
	}
      else if (i == u->n1)
	{
	  if (periodic)
	    {
	      x1[0] = 0.5*(u1(2) - u1(u->n1-1))/dr;
	      x1[1] = 0.5*(u2(2) - u2(u->n1-1))/dr;
	    }
	  else
	    {
	      x1[0] = (u1(i) - u1(i-1))/dr;
	      x1[1] = (u2(i) - u2(i-1))/dr;
	    }
	}
      else
	{
	  x1[0] = 0.5*(u1(i+1) - u1(i-1))/dr;
	  x1[1] = 0.5*(u2(i+1) - u2(i-1))/dr;
	}

/* rr-difference */
      if (periodic && (i == 1 || i == u->n1) )
	{
	  x11[0] = (u1(2) - 2.0*u1(1) + u1(u->n1-1))/dr2;
	  x11[1] = (u2(2) - 2.0*u2(1) + u2(u->n1-1))/dr2;
	}
      else
	{
/* shift i and j away from any boundaries */
	  i = real_max(2, i); i = real_min(u->n1-1, i);

	  x11[0] = (u1(i+1) - 2.0*u1(i) + u1(i-1))/dr2;
	  x11[1] = (u2(i+1) - 2.0*u2(i) + u2(i-1))/dr2;
	}

      length = sqrt(sqr(x1[0]) + sqr(x1[1]));
/* curvature */
      kappa(i0) = (x1[0]*x11[1] - x1[1]*x11[0])/(length*length*length);

      if (kappa(i0) > max_value) max_value = kappa(i0);
    }
  return max_value;
#undef sp
}

static real
smooth_v(real t, hyp_grid_data *info){
  real kappa_max, smooth, fmin, v0;

#define log_cosh(t0)  -info->curv_factor/(2.0*smooth) * \
  log( cosh(smooth*(t0-fmin))/cosh(smooth*(t0-kappa_max)) )
#define smooth_vel(t0)  log_cosh(t0) + 1.0 - v0

/*  return real_max( info->v_min, 1 - info->curv_factor * t);*/
  
  fmin = -10.0; 
  smooth = 0.25; 
  kappa_max = (1.0-info->v_min)/info->curv_factor; 
  v0 = log_cosh(0.0); 

  return smooth_vel(t);
}

static void
time_derivative(grid_function *u, real t, grid_function *u_t, hyp_grid_data *info, real n_sign){
  int i;
  real n_vec[2], v, ds_dt, t_stretch, dummy;
  grid_function *u_tmp;
  real_array_1d *kappa_;

#define u_tmp1(i) compute_index_1d(u_tmp->x_, i)
#define u_tmp2(i) compute_index_1d(u_tmp->y_, i)

  u_tmp  = new_grid_function( u->n1 );
  kappa_ = create_real_array_1d( u->n1 );

/* 2. compute the mean curvature */
  mean_curvature( u, kappa_, info->periodic );
/* stretching in the normal direction? */
  uniform_to_stretch(t, info->normal_stretch, &t_stretch, &ds_dt, &dummy );
/* The time-derivative of (x,y)^T is computed as follows: */
  for (i=1; i<=u->n1; i++)
    {
/* 1. Compute the normals */
      grid_normal( n_vec, u, i, info );
/* 3. Compute the velocity function */
      v = info->width * ds_dt * smooth_v(n_sign*kappa(i), info);
      u_tmp1(i) = n_sign * v * n_vec[0];
      u_tmp2(i) = n_sign * v * n_vec[1];
    }
/* 4. smooth the time-derivative for all interior points */
  for (i=2; i<=u->n1-1; i++)
    { 
       u_t1(i) = (1.0-info->smooth_factor) * u_tmp1(i) +  
	 info->smooth_factor * (u_tmp1(i-1) + u_tmp1(i) + u_tmp1(i+1))/3.0;
       u_t2(i) = (1.0-info->smooth_factor) * u_tmp2(i) +  
	 info->smooth_factor * (u_tmp2(i-1) + u_tmp2(i) + u_tmp2(i+1))/3.0;
     } 

/* smooth along the boundary for the boundary values */
/* low i */
  i = 1;
  if (hyp_bc(1) == 0)
    {
      if (info->periodic)
	{
	  u_t1(1) = (1.0-info->smooth_factor) * u_tmp1(1) +  
	    info->smooth_factor * (u_tmp1(u->n1-1) + u_tmp1(1) + u_tmp1(2))/3.0;
	  u_t2(1) = (1.0-info->smooth_factor) * u_tmp2(1) +  
	    info->smooth_factor * (u_tmp2(u->n1-1) + u_tmp2(1) + u_tmp2(2))/3.0;
	}
      else
	{
	  u_t1(i) = (1.0-info->smooth_factor) * u_tmp1(i) +  
	    info->smooth_factor * (u_tmp1(i) + u_tmp1(i+1))/2.0;
	  u_t2(i) = (1.0-info->smooth_factor) * u_tmp2(i) +  
	    info->smooth_factor * (u_tmp2(i) + u_tmp2(i+1))/2.0;
	}
    } /* end free bc */
  else if (hyp_bc(1) == 1)
    {
      u_t1(i) = 0.0;
      u_t2(i) = (1.0-info->smooth_factor) * u_tmp2(i) +  
	info->smooth_factor * (u_tmp2(i) + u_tmp2(i+1))/2.0;
    } /* end x=constant bc */
  else if (hyp_bc(1) == 2)
    {
      u_t1(i) = (1.0-info->smooth_factor) * u_tmp1(i) +  
	info->smooth_factor * (u_tmp1(i) + u_tmp1(i+1))/2.0;
      u_t2(i) = 0.0;
    } /* end y=constant bc */

/* high i */
  i = u->n1;
  if (hyp_bc(2) == 0)
    {
      if (info->periodic)
	{
	  u_t1(u->n1) = (1.0-info->smooth_factor) * u_tmp1(1) +  
	    info->smooth_factor * (u_tmp1(u->n1-1) + u_tmp1(1) + u_tmp1(2))/3.0;
	  u_t2(u->n1) = (1.0-info->smooth_factor) * u_tmp2(1) +  
	    info->smooth_factor * (u_tmp2(u->n1-1) + u_tmp2(1) + u_tmp2(2))/3.0;
	}
      else
	{
	  u_t1(i) = (1.0-info->smooth_factor) * u_tmp1(i) +  
	    info->smooth_factor * (u_tmp1(i) + u_tmp1(i-1) )/2.0;
	  u_t2(i) = (1.0-info->smooth_factor) * u_tmp2(i) +  
	    info->smooth_factor * (u_tmp2(i) + u_tmp2(i-1) )/2.0;
	}
    }
  else if (hyp_bc(2) == 1)
    {
      u_t1(i) = 0.0;
      u_t2(i) = (1.0-info->smooth_factor) * u_tmp2(i) +  
	info->smooth_factor * (u_tmp2(i) + u_tmp2(i-1) )/2.0;
  }
  else if (hyp_bc(2) == 2)
    {
      u_t1(i) = (1.0-info->smooth_factor) * u_tmp1(i) +  
	info->smooth_factor * (u_tmp1(i) + u_tmp1(i-1) )/2.0;
      u_t2(i) = 0.0;
  }

/* free the memory */
  delete_grid_function( u_tmp );
  delete_real_array_1d( kappa_ );

}

hyp_grid_data *
init_hyp_grid(input_output *io_ptr, generic_mapping *grid_ptr, generic_curve *curve_ptr){
  hyp_grid_data *info;
  
/* allocate memory */
  info = (hyp_grid_data *) malloc( sizeof(hyp_grid_data) );

/* general data */
  grid_ptr->grid_type = 8; 
  grid_ptr->grid_mapping = hyp_grid_mapping; 
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->mapping_data_ptr = info; 
  grid_ptr->cleanup_mapping_data = delete_hyp_grid;
  grid_ptr->set_mapping = set_hyp_grid; 
  grid_ptr->write_specific_mapping = write_hyp_map; 
  grid_ptr->copy_mapping_data = copy_hyp_map;
  
/* default number of grid points */
  grid_ptr->r_points = 30;
  grid_ptr->s_points = 5;

/* and periodicity */
  grid_ptr->r_period = curve_ptr->periodic;
  grid_ptr->s_period = 0;

/* specific data */
  info->curve_ptr = curve_ptr;
  info->width = 0.5; 
  info->normal_stretch = NULL;
  info->tangential_stretch = NULL;
  info->v_min = 0.25;
  info->curv_factor = 0.25;
  info->smooth_factor = 0.5;
  info->time_factor = 1.0;
  info->status = 0;
  info->hyp_bc_ = create_int_array_1d(2);
  info->flip = 0;
  info->periodic = grid_ptr->r_period;
  info->c_grid = FALSE; /* the default is NOT a C-grid */
  info->n1 = grid_ptr->r_points;
  info->n2 = grid_ptr->s_points;
  info->n_cut = 0;  
  info->n_curve = info->n1;

  if (info->periodic)
    printf("The curve is marked as periodic!\n");

/* both-sided stuff */
  info->both_sides = FALSE;

/* initialize bc */
  hyp_bc(1) = hyp_bc(2) = 0;

/* initially, copy the bounding box from the curve */
  grid_ptr->x_min_0 = curve_ptr->x_min;
  grid_ptr->x_max_0 = curve_ptr->x_max; 
  grid_ptr->y_min_0 = curve_ptr->y_min; 
  grid_ptr->y_max_0 = curve_ptr->y_max; 

/* no surface arrays to begin with */
  info->x_ = NULL;
  info->y_ = NULL;

/* compute the translated/rotated/scaled bounding box */
  reset_view_port( grid_ptr );

  return info;
}

static void *
delete_hyp_grid( void *mapping_data_ptr ){
  hyp_grid_data *info;
  
/* cast the pointer */
  info = (hyp_grid_data *) mapping_data_ptr;

/* free the boundary condition array */
  info->hyp_bc_ = delete_int_array_1d(info->hyp_bc_);

/* free the stretchings */
  info->tangential_stretch = delete_generic_stretching( info->tangential_stretch );
  info->normal_stretch = delete_generic_stretching( info->normal_stretch );

/* free the arrays */
  info->x_ = delete_real_array_2d( info->x_ );
  info->y_ = delete_real_array_2d( info->y_ );

/* free the structure */
  free(info);

  return NULL;
}

#define R0 ((i_loc-1)   * info->d1)
#define R1 (i_loc * info->d1)
#define ALPHA_0(r)   ((R1 - r)/info->d1)
#define D_ALPHA_0(r) (-1.0/info->d1)
#define ALPHA_1(r)   ((r - R0)/info->d1)
#define D_ALPHA_1(r) (1.0/info->d1)

#define S0 ((j_loc-1)  * info->d2)
#define S1 (j_loc * info->d2)
#define BETA_0(s)   ((S1 - s)/info->d2)
#define D_BETA_0(s) (-1.0/info->d2)
#define BETA_1(s)   ((s - S0)/info->d2)
#define D_BETA_1(s) (1.0/info->d2)

#define LIN(v,r,s) (ALPHA_0(r) * BETA_0(s) * v(i_loc  ,j_loc  ) +  \
                    ALPHA_1(r) * BETA_0(s) * v(i_loc+1,j_loc  ) +  \
                    ALPHA_0(r) * BETA_1(s) * v(i_loc  ,j_loc+1) +  \
                    ALPHA_1(r) * BETA_1(s) * v(i_loc+1,j_loc+1))

#define LIN_R(v,r,s) (D_ALPHA_0(r) * BETA_0(s) * v(i_loc  ,j_loc   ) +  \
                      D_ALPHA_1(r) * BETA_0(s) * v(i_loc+1,j_loc   ) +  \
                      D_ALPHA_0(r) * BETA_1(s) * v(i_loc  ,j_loc+1 ) +  \
                      D_ALPHA_1(r) * BETA_1(s) * v(i_loc+1,j_loc+1 ) )

#define LIN_S(v,r,s) (ALPHA_0(r) * D_BETA_0(s) * v(i_loc  ,j_loc  ) +  \
                      ALPHA_1(r) * D_BETA_0(s) * v(i_loc+1,j_loc  ) +  \
                      ALPHA_0(r) * D_BETA_1(s) * v(i_loc  ,j_loc+1) +  \
                      ALPHA_1(r) * D_BETA_1(s) * v(i_loc+1,j_loc+1) )


static int 
hyp_grid_mapping( grid_point *gp, void *data_ptr){
  hyp_grid_data *info;
  int i_loc, j_loc;

/* cast the data pointer to the right type */
  info = (hyp_grid_data *) data_ptr;

/* get the nearest grid point below and to the left of (r,s) */
  i_loc = int_max( 1, 1 + (int) ( (info->n1-1)*gp->r ) );
  j_loc = int_max( 1, 1 + (int) ( (info->n2-1)*gp->s ) );

/* prevent overflow */
  i_loc = int_min( i_loc, info->n1-1 );
  j_loc = int_min( j_loc, info->n2-1 );

/* evaluate the tri-linear mapping */
  gp->x  = LIN(x_hyp, gp->r, gp->s);
  gp->y  = LIN(y_hyp, gp->r, gp->s);

  gp->xr = LIN_R(x_hyp, gp->r, gp->s);
  gp->yr = LIN_R(y_hyp, gp->r, gp->s);

  gp->xs = LIN_S(x_hyp, gp->r, gp->s);
  gp->ys = LIN_S(y_hyp, gp->r, gp->s);

  return 0;
}

static void 
write_hyp_map( int32 dir, void *data_ptr ){
  hyp_grid_data *info;
  int32 stretch_dir;

/* cast the data pointer to the right type */
  info = (hyp_grid_data *) data_ptr;

/* save the name of the curve */
  hput_string( info->curve_ptr->curve_name, "curve_name", dir );

  hput_real( info->width, "width", dir );
  hput_real( info->curv_factor, "curv-factor", dir );
  hput_real( info->v_min, "v_min", dir );
  hput_real( info->smooth_factor, "smooth-factor", dir );
  hput_real( info->time_factor, "time-factor", dir );

/* save the sizes */
  hput_int( info->n1, "n1", dir );
  hput_int( info->n2, "n2", dir );
  hput_int( info->n_curve, "n-curve", dir );
  hput_int( info->n_cut, "n-cut", dir );

  hput_int( info->flip, "flip", dir );
  hput_int( info->both_sides, "both-sides", dir );
  hput_int( info->periodic, "periodic", dir );
  hput_int( info->c_grid, "c-grid", dir );

  hput_real( info->d1, "d1", dir );
  hput_real( info->d2, "d2", dir );
  hput_real( info->x_branch, "x-branch", dir );
  hput_real( info->y_branch, "y-branch", dir );

/* save the arrays */
  hput_real_array_2d( info->x_, "x", dir );
  hput_real_array_2d( info->y_, "y", dir );
  hput_int_array_1d( info->hyp_bc_, "hyp-bc", dir );

  hput_int( info->status, "status", dir );

/* r-stretching */
  if ((stretch_dir = create_dir("tangential-stretching", 
				"tangential-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( stretch_dir, info->tangential_stretch );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }
  else
    printf("Error: write_hyp_grid: unable to create a directory for "
	   "`tangential-stretching'.\n");

/* s-stretching */
  if ((stretch_dir = create_dir("normal-stretching", 
				"normal-stretching", dir)) != -1){
/* save the generic and specific curve data */
    write_stretch( stretch_dir, info->normal_stretch );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }
  else
    printf("Error: write_hyp_grid: unable to create a directory for "
	   "`normal-stretching'.\n");
}

int 
read_hyp_map( int32 dir, generic_mapping *grid_ptr, generic_curve *first_c ){
  hyp_grid_data *info;
  char *curve_name;
  int32 stretch_dir;

  info = (hyp_grid_data *) malloc( sizeof(hyp_grid_data) );

/* read the name */
  curve_name = hget_string( "curve_name", dir );

/* look for this curve in the curve list */
  for (info->curve_ptr = first_c; 
       info->curve_ptr != NULL && 
       strcmp(info->curve_ptr->curve_name, curve_name) != 0;
       info->curve_ptr = info->curve_ptr->next);

/* was the curve found? */
  if (info->curve_ptr == NULL){
    printf("ERROR: read_hyp_map. Could not find curve with name: %s\n",
	   curve_name);
    free( curve_name );
    free( info );
    return FALSE;
  }
  free( curve_name );

  hget_real( &(info->width), "width", dir );
  hget_real( &(info->curv_factor), "curv-factor", dir );
  hget_real( &(info->v_min), "v_min", dir );
  hget_real( &(info->smooth_factor), "smooth-factor", dir );
  hget_real( &(info->time_factor), "time-factor", dir );

/* read the sizes */
  hget_int( &(info->n1), "n1", dir );
  hget_int( &(info->n2), "n2", dir );
  if (!hget_int( &(info->n_curve), "n-curve", dir ))
    info->n_curve = info->n1;
  if (!hget_int( &(info->n_cut), "n-cut", dir ))
    info->n_cut = 0;

  hget_int( &(info->flip), "flip", dir );
  if (!hget_int( &(info->both_sides), "both-sides", dir ))
    info->both_sides = FALSE;
  if (!hget_int( &(info->periodic), "periodic", dir ))
    info->periodic = FALSE;
  if (!hget_int( &(info->c_grid), "c-grid", dir ))
    info->c_grid = FALSE;

  hget_real( &(info->d1), "d1", dir );
  hget_real( &(info->d2), "d2", dir );
  if (!hget_real( &(info->x_branch), "x-branch", dir ))
    info->x_branch = 1.0;
  if (!hget_real( &(info->y_branch), "y-branch", dir ))
    info->y_branch = 0.0;

/* read the arrays */
  info->x_ = hget_real_array_2d( "x", dir );
  info->y_ = hget_real_array_2d( "y", dir );
  info->hyp_bc_ = hget_int_array_1d( "hyp-bc", dir );

  hget_int( &(info->status), "status", dir );

/* read the stretchings */
  if ((stretch_dir = locate_dir("tangential-stretching", dir)) != -1){
/* save the generic and specific curve data */
    info->tangential_stretch = read_stretching( stretch_dir );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }

  if ((stretch_dir = locate_dir("normal-stretching", dir)) != -1){
/* save the generic and specific curve data */
    info->normal_stretch     = read_stretching( stretch_dir );
/* release the sub-directory */
    Vdetach(stretch_dir);
  }

/* setup the pointers to all specific functions */
  grid_ptr->grid_mapping = hyp_grid_mapping; 
  grid_ptr->inverse_grid_mapping = NULL;
  grid_ptr->mapping_data_ptr = info; 
  grid_ptr->cleanup_mapping_data = delete_hyp_grid;
  grid_ptr->set_mapping = set_hyp_grid; 
  grid_ptr->write_specific_mapping = write_hyp_map; 
  grid_ptr->copy_mapping_data = copy_hyp_map;

  return OK;
}

static void *
copy_hyp_map( void *mapping_data_ptr ){
  hyp_grid_data *info, *old_info;
  int i, j;
  
  info = (hyp_grid_data *) malloc( sizeof(hyp_grid_data) );

  old_info = (hyp_grid_data *)mapping_data_ptr;

/* only make a shallow copy of the curve. */
  info->curve_ptr = old_info->curve_ptr;

/* remember that this curve is now used by both mappings */
  info->curve_ptr->used_by_grid++;

  info->width = old_info->width;
  info->curv_factor = old_info->curv_factor;
  info->v_min = old_info->v_min;
  info->smooth_factor = old_info->smooth_factor;
  info->time_factor = old_info->time_factor;

  info->n1 = old_info->n1;
  info->n2 = old_info->n2;
  info->flip = old_info->flip;

  info->d1 = old_info->d1;
  info->d2 = old_info->d2;

#define old_x_hyp(i,j) compute_index_2d(old_info->x_,i,j)
#define old_y_hyp(i,j) compute_index_2d(old_info->y_,i,j)
  info->x_ = create_real_array_2d( info->n1, info->n2 );
  info->y_ = create_real_array_2d( info->n1, info->n2 );

/* copy the array elements */
  for (i=1; i<=info->n1; i++)
    for (j=1; j<=info->n2; j++){
      x_hyp(i,j) = old_x_hyp(i,j);
      y_hyp(i,j) = old_y_hyp(i,j);
    }

#undef old_x_hyp
#undef old_y_hyp

#define old_hyp_bc(i) compute_index_1d(old_info->hyp_bc_,i)
  info->hyp_bc_ = create_int_array_1d(2);

  /* copy the elements */
  hyp_bc(1) = old_hyp_bc(1);
  hyp_bc(2) = old_hyp_bc(2);
#undef old_hyp_bc

  info->status = old_info->status;

/* make a copy of the stretching functions */
  info->tangential_stretch = copy_stretching( old_info->tangential_stretch );
  info->normal_stretch = copy_stretching( old_info->normal_stretch );

  return (void *)info;
}

static int
newton_search(real xp, real yp, real *r_loc, real *n_loc, generic_curve *project_curve){
  curve_point cp;
  int i, max_iter=10;
  real t, a11, a12, a21, a22, b1, b2, tol=1.e-10, det, dr, dt, res;
  
  for (i=0; i<max_iter; i++)
    {
      cp.r = *r_loc;
      t = *n_loc;
      curve_function(&cp, project_curve);
      b1 = xp - (cp.x - t*cp.yr);
      b2 = yp - (cp.y + t*cp.xr);
      if (fabs(b1) + fabs(b2) < tol) break;

      a11 = cp.xr;
      a12 = -cp.yr;
      a21 = cp.yr;
      a22 = cp.xr;

      det = a11*a22 - a12*a21;
      dr = ( a22*b1 - a12*b2)/det;
      dt = (-a21*b1 + a11*b2)/det;
      
      *r_loc += dr;
      *n_loc += dt;
    }
  if ((res=fabs(b1) + fabs(b2)) >= tol) 
    {
      printf("Warning newton_search converged poorly, res=%e\n", res);
      return 0;
    }
  else
    return 1;

}


static int
invert_curve(real xp, real yp, real *r_loc, real *n_loc, generic_curve *project_curve){
  real dr, r_close, tmp_dist, min_dist;
  int i, n_proj=100;
  curve_point cp;

  dr = 1.0/((real) n_proj - 1);
/* find the closest point on the discretized curve */
  r_close = -1.0;
  min_dist = 1.e10;
  for (i=0; i<n_proj; i++)
    {
      cp.r = i*dr;
      curve_function( &cp, project_curve );
      if ((tmp_dist = sqr(cp.x-xp) + sqr(cp.y-yp)) < min_dist)
	{
	  min_dist = tmp_dist;
	  r_close = i*dr;
	}
    } /* end for */

/* do a Newton iteration to invert the curve */
  *r_loc = r_close;
  *n_loc = 0.0;
  return newton_search(xp, yp, r_loc, n_loc, project_curve);
}

static void
compute_projection(input_output *io_ptr, generic_mapping *grid_ptr, 
		   generic_curve *project_curve ){
  char *prompt;
  int icom, quit=0, side = 0, n2, i, j;
  real *x_gap, *y_gap;
  hyp_grid_data *info;
  curve_point cp;
  real s, w, r_loc, n_loc;
  const real alpha=5.0;

#define x_gap(i) compute_index_1d(x_gap_,i)
#define y_gap(i) compute_index_1d(y_gap_,i)

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
      side = 1;
      quit = 1;
      break;

    case RIGHT:
      side = 2;
      quit = 1;
      break;

    case CANCEL:
      return;
      break;

    default:
      break;

    }
  } while (!quit);

/* do the projection for the selected face */

/* get the specific data for the hyperbolic mapping */
  info = (hyp_grid_data *) grid_ptr->mapping_data_ptr; 

/* left or right face */
  n2 = info->n2; 
/* pick side */
  if (side == 1) 
    i = 1; 
  else 
    i = info->n1; 

/* allocate space for the projected curve */
  x_gap = (real *) malloc( (n2+1)*sizeof(real) ); 
  y_gap = (real *) malloc( (n2+1)*sizeof(real) ); 

/* do the projection */
  printf("Projecting...\n"); 
  j = 1;
/* invert the curve for (x_hyp(i,j), y_hyp(i,j)) */
  if (invert_curve(x_hyp(i,j), y_hyp(i,j), &r_loc, &n_loc, project_curve)){
/* evaluate the curve */
    cp.r = r_loc;
    curve_function(&cp, project_curve); 
/* save the gap */
    x_gap[j] = cp.x - x_hyp(i,j); 
    y_gap[j] = cp.y - y_hyp(i,j); 

    for (j=2; j<=n2; j++){ 
      if (newton_search(x_hyp(i,j), y_hyp(i,j), &r_loc, &n_loc, project_curve))
	{
/* evaluate the curve */
	  cp.r = r_loc;
	  curve_function(&cp, project_curve); 
/* save the gap */
	  x_gap[j] = cp.x - x_hyp(i,j); 
	  y_gap[j] = cp.y - y_hyp(i,j); 
	}
      else
	{
	  x_gap[j] = 0.0;
	  y_gap[j] = 0.0;
	}
    } /* end for */
  }
  else{ 
    printf("compute_projection: Warning: unable to project the point (%e, %e) " 
	   "onto the curve `%s'\n", x_hyp(i,j), y_hyp(i,j),  
	   project_curve->curve_name); 
    for (j=1; j<=n2; j++)
      {
	x_gap[j] = 0.0; 
	y_gap[j] = 0.0; 
      }
  }

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

/* Move the i=1 grid line eventhough it defines the shape of the boundary. In this way */
/* the grid can stay smooth all the way up to the boundary when it is stretched in  */
/* the j-direction */

    for (j=1; j<=info->n2; j++){ 
      x_hyp(i,j) = x_hyp(i,j) + w*x_gap[j]; 
      y_hyp(i,j) = y_hyp(i,j) + w*y_gap[j]; 
    } 
  } /* end for i */

/* cleanup */
  free( x_gap ); 
  free( y_gap ); 

/* end left or right face */

}
