static int inside_region(real xp, real yp, component_grid *grid_ptr){
  int i, j, intersect;
  const real eps=1.0e-7;
  real x_star;

  intersect = 0;

/* a point is outside if it is outside of the boudning box */
  if ( xp < grid_ptr->x_min || xp > grid_ptr->x_max || 
      yp < grid_ptr->y_min || yp > grid_ptr->y_max)
    return 0;

  if ( !grid_ptr->s_period ){
    j = range(1,2);
    for( i=range(1,1); i<range(2,1); i++)
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
/*    printf("Inside_region after range(1,2): intersect = %d\n",intersect);*/
  
    j = range(2,2);
    for( i=range(1,1); i<range(2,1); i++)
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
/*    printf("Inside_region after range(2,2): intersect = %d\n",intersect);*/
  }

  if ( !grid_ptr->r_period ){
    i = range(1,1);
    for( j=range(1,2); j<range(2,2); j++)
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

/*    printf("Inside_region after range(1,1): intersect = %d\n",intersect);*/
    
    i = range(2,1);
    for( j=range(1,2); j<range(2,2); j++)
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
/*  printf("Inside_region after range(2,1): intersect = %d\n",intersect);*/
  }  

/* The point (xp,yp) was outside if intersect is even */  
  return ((intersect % 2 == 0)? 0 : 1);
}


