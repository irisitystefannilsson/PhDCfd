/**
 **  Lists, grids, interpolations, ...
 **/
/*bug: underscore explicitly added to fortran names, this is not portable! */

#include "mappings.h"
 
#ifdef F_NO_UNDERSCORE
extern void
dslucs( int *n, real *b, real *x, int *nelt, int *ia, int *ja, real *a, int *isym,
        int *itol, real *tol, int *itmax, int *iter, real *err, int *ierr, int *iunit,
        real *rwork, int *lenw, int *iwork, int *leniw); /* I love FORTRAN... NOT */
#else
extern void
dslucs_( int *n, real *b, real *x, int *nelt, int *ia, int *ja, real *a, int *isym,
        int *itol, real *tol, int *itmax, int *iter, real *err, int *ierr, int *iunit,
        real *rwork, int *lenw, int *iwork, int *leniw); /* I love FORTRAN... NOT */
#endif
 pointlist aspointlist(name,npoints,x,y)
	char	*name;
	int		npoints;
	real	*x;
	real	*y;
	
{
	pointlist	answer;
	
	answer.name = name;
	answer.npoints = npoints;
	answer.x = x;
	answer.y = y;
	return(answer);
}
	

complexcurve pointlist2complex(ws)
	pointlist	ws;
{
	complexcurve	cc;
	complex			*z;
	real			*x, *y;
	
	cc.npoints = ws.npoints;
	cc.z = (complex *)malloc(ws.npoints * sizeof(complex));
	x = ws.x;
	y = ws.y;
	for (z = cc.z; z < cc.z + ws.npoints; z++)
	{
		z->re = *x++;
		z->im = *y++;
	}
	return(cc);
}

real length(x,y)
	real	x,y;
{
	real	l;
	
	l = sqrt(pow(x,2.0)+pow(y,2.0));
	return(l);
}

real tan_angle(x1,y1,x2,y2)
	real	x1,y1,x2,y2;
{
	real	x, tb;
	
	x = (x1*y2-y1*x2)/length(x1,y1)/length(x2,y2);
	tb = x/sqrt(1-pow(x,2.0));
	return(tb);
}

CPCspline 
mkCPCspline(pointlist points){
	real	*a, *r, *s;
	int		*ia, *ja;
	int		i,ie,next,prev;
	real	*l;
	int		n, nn, ndim;
	CPCspline	spline;
	real	tanbeta_12;
	real	*x;
	real	*y,yprev;
	
	int	nelt,isym,itol,itmax,iter,ierr,iunit,lenw,leniw,*iwork;
	real    tol,err,*rwork;
	
	spline.name = (char *)malloc(1);
	strcpy(spline.name,"\0");
	spline.npoints = points.npoints;
	spline.ndimensions = 2;
	spline.ulist = NULL;
	spline.B = (real *) malloc((spline.npoints-1) * spline.ndimensions * 4 * sizeof(real));
	
	nn = points.npoints-1;
	ndim = spline.ndimensions;
	a = (real *)calloc(nn*3,sizeof(real));
	ia = (int *)calloc(nn*3,sizeof(int));
	ja = (int *)calloc(nn*3,sizeof(int));
	r = (real *)calloc(nn,sizeof(real));
	s = (real *)calloc(nn,sizeof(real));memset(s,0,nn*sizeof(real));
	l = (real *)calloc(nn,sizeof(real));
	x = points.x;
	y = points.y;
	
	for (n=0; n < nn; n++){
	  next=n+1;
	  l[n] = length(x[next]-x[n],y[next]-y[n]);
	}
	for (ie=0,n=0; n < nn; n++){
	  if (n==nn-1) 
	    next=0;
	  else
	    next=n+1;

	  if (n==0){ 
	    prev=nn-1;
	    yprev=y[prev]-2*((double)M_PI);
	  }
	  else{
	    prev=n-1;
	    yprev=y[prev];
	  }
		
	  tanbeta_12 = tan_angle(x[n]-x[prev],y[n]-yprev,x[n+1]-x[n],y[n+1]-y[n]);
	  a[ie] = 1/l[prev];
	  ia[ie]=n+1;
	  ja[ie++]=prev+1;

	  a[ie] = 1/l[n];
	  ia[ie]=n+1;
	  ja[ie++]=next+1;

	  a[ie] = 2 * (a[ie-1] + a[ie-2]);
	  ia[ie]=n+1;
	  ja[ie++]=n+1;

	  r[n] = 6 * tanbeta_12;
	}

	nelt = nn*3;
	isym = 0;
	itol = 2;
	tol = 1.0e-9;
	itmax = 600;
	iunit = 0;
	lenw = nelt+nn*11;
	leniw = nelt+7*nn+14;
	rwork = (real *)calloc(lenw,sizeof(real));
	iwork = (int *)calloc(leniw,sizeof(int));
#ifdef F_NO_UNDERSCORE
	dslucs(&nn,r,s,&nelt,ia,ja,a,&isym,&itol,&tol,&itmax,&iter,&err,&ierr,
	       &iunit,rwork,&lenw,iwork,&leniw);
#else
	dslucs_(&nn,r,s,&nelt,ia,ja,a,&isym,&itol,&tol,&itmax,&iter,&err,&ierr,
	       &iunit,rwork,&lenw,iwork,&leniw);
#endif
	
	for (i=0,n=0; n < nn; n++)
	{
	    next=n+1;
		if (n==nn-1) next=0;
		spline.B[i++] = x[n];
		spline.B[i++] = x[next];
		spline.B[i++] = (x[next]-x[n]+(2*s[n]+s[next])*(y[n+1]-y[n])/(6*l[n]));
		spline.B[i++] = (x[next]-x[n]-(s[n]+2*s[next])*(y[n+1]-y[n])/(6*l[n]));
		spline.B[i++] = y[n];
		spline.B[i++] = y[n+1];
		spline.B[i++] = (y[n+1]-y[n]+(2*s[n]+s[next])*(x[next]-x[n])/(6*l[n]));
		spline.B[i++] = (y[n+1]-y[n]-(s[n]+2*s[next])*(x[next]-x[n])/(6*l[n]));
	}

	free(a);
	free(ia);
	free(ja);
	free(r);
	free(s);
	free(l);
	free(rwork);
	free(iwork);
	
	return(spline);
}

void reparametrizeCPCspline(spline,parlist)
	CPCspline	*spline;
	real		*parlist;
{
	if (spline->ulist == NULL)
		spline->ulist = (real *) malloc(spline->npoints * sizeof(real));
	memcpy(spline->ulist,parlist,spline->npoints * sizeof(real));
}

real F1(u)
	real	u;
{
	real	f;
	f = 2*pow(u,3.0)-3*pow(u,2.0)+1;
	return(f);
}

real F2(u)
	real	u;
{
	real	f;
	f = -2*pow(u,3.0)+3*pow(u,2.0);
	return(f);
}

real F3(u)
	real	u;
{
	real	f;
	f = pow(u,3.0)-2*pow(u,2.0)+u;
	return(f);
}

real F4(u)
	real	u;
{
	real	f;
	f = pow(u,3.0)-pow(u,2.0);
	return(f);
}

real dF1(u)
	real	u;
{
	real	f;
	f = 6*pow(u,2.0)-6*u;
	return(f);
}

real dF2(u)
	real	u;
{
	real	f;
	f = -6*pow(u,2.0)+6*u;
	return(f);
}

real dF3(u)
	real	u;
{
	real	f;
	f = 3*pow(u,2.0)-4*u+1;
	return(f);
}

real dF4(u)
	real	u;
{
	real	f;
	f = 3*pow(u,2.0)-2*u;
	return(f);
}

real *
CPCp(real	*pm, real u, CPCspline spline){
	real	*B;
        real	dfdu;
        real	fn;
	int	i, j;
	int	n,ndim;
        real	u0;
	
	ndim = spline.ndimensions;
	if (spline.ulist == NULL)
	{
		i = floor(u);
		u = u-i;
	}
	else
	{
		u0=u;
		for (i=0;spline.ulist[i+1] < u;i++);
		u = (u-spline.ulist[i])/(spline.ulist[i+1]-spline.ulist[i]);
		B = spline.B + i*ndim*4 + 4;
		for (n=0;n<5;n++)
		{
			fn = u0-(B[0]*F1(u)+B[1]*F2(u)+B[2]*F3(u)+B[3]*F4(u));
			dfdu = B[0]*dF1(u)+B[1]*dF2(u)+B[2]*dF3(u)+B[3]*dF4(u);
			u=u+fn/dfdu;
		}
	}
	for (j=0;j<ndim;j++)
	{
		B = spline.B + i*ndim*4 + j*4;
		pm[j] = B[0]*F1(u)+B[1]*F2(u)+B[2]*F3(u)+B[3]*F4(u);
	}
	return(pm);
}
