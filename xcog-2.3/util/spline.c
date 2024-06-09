#include "real.h"
#include "spline.h"

void spline(real x[], real y[], int n, real yp1, real ypn, real y2[]){
  int i,k;
  real p,qn,sig,un,*u;

  u=vector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free_vector(u,1,n-1);
}

/* observe: 
   first argument is the nodepoints for the independent variable, 
   second argument is the nodepoints for the dependent variable,
   third argument is the spline coefficients 
*/
void cseval(real xa[], real ya[], real y2a[], int n, real x, real *y
	    , real *yp, real *ypp){
  static int klo,khi,k,called=0;
  static real h;
  real a,b,ap,bp;

/* can we use the preious interval or not? */
  if (!called || khi > n || x < xa[klo] || x > xa[khi]){
    called=1;
    klo=1;
    khi=n;
    while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if (xa[k] > x) khi=k;
      else klo=k;
    }
  }

  h=xa[khi]-xa[klo];
  if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
/* evaluate the spline */
  a=(xa[khi]-x)/h;  b=(x-xa[klo])/h;
/*  *y = a*ya[klo] + b*ya[khi]; */
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;

/* evaluate the derivative */
  ap = -1.0/h; bp = 1.0/h;
/*  *yp = ap*ya[klo] + bp*ya[khi]; */
  *yp = ap*ya[klo] + bp*ya[khi] + 
    ( ap*(3.0*a*a - 1.0)*y2a[klo] + bp*(3.0*b*b - 1.0)*y2a[khi] )*(h*h)/6.0; 

/* second derivative */
  *ypp = a*y2a[klo] + b*y2a[khi]; 

}
