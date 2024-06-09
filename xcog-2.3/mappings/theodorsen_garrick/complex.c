#include "mappings.h"

complex cneg(complex z)
{	
	complex	zn;
	
	zn.re = -z.re;
	zn.im = -z.im;
	return(zn);
}

complex cadd(complex z,complex w)
{	
	complex	zn;
	
	zn.re = z.re + w.re;
	zn.im = z.im + w.im;
	return(zn);
}

complex cconjugate(complex z)
{	
	complex	zn;
	
	zn.re = z.re;
	zn.im = -z.im;
	return(zn);
}

real carg(complex z)
{	
	real	a;
	
	a = atan2(z.im,z.re);
	return(a);
}

real c_abs(complex z)
{	
	real	a;
	
	a = sqrt(pow(z.re,2.0)+pow(z.im,2.0));
	return(a);
}
