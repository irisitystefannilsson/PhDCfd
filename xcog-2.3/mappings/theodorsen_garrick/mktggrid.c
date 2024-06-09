#include "mappings.h"

#define PI 3.14159265358979

#ifdef F_NO_UNDERSCORE
extern void karmanntrefftz();
extern void theodorsengarrick();
extern void translate();
#else
extern void karmanntrefftz_();
extern void theodorsengarrick_();
extern void translate_();
#endif

void mktggrid(real *x, real *y, real r, real s, tgcoefflist tgcoeff)
{	
	real		angle;
	real		radius;
	int		m;
	complex		z;

/* Make a grid-point by forward transformations */
	m = 1;
        radius = 1.+s;
        angle = r*2.*PI;
        z.re = radius*cos(angle);
        z.im = radius*sin(angle);
#ifdef F_NO_UNDERSCORE
	theodorsengarrick(
#else
	theodorsengarrick_(
#endif
		&z,
		&m,
		tgcoeff.a,
		tgcoeff.b,
		&(tgcoeff.npoints));
#ifdef F_NO_UNDERSCORE
	translate(&m,&z,&(tgcoeff.c));
#else
	translate_(&m,&z,&(tgcoeff.c));
#endif
#ifdef F_NO_UNDERSCORE
	karmanntrefftz(
#else
	karmanntrefftz_(
#endif
		&z,
		&m,
		&(tgcoeff.beta),
		&(tgcoeff.z0),
		&(tgcoeff.z1));
	*x = z.re;
	*y = z.im;
} 
