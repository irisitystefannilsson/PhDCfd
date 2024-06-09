typedef struct {
	real	re;
	real	im;
} complex;

complex cneg(complex z);
complex cadd(complex z,complex w);
complex cconjugate(complex z);
real carg(complex z);
real c_abs(complex z);
