typedef struct{
	char	*name;
	int		npoints;
	real	*x;
	real	*y;
}  pointlist;

typedef struct {
	char	*name;
	int		n1;
	int		n2;
	complex	*z;
} complexgrid;
	
typedef struct {
	int		npoints;
	complex	*z;
} complexcurve;

typedef struct {
	char	*name;
	int		npoints;
	int		ndimensions;
	real	*ulist;
	real	*B;
} CPCspline;

typedef struct {
	char	*name;
	int	npoints;
	real	*a;
	real	*b;
	real	te;
	real	beta;
	complex	z0;
	complex	z1;
	complex	c;
	real	rte;
	real	dtau;
	real	dxte;
	real	dyte;
} tgcoefflist;

typedef struct {
	char	*name;
	int		n1;
	int		n2;
	real	*p;
	complex	*u;
	real	alfa;
	real	t;
} complexflow;

typedef struct {
	char	*name;
	int		n1;
	int		n2;
	double	*fi;
} potential;

typedef struct {
	char	*name;
	int		n1;
	int		n2;
	int		dim;
	double	*x;
} v1grid;

typedef struct tensortype{
	char	*name;
	char	*comp;
	int		n1;
	int		n2;
	int		dim;
	int		rank;
	double	*tc;
} tensor;

typedef struct {
	char	*name;
	int		n1;
	int		n2;
	int		dim;
	int		*pmask;
	int		*umask;
} v1mask;

typedef struct {
	char	*name;
	int		n1;
	int		n2;
	int		*msk;
} mask;

typedef struct {
	char	*name;
	char	*comp;
	int		n1;
	int		n2;
	int		dim;
	double	*p;
	double	*u;
	double	re;
	double	t;
	double	dt;
	int 	it;
} v1flow;


typedef struct {
	char	*gridname;
	char	*flowname;
	char	*maskname;
	int		imax;
	double	et;
	double	gt;
} v1param;

/* Function prototypes */
void 
mktgcoeff(int npoints, real *points, tgcoefflist *tgcoeff);
void 
mktggrid(real *x, real *y, real r, real s, tgcoefflist tgcoeff);
pointlist 
aspointlist(char *name, int npoints, real *x, real *y);
real *
CPCp(real *pm,real u,CPCspline spline);
real 
F1(real u);
real 
F2(real u);
real 
F3(real u);
real 
F4(real u);
real 
length(real x,real y);
CPCspline 
mkCPCspline(pointlist points);
void 
reparametrizeCPCspline(CPCspline *spline,real *parlist);
real 
tan_angle(real x1,real y1,real x2,real y2);
