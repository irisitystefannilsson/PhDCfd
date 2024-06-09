#include "math.h"
#include "real.h"
#include "log_cosh.h"

real log_cosh(real slope, real r0, real r1, real sharpness, real r){
  real v, lc_0, lc_1;
  const real big_c = 20.0;

/* this function needs to be protected against overflow */
  if (fabs(sharpness*(r-r0)) < big_c)
    lc_0 = log(cosh(sharpness*(r-r0)));
  else
    lc_0 = fabs(sharpness*(r-r0)) - log(2.0);

  if (fabs(sharpness*(r-r1)) < big_c)
    lc_1 = log(cosh(sharpness*(r-r1)));
  else
    lc_1 = fabs(sharpness*(r-r1)) - log(2.0);

  v = 0.5*slope/sharpness * (lc_0 - lc_1);
  return v;
}
