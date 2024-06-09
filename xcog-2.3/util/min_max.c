#include "real.h"
#include "min_max.h"

int int_min(int a, int b){
  if (a < b)
    return a;
  else
    return b;
}

int int_max(int a, int b){
  if (a > b)
    return a;
  else
    return b;
}

real real_min(real a, real b){
  if (a < b)
    return a;
  else
    return b;
}

real real_max(real a, real b){
  if (a > b)
    return a;
  else
    return b;
}
