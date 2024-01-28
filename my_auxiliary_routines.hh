// -*-c++-*-
#ifndef MY_AUXILIARY_ROUTINES_H

#define MY_AUXILIARY_ROUTINES_H

#include "A++.h"

/*
  First some global variables and function declarations
*/

const double viscosityG = 1;



doubleArray
inlet_U(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
inlet_V(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
inlet_U_dt(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
inlet_V_dt(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
zero(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
unity(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
bc_x(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
bc_y(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
F_u(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
F_v(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
F_u_x(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
F_v_y(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
U_b(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
V_b(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
P_b(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
dudx(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
dudy(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
dvdx(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
dvdy(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
dudt(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleSerialArray
dudt(const doubleSerialArray & x, const doubleSerialArray & y, const double t, int grid=-1, int side=-1);

doubleArray
dvdt(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleSerialArray
dvdt(const doubleSerialArray & x, const doubleSerialArray & y, const double t, int grid=-1, int side=-1);

doubleArray
laplaceP(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
p_x(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
p_y(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
sinx(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
sinxsiny(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
cosxsiny(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
sinxcosy(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
cosxcosy(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
cosx(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
siny(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
cosy(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
linx(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

doubleArray
liny(const doubleArray & x, const doubleArray & y, const double t, int grid=-1, int side=-1);

int
ArgIsNumber(char* arg, int argc, char** argv);

#endif // #ifndef MY_AUXILIARY_ROUTINES_H
