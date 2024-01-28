#include <mpi.h>

#include "my_auxiliary_routines.hh"

#include <algorithm>


doubleArray
inlet_U(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  //  doubleArray retA = //sign(0.25-abs(y-0.25), sqrt(abs(0.25-abs(y-0.25))));
    //    sign(y,sqrt(abs(y)));
    // return retA;
  if (grid==0 && side==0)
    {

      doubleArray retA = 
	min(1., 30*(pow(t,5.)/5. + pow(t,3.)/3. - pow(t,4.)/2)) *cos(y*0);
      
      return retA;
    }
  
  doubleArray retB = sin(y*0);
  return retB;
}

doubleArray
inlet_V(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray retA = sin(x*0);
  
  return retA;
}

doubleArray
inlet_U_dt(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray retA = sin(x*0);
  
  return retA;
  /*  if (grid==0 && side==0 && t <0)
  {
    doubleArray retA = 
      30*pow(-1.+t,2.)*pow(t,2.) *cos(y*0);
    return retA;
  }

  doubleArray retB = sin(x*0);
  return retB;*/
}

doubleArray
inlet_V_dt(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray retA = sin(x*0);
  
  return retA;
}

doubleArray
zero(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray retA = sin(x*0);

  return retA;
}

doubleArray
unity(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray retA = cos(x*0);
  
  return retA;
}

doubleArray
bc_x(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = sin(x*0);
  
  return resArray;
}

doubleArray
bc_y(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = sin(x*0);
  
  return resArray;
}

doubleArray
F_u(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    ((cos(x)*(sin(x) + 2*exp(2*t)*(-1 + viscosityG)*sin(y)))/exp(4*t));
  //    sin(x*0);
  return resArray;
}

doubleArray
F_v(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    ((cos(y)*(-2*exp(2*t)*(-1 + viscosityG)*sin(x) + sin(y)))/exp(4*t));
  //    sin(x*0);
  return resArray;
}

doubleArray
F_u_x(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    //  2.*sin(x)*sin(y);
    sin(x*0);
  return resArray;
}

doubleArray
F_v_y(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    //    -2.*sin(y)*sin(x);
    sin(x*0);
  return resArray;
}

doubleArray
U_b(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  //const double overlap = 0.25;
  //const double R = 0.5;
  doubleArray resArray(x);
  
  resArray = -cos(x)*sin(y)*exp(-2*t);
  //  doubleArray radius = (sqrt(pow(x,2.)+pow(y,2.))-R)/overlap;
  //  if (grid == 0)
  //    where (radius <= 0)
  //      resArray = 0;
  //  else if (grid == 1)
  //    where (radius >= 1)
  //      resArray = 0;
    //    -cos(x)*sin(y);
    //    cos(x*0)*t;
    /*x*x*x-x*x+x;*/
  return resArray;
}

doubleArray
V_b(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  //const double overlap = 0.25;
  //const double R = 0.5;
  doubleArray resArray = 
    sin(x)*cos(y)*exp(-2*t);
  //  doubleArray radius = (sqrt(pow(x,2.)+pow(y,2.))-R)/overlap;
  //  if (grid == 0)
  //    where (radius <= 0)
  //      resArray = 0;
  //  else if (grid == 1)
  //    where (radius >= 1)
  //      resArray = 0;

  return resArray;
}

doubleArray
P_b(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    -0.25*(cos(2*x)+cos(2*y))*exp(-4*t);
    //    -0.25*(cos(2*x)+cos(2*y));
    //    -((-2. + t)*t*pow(x,2.))/2. - ((2. - 2.*t + pow(t,2.))*pow(y,2.))/2.;
    //    sin(x*0);
  return resArray;
}

doubleArray
dudx(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    sin(x)*sin(y)*exp(-2*t);
  return resArray;
}

doubleArray
dudy(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    -cos(x)*cos(y)*exp(-2*t);
  return resArray;
}

doubleArray
dvdx(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    cos(x)*cos(y)*exp(-2*t);
  return resArray;
}

doubleArray
dvdy(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    -sin(x)*sin(y)*exp(-2*t);
  return resArray;
}

doubleArray
dudt(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    //      -(sin(t)*pow(sin(x),2.)*sin(y));
    2.*cos(x)*sin(y)*exp(-2*t);
  //    sin(x*0);
  return resArray;
}

doubleSerialArray
dudt(const doubleSerialArray & x, const doubleSerialArray & y, const double t, int grid, int side)
{
  doubleSerialArray resArray = 
    //      -(sin(t)*pow(sin(x),2.)*sin(y));
    2.*cos(x)*sin(y)*exp(-2*t);
  //    sin(x*0);
  return resArray;
}
doubleArray
dvdt(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    -2.*cos(y)*sin(x)*exp(-2*t);
  //    sin(x*0);
  return resArray;
}
doubleSerialArray
dvdt(const doubleSerialArray & x, const doubleSerialArray & y, const double t, int grid, int side)
{
  doubleSerialArray resArray = 
    -2.*cos(y)*sin(x)*exp(-2*t);
  //    sin(x*0);
  return resArray;
}

doubleArray
laplaceP(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    -8.*viscosityG*cos(t)*cos(x)*cos(y)*sin(x)*sin(y);

  return resArray;
}

doubleArray
p_x(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    viscosityG*cos(t)*pow(cos(x),2.)*cos(y)*sin(y) - 
    viscosityG*cos(t)*cos(y)*pow(sin(x),2.)*sin(y);
  return resArray;
}

doubleArray
p_y(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    viscosityG*cos(t)*cos(x)*pow(cos(y),2.)*sin(x) - 
    viscosityG*cos(t)*cos(x)*sin(x)*pow(sin(y),2.);
  return resArray;
}

doubleArray
sinx(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    sin(x);

  return resArray;
}

doubleArray
sinxsiny(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    sin(x)*sin(y);

  return resArray;
}

doubleArray
cosxsiny(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    cos(x)*sin(y);

  return resArray;
}

doubleArray
sinxcosy(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    sin(x)*cos(y);

  return resArray;
}

doubleArray
cosxcosy(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    cos(x)*cos(y);

  return resArray;
}

doubleArray
cosx(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    cos(x);

  return resArray;
}

doubleArray
siny(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    sin(y);

  return resArray;
}

doubleArray
cosy(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = 
    cos(y);

  return resArray;
}

doubleArray
liny(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = y;

  return resArray;
}

doubleArray
linx(const doubleArray & x, const doubleArray & y, const double t, int grid, int side)
{
  doubleArray resArray = x;

  return resArray;
}

int
ArgIsNumber(char* arg, int argc, char** argv)
{
  if (argc > 1)
    {
      for (int i=1; i<argc; i++)
	{
	  if (!strcmp(arg, argv[i]))
	    return i;
	}
    }

  return 0;
}
	  
