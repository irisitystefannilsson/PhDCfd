// -*-c++-*-
// We first define some functions
// used to implement boundary- and
// initial conditions for the gridFunctions

// This is used to set the correct
// boundary conditions for the
// pressure equation
namespace ins
{
  static int grid = 0;
  static int side = 0;
  enum drivers {constantInflow, parabolicInflow, drivenLid};

  void
  setPressureBCs(gridFunction & RHS, 
                 gridFunction & pressure, 
                 gridFunction & U, 
                 gridFunction & V, 
                 CompositeGrid *Kimera);

  // rhs_px is used in setPressureBCs...
  doubleArray
  rhs_px(Index ix, 
         Index iy, 
         int k, 
         double time, 
         gridFunction & U, 
         gridFunction & V);
  
  // ...as is rhs_py
  doubleArray 
  rhs_py(Index ix, 
         Index iy, 
         int k, 
         double time, 
         gridFunction & U, 
         gridFunction & V);
  
  // Cd is the multiplying factor
  // used for the divergence-sink in
  // the pressure equation
  doubleArray 
  Cd(int grid, 
     CompositeGrid *Kimera, 
     double dt, 
     double viscosity, 
     Index ix, 
     Index iy);
  
  // twilight-flow solution for testing of convergence
  doubleArray
  twilight_U(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid=-1, 
             int side=-1);
  doubleArray
  twilight_V(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid=-1, 
             int side=-1);
  doubleArray
  twilight_P(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid=-1, 
             int side=-1);
  
  
  // inlet_[U,V][_dt] gives
  // values for the velocity field and
  // its time-derivative and boundaries
  doubleArray
  inlet_U(const doubleArray & x, 
          const doubleArray & y, 
          const double t, 
          int grid=-1, 
          int side=-1);
  
  doubleArray
  inlet_V(const doubleArray & x, 
          const doubleArray & y, 
          const double t, 
          int grid=-1, 
          int side=-1);
  
  doubleArray
  inlet_U_dt(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid=-1, 
             int side=-1);
  
  doubleArray
  inlet_V_dt(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid=-1, 
             int side=-1);
  
  // zero returns doubleArray
  // with all zeros
  doubleArray
  zero(const doubleArray & x, 
       const doubleArray & y, 
       const double t, 
       int grid=-1, 
       int side=-1);
  
  // unity returns doubleArray
  // with all ones
  doubleArray
  unity(const doubleArray & x, 
        const doubleArray & y, 
        const double t, 
        int grid=-1, 
        int side=-1);
  
  // compute how well we resolve
  // the spatial scales of the
  // flow
  doubleArray
  smallestLengthScale(int grid, 
                      CompositeGrid *Kimera, 
                      gridFunction & U, 
                      gridFunction & V, 
                      double viscosity, 
                      Index ix, 
                      Index iy);
  
  // Time to accelerate the inflow 
  // boundary condition
  const double accelerationTime = 1.;
  
  // Constants used for non-symmetric 
  // perturbation of flow
  const double pertTime = 0.1;
  const double pert = 0.1;
  
  // Viscosity of the fluid
  double viscosityG = 1;
  
  // hdf5-File to store results in
#define RESULTFILE "NSflow.h5"
  
  //
  // Then we define the functions declared above
  //
	  
  doubleArray
  twilight_U(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid, 
             int side)
  {
    doubleArray resArray(x);
    resArray = -cos(x)*sin(y)*exp(-2*t);
    return resArray;
  }
  doubleArray
  twilight_U_dt(const doubleArray & x, 
                const doubleArray & y, 
                const double t, 
                int grid, 
                int side)
  {
    doubleArray resArray(x);
    resArray = 2*exp(-2*t)*cos(x)*sin(y);
    return resArray;
  }
  doubleArray
  twilight_V(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid, 
             int side)
  {
    doubleArray resArray(x);
    resArray = sin(x)*cos(y)*exp(-2*t);
    return resArray;
  }
  doubleArray
  twilight_V_dt(const doubleArray & x, 
                const doubleArray & y, 
                const double t, 
                int grid, 
                int side)
  {
    doubleArray resArray(x);
    resArray = -2*exp(-2*t)*sin(x)*cos(y);
    return resArray;
  }
  doubleArray
  twilight_P(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid, 
             int side)
  {
    doubleArray resArray(x);
    resArray = -0.25*(cos(2*x)+cos(2*y))*exp(-4*t);
    return resArray;
  }
  doubleArray
  twilight_P_x(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid, 
             int side)
  {
    doubleArray resArray(x);
    resArray = 5.0E-1*exp(-4*t)*sin(2*x);
    return resArray;
  }
  doubleArray
  twilight_P_y(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid, 
             int side)
  {
    doubleArray resArray(x);
    resArray = 5.0E-1*exp(-4*t)*sin(2*y);
    return resArray;
  }
  doubleArray
  inlet_U(const doubleArray & x, 
          const doubleArray & y, 
          const double t, 
          int grid, 
          int side)
  {
    if (grid==ins::grid && side==ins::side)
      {
        if (t <= 1.)
          {
            doubleArray retA(y);
	    retA = 
              30*(pow(t,5.)/5. + pow(t,3.)/3. - pow(t,4.)/2)*cos(y*0);
            return retA;
          }
        else if (t > 1. && t <= 1 + pertTime)
          {
            double tt = (t - 1.)/pertTime;
            doubleArray retB = 
              1.*cos(y*0) + pert*30*(pow(tt,5.)/5. + pow(tt,3.)/3. - pow(tt,4.)/2) *cos(y);
            return retB;
            
          }
        else if (t > 1. + pertTime && t <= 1 + 2.*pertTime)
          {
            double tt = (t - 1. - pertTime)/pertTime;
            doubleArray retC = 
              (1. + pert)*cos(y*0) - pert*30*(pow(tt,5.)/5. + pow(tt,3.)/3. - pow(tt,4.)/2) *cos(y);
            return retC;
            
          }
        else
          {
            doubleArray retD = cos(y*0);
            return retD;
          }
      }
    
    doubleArray retB = sin(y*0);
    return retB;
  }
  
  doubleArray
  inlet_V(const doubleArray & x, 
          const doubleArray & y, 
          const double t, 
          int grid, 
          int side)
  {
    doubleArray retA = sin(x*0);
    
    return retA;
  }
  
  doubleArray
  inlet_U_dt(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid, 
             int side)
  {
    if (grid == ins::grid && side == ins::side && t < accelerationTime)
      {
        doubleArray retA = 
          30*pow(-1.+t,2.)*pow(t,2.) *cos(y*0);
        return retA;
      }
  
    doubleArray retB = sin(x*0);
    return retB;
  }

  doubleArray
  inlet_V_dt(const doubleArray & x, 
             const doubleArray & y, 
             const double t, 
             int grid, 
             int side)
  {
    doubleArray retA = sin(x*0);
  
    return retA;
  }

  doubleArray
  lidDriver(const doubleArray & x, 
       const doubleArray & y, 
       const double t, 
       int grid, 
       int side)
  {
    if (grid == 0 && side == 3)
      {
        if (t <= 1.)
          {
            doubleArray retA = 
              30*(pow(t,5.)/5. + pow(t,3.)/3. - pow(t,4.)/2) *cos(y*0);
            //            std::cout << "In lid driver now!!" << std::endl;
            return retA;
          }
        doubleArray retA = cos(0*y);
        return retA;
      }
    doubleArray retA = sin(x*0);
    
    return retA;
  }

  doubleArray
  zero(const doubleArray & x, 
       const doubleArray & y, 
       const double t, 
       int grid, 
       int side)
  {
    doubleArray retA = sin(x*0);

    return retA;
  }

  doubleArray
  unity(const doubleArray & x, 
        const doubleArray & y, 
        const double t, 
        int grid, 
        int side)
  {
    doubleArray retA = cos(x*0);
  
    return retA;
  }

  void
  setPressureBCs(gridFunction & RHS, 
                 gridFunction & pressure, 
                 gridFunction & U, 
                 gridFunction & V, 
                 CompositeGrid *Kimera)
  {
    int k;

    double dt = U.getTime()-pressure.getTime();

    for (k=0; k<Kimera->nrGrids(); k++)
      {
        Index ixb = pressure.getDiscrBounds(0,k);
        Index low_r(pressure.getDiscrBounds(0,k).getBase()), hi_r(pressure.getDiscrBounds(0,k).getBound());
        Index iyb = pressure.getDiscrBounds(1,k);
        Index low_s(pressure.getDiscrBounds(1,k).getBase()), hi_s(pressure.getDiscrBounds(1,k).getBound()); 
      
        for (int side=Side(lowR); side<=Side(hiS); side++)
          {
            if (pressure.getBoundaryType(k,(Side) side) == INTERPOLATION) // interpolation
              {
              }
            else if (pressure.getBoundaryType(k,(Side) side) == DIRICHLET) // Dirichlet
              {
              }
            else if (pressure.getBoundaryType(k,(Side) side) == NEUMANN) // Neumann
              {
                if (side==0)
                  RHS(k)(low_r-1,iyb) = Kimera->normalVector_R_Component(k,side,ixb,low_r,hi_r,iyb,low_s,hi_s)*rhs_px(low_r,iyb,k,dt,U,V) + Kimera->normalVector_S_Component(k,side,ixb,low_r,hi_r,iyb,low_s,hi_s)*rhs_py(low_r,iyb,k,dt,U,V);
                else if (side==1)
                  RHS(k)(hi_r+1,iyb) = Kimera->normalVector_R_Component(k,side,ixb,low_r,hi_r,iyb,low_s,hi_s)*rhs_px(hi_r,iyb,k,dt,U,V) + Kimera->normalVector_S_Component(k,side,ixb,low_r,hi_r,iyb,low_s,hi_s)*rhs_py(hi_r,iyb,k,dt,U,V);
                else if (side==2)
                  RHS(k)(ixb,low_s-1) = Kimera->normalVector_R_Component(k,side,ixb,low_r,hi_r,iyb,low_s,hi_s)*rhs_px(ixb,low_s,k,dt,U,V) + Kimera->normalVector_S_Component(k,side,ixb,low_r,hi_r,iyb,low_s,hi_s)*rhs_py(ixb,low_s,k,dt,U,V);
                else if (side==3)
                  RHS(k)(ixb,hi_s+1) = Kimera->normalVector_R_Component(k,side,ixb,low_r,hi_r,iyb,low_s,hi_s)*rhs_px(ixb,hi_s,k,dt,U,V) + Kimera->normalVector_S_Component(k,side,ixb,low_r,hi_r,iyb,low_s,hi_s)*rhs_py(ixb,hi_s,k,dt,U,V);
              }
            else if (pressure.getBoundaryType(k,(Side) side) == PERIODIC) // Periodic
              {
              }
          }
      }
    for (k=0; k<Kimera->nrGrids(); k++)
      {
        Index ixb = Kimera->getBounds(0,k);
        Index low_r(pressure.getDiscrBounds(0,k).getBase()), hi_r(pressure.getDiscrBounds(0,k).getBound());
        Index iyb = Kimera->getBounds(1,k);
        Index low_s(pressure.getDiscrBounds(1,k).getBase()), hi_s(pressure.getDiscrBounds(1,k).getBound()); 
      
        for (int side=Side(lowR); side<=Side(hiS); side++)
          {
            if (pressure.getBoundaryType(k,(Side) side) == INTERPOLATION) // interpolation
              {
              }
            else if (pressure.getBoundaryType(k,(Side) side) == DIRICHLET) // Dirichlet
              {
                pressure.dirichlet(side,k,RHS);
              }
            else if (pressure.getBoundaryType(k,(Side) side) == NEUMANN) // Neumann
              {
              }
          }
      }
  }

  doubleArray
  rhs_px(Index ix, 
         Index iy, 
         int k, 
         double dt, 
         gridFunction & U, 
         gridFunction & V)
  {
    doubleArray resArray = 
      //      U.forcing(& twilight_P_x, k, 0,ix,iy);
      - U.boundary_d_dt(k, 0, ix, iy)
      - U(k)(ix,iy)*U.x(k,ix,iy) - V(k)(ix,iy)*U.y(k,ix,iy) 
      - viscosityG*(V.xy(k,ix,iy) - U.yy(k,ix,iy)) 
      ;

    return resArray;
  }
  
  doubleArray 
  rhs_py(Index ix, 
         Index iy, 
         int k, 
         double dt, 
         gridFunction & U, 
         gridFunction & V)
  {
    doubleArray resArray = 
      //      U.forcing(& twilight_P_y, k, 0,ix,iy)
      - V.boundary_d_dt(k, 0, ix, iy) 
      - U(k)(ix,iy)*V.x(k,ix,iy) - V(k)(ix,iy)*V.y(k,ix,iy) 
      - viscosityG*(U.xy(k,ix,iy) - V.xx(k,ix,iy)) 
      ;

    return resArray;
  }

  doubleArray 
  Cd(int grid, 
     CompositeGrid *Kimera, 
     double dt, 
     double viscosity, 
     Index ix, 
     Index iy)
  {
    double Ct = 0.25, C0 = 1.;
    if (Kimera->gType(grid) == 1)
      {
        doubleArray result1 = 
          C0*viscosity*(
                        + 1. / 
                        pow(Kimera->gridSize(grid, 0, 0, ix, iy), 2.) 
                        + 1. / 
                        pow(Kimera->gridSize(grid, 1, 0, ix, iy), 2.)
                        );
        result1 = min(result1, Ct / dt);
        return result1;
      }
    else
      {
        doubleArray result2 = 
          C0*viscosity*(
                        + 1. / (
                                + pow(Kimera->gridSize(grid, 0, 0, ix, iy), 2.) 
                                + pow(Kimera->gridSize(grid, 1, 1, ix, iy), 2.)
                                )
                        + 1. / (
                                + pow(Kimera->gridSize(grid, 0, 1, ix, iy), 2.) 
                                + pow(Kimera->gridSize(grid, 1, 0, ix, iy), 2.) 
                                )
                        );
        result2 = min(result2, Ct / dt);
        return result2;
      }
  }

  doubleArray
  smallestLengthScale(int grid, 
                      CompositeGrid *Kimera, 
                      gridFunction & U, 
                      gridFunction & V, 
                      double viscosity, 
                      Index ix, 
                      Index iy)
  {
    double small = 0.08;
    doubleArray result(ix,iy);

    result = sqrt(viscosity / max(pow(U.x(grid,ix,iy),2.)+pow(U.y(grid,ix,iy),2.),pow(V.x(grid,ix,iy),2.)+pow(V.y(grid,ix,iy),2.), small));

    if (Kimera->gType(grid) == 1)
      result /= max(Kimera->gridSize(grid, 0, 0, ix, iy),Kimera->gridSize(grid, 1, 0, ix, iy));
    else
      result /= max(max(Kimera->gridSize(grid, 0, 0, ix, iy),Kimera->gridSize(grid, 0, 1, ix, iy)),max(Kimera->gridSize(grid, 1, 0, ix, iy),Kimera->gridSize(grid, 1, 1, ix, iy)) );

    return result;
  }

}
