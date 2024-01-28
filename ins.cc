#include <mpi.h>
#include "GridFunc.hh"
#include "OGEquation.hh"


#include <algorithm>

#define RESULTFILE "twilight.h5"

#include "my_auxiliary_routines.hh"

void
setPressureBCs(gridFunction & RHS, 
	       gridFunction & pressure, 
	       gridFunction & U, 
	       gridFunction & V, 
	       CompositeGrid *Kimera);

doubleArray
rhs_px(Index ix, 
       Index iy, 
       int k, 
       double time, 
       gridFunction & U, 
       gridFunction & V);

doubleArray 
rhs_py(Index ix, 
       Index iy, 
       int k, 
       double time, 
       gridFunction & U, 
       gridFunction & V);

doubleArray 
Cd(int grid, 
   CompositeGrid *Kimera, 
   double dt, 
   double viscosity, 
   Index ix, 
   Index iy);

doubleArray
smallestLengthScale(int grid, 
		    CompositeGrid *Kimera, 
		    gridFunction & U, 
		    gridFunction & V, 
		    double viscosity, 
		    Index ix, 
		    Index iy);


/*
  main program starts here
*/

int 
main(int argc, char** argv)
{
  ios::sync_with_stdio();
  
  int Number_of_Processors = 0;
  PetscBool flg;
  
  Optimization_Manager::Initialize_Virtual_Machine("", Number_of_Processors, argc, argv);

  // initialize PETSc
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  Partitioning_Type::SpecifyDefaultInternalGhostBoundaryWidths(1, 1);

  int myid = Communication_Manager::My_Process_Number;
  
  MPI_Comm myComm;
  MPI_Group myGroup;
  
  MPI_Comm_group(MPI_COMM_WORLD, &myGroup);
  MPI_Comm_create(MPI_COMM_WORLD, myGroup, &myComm);
  
  Index::setBoundsCheck(on);

  /*
    Kimera is a composite grid built from P++ arrays
  */
  int dim[2], dimsize = 2;
  CompositeGrid *Kimera = 0x0;
  PetscOptionsGetIntArray(PETSC_NULL, PETSC_NULL, "--square", dim, &dimsize, &flg);
  if (flg == PETSC_TRUE)
    Kimera = new CompositeGrid(dim[0], dim[1]);
  else
    {
      char filename[100];
      PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "--filename", filename, 100, &flg);
      if (flg == PETSC_TRUE)
	Kimera = new CompositeGrid(filename);
      else
	{
	  cout << "Hey!, we need a grid from somewhere \n\n";
	  exit(1);
	}
    }
  double starttime = 0.0;
  bool restart = false;
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "--restart", &starttime, &flg);
  if (flg == PETSC_TRUE) restart = true;

  /*
    Declare some gridfunctions attached to the Kimera grid
  */
  gridFunction U((*Kimera), starttime), V((*Kimera), starttime);
  gridFunction oldU((*Kimera), starttime), oldV((*Kimera), starttime);
  gridFunction tempU((*Kimera), starttime), tempV((*Kimera), starttime);
  gridFunction pressure((*Kimera), starttime), oldpressure((*Kimera), starttime);
  gridFunction RHS((*Kimera), starttime), div((*Kimera), starttime);
  
  /*
    Equations that need to be implicitly solved during the simulation
  */
  sparseSolver sps = Aztec;

  OGEquation pressurePoisson(&(*Kimera), sps);
  OGEquation U_momentum(&(*Kimera), PETSc), V_momentum(&(*Kimera), PETSc);
  

  //  U.setUseConservativeDifferenceOperators();
  //  V.setUseConservativeDifferenceOperators();

  /*
    Initialize U, V, and Pressure and their boundary value functions
  */
  U.setBC(xComponentOfVector); V.setBC(yComponentOfVector);
  oldU.setBC(xComponentOfVector); oldV.setBC(yComponentOfVector);
  tempU.setBC(xComponentOfVector); tempV.setBC(yComponentOfVector);

  functionType pressureType = scalar;

  PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "--allDir", &flg);
  if (flg == PETSC_TRUE) pressureType = DirichletAllOver;

  RHS.setBC(scalar); pressure.setBC(pressureType); div.setBC(DirichletAllOver);

  
  pressure.setFlagValues();
 
  U.setFlagValues(); V.setFlagValues(); 
  div.setFlagValues();
  
  pressurePoisson.setUnknown(pressure); 
  U_momentum.setUnknown(U); V_momentum.setUnknown(V);
  
  bool allNeumann = false;
  if (pressureType != DirichletAllOver)
    {
      PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "--allNeumann", &flg);
      if (flg == PETSC_TRUE) allNeumann = true;
    }

  U_momentum.setOperator(Poisson_with_TD_diagonal); 
  V_momentum.setOperator(Poisson_with_TD_diagonal); 
  pressurePoisson.setOperator(Poisson, allNeumann);

  U.set_Dirichlet(& U_b); U.set_bc_x(& zero); U.set_bc_y(& zero);
  V.set_Dirichlet(& V_b); V.set_bc_x(& zero); V.set_bc_y(& zero);
  oldU.set_Dirichlet(& U_b); oldU.set_bc_x(& zero); oldU.set_bc_y(& zero);
  oldV.set_Dirichlet(& V_b); oldV.set_bc_x(& zero); oldV.set_bc_y(& zero);
  pressure.set_Dirichlet(& P_b); pressure.set_bc_x(& zero); pressure.set_bc_y(& zero);
  RHS.set_Dirichlet(& P_b);

  if (restart)
    {
      U.readFromHDF5File(RESULTFILE, "uAn");
      V.readFromHDF5File(RESULTFILE, "vAn");
    }
  else
    {
      U.initializeGridData();
      V.initializeGridData();
    }
  pressure.initializeGridData();
  oldU.initializeGridData();
  oldV.initializeGridData();
  RHS.initializeGridData();
  //  div.initializeGridData();

  int maxNits = 100;
  int printIntervall = 20;
  double cfl = 0.5;
  
  double stoptime = 0.1, startTime = MPI_Wtime();
  
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "--stoptime", &stoptime, &flg);

  PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "--maxNits", &maxNits, &flg);

  PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "--savematrix", &flg);
  if (flg == PETSC_TRUE)
    {
      pressurePoisson.writeMatrixOperatorToFile("laplaceOp");
      U_momentum.writeMatrixOperatorToFile("umom");
      V_momentum.writeMatrixOperatorToFile("vmom");
    }
  
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-cfl", &cfl, &flg);
  
  /*
    Now we integrate the INS-equation in time
  */
  

  if (myid == 0) cout << "\nStarting time-stepping now\n";
  
  //  div.initializeGridData();

  double st = 0, dt = 0, oldDt = 0;
  bool noBC = false;
  int k;
  
  //  pressurePoisson.setOption(ksp_method, Gmres);
  //  pressurePoisson.setOption(precond_method, Asm);
  //  pressurePoisson.setOption(sub_block_method, Ilu);


  dt = U.get_CFL_Dt(cfl, U, V, (*Kimera));
  //  dt = U.getDt(cfl, U, V, viscosity*0.,(*Kimera));
  oldpressure.addTime(-dt);
  /*
    Start time-stepping
  */

  for (int it=0; it<maxNits; it++)
    {
      oldDt = dt;
      //      dt = U.get_CFL_Dt(cfl, U, V, (*Kimera));
      dt = min(dt, stoptime-U.getTime());
      
      if (myid == 0)
	cout << "Dt = " << dt << endl;
      if (it % printIntervall == 0 && myid ==0)
	{
	  cout << "\nIteration nmbr " << it << ", time is " << U.getTime() << endl;
	  cout << " : dt = " << dt << endl;
	}
      
      tempU = U;
      tempV = V;
      oldpressure = pressure;

      //      U.setDivergenceToZeroAtBoundaries(V, xComponentOfVector);
      //      V.setDivergenceToZeroAtBoundaries(U, yComponentOfVector);
      
      for (k=0; k<(*Kimera).nrGrids(); k++)
	{
	  Index ix = pressure.getDiscrBounds(0,k);
	  Index iy = pressure.getDiscrBounds(1,k);
	  
	  div(k)(ix,iy) = U.x(k,ix,iy)+V.y(k,ix,iy);
	}
      div.finishBCs();

      double divL_h = div.l_h();
      if (it % 4 == 0 && myid ==0)
	{
	  cout << " : l_h norm of divergence is " << divL_h << endl;
	}
      double divL_2 = div.l_2();
      if (it % 4 == 0 && myid ==0)
	{
	  cout << " : l_2 norm of divergence is " << divL_2 << endl;
	}
      double divL_inf = div.l_inf();
      if (it % 4 == 0 && myid ==0)
	{
	  cout << " : l_inf norm of divergence is " << divL_inf << endl;
	}
      if (it % printIntervall == 0 && myid ==0)
	{
	  cout << "---PREDICTOR STEP---" << endl;
	}
      
      for (k=0; k<(*Kimera).nrGrids(); k++)
	{
	  Index ix = pressure.getDiscrBounds(0,k);
	  Index iy = pressure.getDiscrBounds(1,k);
	  
	  RHS(k)(ix,iy) = 
	    //	    + Cd(k,&(*Kimera),dt,viscosity,ix,iy)*div(k)(ix,iy)
 	    - 2.*U.y(k,ix,iy)*V.x(k,ix,iy)
 	    - U.x(k,ix,iy)*U.x(k,ix,iy) - V.y(k,ix,iy)*V.y(k,ix,iy)
            + pressure.forcing(& F_u_x, k, 0)
	    + pressure.forcing(& F_v_y, k, 0)
	    ;
	}
      setPressureBCs(RHS, pressure, U, V, &(*Kimera));

      if (it % printIntervall == 0 && myid == 0)
	{
	  st = MPI_Wtime();
	  cout << " Solving for pressure now " 
	       << endl;
	}
      pressurePoisson.solve(RHS, &(*Kimera), 0, allNeumann, noBC);
      //      pressure.initializeGridData();
      //      pressurePoisson.solve(RHS, &(*Kimera), 0);
      if (it % printIntervall == 0 && myid == 0)
	cout << "  : solution time for pressure-equation was " 
	     << MPI_Wtime() - st 
	     << endl;

      pressure.addTime(dt);
      oldpressure.addTime(dt);
      if (it == 0)
	{
	  oldpressure = pressure;
	  oldU = U;
	  oldV = V;
	}

      for (k=0; k<(*Kimera).nrGrids(); k++)
	{
	  Index ui = U.getDiscrBounds(0,k);   Index vi = V.getDiscrBounds(0,k);
	  Index uj = U.getDiscrBounds(1,k);   Index vj = V.getDiscrBounds(1,k);
	  
	  U(k)(ui,uj) = U(k)(ui,uj) +
	    dt*( 
		- ( (1.+dt/(2.*oldDt))*(U(k)(ui,uj)*U.x(k) + V(k)(ui,uj)*U.y(k))
		    - (dt/(2.*oldDt))*(oldU(k)(ui,uj)*oldU.x(k) + oldV(k)(ui,uj)*oldU.y(k)) )
		+ 0.5*( viscosityG*U.laplacian(k) )
		+ U.forcing(& zero, k, 0.5*dt)
		- ( (1.+dt/(2.*oldDt))*pressure.x(k,ui,uj) - (dt/(2.*oldDt))*oldpressure.x(k,ui,uj) )
		);
	  V(k)(vi,vj) = V(k)(vi,vj) + 
	    dt*( 
		- ( (1.+dt/(2.*oldDt))*(tempU(k)(vi,vj)*V.x(k) + V(k)(vi,vj)*V.y(k))
		    - (dt/(2.*oldDt))*(oldU(k)(vi,vj)*oldV.x(k) + oldV(k)(vi,vj)*oldV.y(k)) )
		+ 0.5*( viscosityG*V.laplacian(k) )
		+ V.forcing(& zero, k, 0.5*dt)
		- ( (1.+dt/(2.*oldDt))*pressure.y(k,vi,vj) - (dt/(2.*oldDt))*oldpressure.y(k,vi,vj) )
		);
	  
	}
      U.addTime(dt);
      V.addTime(dt);
      
      if (it % printIntervall == 0 && myid == 0)
	{
	  st = MPI_Wtime();
	  cout << " Solving for U now " << endl;
	}
      U_momentum.solve(U, &(*Kimera), 0.5*viscosityG*dt);
      if (it % printIntervall == 0 && myid == 0)
	cout << "  : solution time for U-equation was " 
	     << MPI_Wtime() - st 
	     << endl;
      
      if (it % printIntervall == 0 && myid == 0)
	{
	  st = MPI_Wtime();
	  cout << " Solving for V now " << endl;
	}
      V_momentum.solve(V, &(*Kimera), 0.5*viscosityG*dt);
      if (it % printIntervall == 0 && myid == 0)
	cout << "  : solution time for V-equation was " 
	     << MPI_Wtime() - st 
	     << endl;
      if (it % printIntervall == 0 && myid ==0)
	{
	  cout << "---CORRECTOR STEP---" << endl;
	}
      oldpressure = pressure;
      
      oldU = tempU;
      oldV = tempV;
      tempU = U;
      tempV = V;

      //      U.setDivergenceToZeroAtBoundaries(V, xComponentOfVector);
      //      V.setDivergenceToZeroAtBoundaries(U, yComponentOfVector);

      for (k=0; k<(*Kimera).nrGrids(); k++)
	{
	  Index ix = pressure.getDiscrBounds(0,k);
	  Index iy = pressure.getDiscrBounds(1,k);
	  
	  div(k)(ix,iy) = U.x(k,ix,iy)+V.y(k,ix,iy);
	}
      div.finishBCs();
      for (k=0; k<(*Kimera).nrGrids(); k++)
	{
	  Index ix = pressure.getDiscrBounds(0,k);
	  Index iy = pressure.getDiscrBounds(1,k);
	  
	  RHS(k)(ix,iy) = 
	    //	    + Cd(k,&(*Kimera),dt,viscosityG,ix,iy)*div(k)(ix,iy)
 	    - 2.*U.y(k,ix,iy)*V.x(k,ix,iy)
 	    - U.x(k,ix,iy)*U.x(k,ix,iy) - V.y(k,ix,iy)*V.y(k,ix,iy)
	    + pressure.forcing(& F_u_x, k, 0)
	    + pressure.forcing(& F_v_y, k, 0)
	    ;
	}
      setPressureBCs(RHS, pressure, U, V, &(*Kimera));
      if (it % printIntervall == 0 && myid == 0)
	{
	  st = MPI_Wtime();
	  cout << " Solving for pressure now " << endl;
	}
      pressurePoisson.solve(RHS, &(*Kimera), 0, allNeumann, noBC);
      //      pressure.initializeGridData();
      //      pressurePoisson.solve(RHS, &(*Kimera), 0);
      if (it % printIntervall == 0 && myid == 0)
	cout << "  : solution time for pressure-equation was " 
	     << MPI_Wtime() - st 
	     << endl;

      for (k=0; k<(*Kimera).nrGrids(); k++)
	{
	  Index ui = U.getDiscrBounds(0,k);    Index vi = V.getDiscrBounds(0,k);
	  Index uj = U.getDiscrBounds(1,k);    Index vj = V.getDiscrBounds(1,k);
	  
	  U(k)(ui,uj) = oldU(k)(ui,uj) +
	    dt*( 
		- 0.5*( 1.*(U(k)(ui,uj)*U.x(k) + V(k)(ui,uj)*U.y(k))
			+ 1.*(oldU(k)(ui,uj)*oldU.x(k) + oldV(k)(ui,uj)*oldU.y(k)) )
		+ 0.5*( viscosityG*oldU.laplacian(k) )
		+ U.forcing(& zero, k, -0.5*dt)
		- 0.5*( 1.*pressure.x(k,ui,uj) + 1.*oldpressure.x(k,ui,uj) )
		);
	  V(k)(vi,vj) = oldV(k)(vi,vj) + 
	    dt*( 
		- 0.5*( 1.*(tempU(k)(vi,vj)*V.x(k) + tempV(k)(vi,vj)*V.y(k))
			+ 1.*(oldU(k)(vi,vj)*oldV.x(k) + oldV(k)(vi,vj)*oldV.y(k)) )
		+ 0.5*( viscosityG*oldV.laplacian(k) )
		+ V.forcing(& zero, k, -0.5*dt)
		- 0.5*( 1.*pressure.y(k,vi,vj) + 1.*oldpressure.y(k,vi,vj) )
		);
	  
	}
      if (it % printIntervall == 0 && myid == 0)
	{
	  st = MPI_Wtime();
	  cout << " Solving for U now " << endl;
	}
      U_momentum.solve(U, &(*Kimera), 0.5*viscosityG*dt);
      if (it % printIntervall == 0 && myid == 0)
	cout << "  : solution time for U-equation was " 
	     << MPI_Wtime() - st 
	     << endl;
      
      if (it % printIntervall == 0 && myid == 0)
	{
	  st = MPI_Wtime();
	  cout << " Solving for V now " << endl;
	}
      V_momentum.solve(V, &(*Kimera), 0.5*viscosityG*dt);
      if (it % printIntervall == 0 && myid == 0)
	cout << "  : solution time for V-equation was " 
	     << MPI_Wtime() - st 
	     << endl;

      // check if simulation can be stopped
      if (U.getTime() >= stoptime) break;
      
      if (it % printIntervall == 0 && myid ==0)
	{
	  cout << "|----------------------------------------------------------|" 
	       << endl 
	       << endl;
	}
    }

  //  V.setDivergenceToZeroAtBoundaries(U, yComponentOfVector);
  //  U.setDivergenceToZeroAtBoundaries(V, xComponentOfVector);

  for (k=0; k<(*Kimera).nrGrids(); k++)
    {
      Index ix = pressure.getDiscrBounds(0,k);
      Index iy = pressure.getDiscrBounds(1,k);
      
      div(k)(ix,iy) = U.x(k,ix,iy)+V.y(k,ix,iy);
    }  
  div.finishBCs();
  div.interpolate();

  for (k=0; k<(*Kimera).nrGrids(); k++)
    {
      Index ix = pressure.getDiscrBounds(0,k);
      Index iy = pressure.getDiscrBounds(1,k);
      
      RHS(k)(ix,iy) = 
  	- 2.*U.y(k,ix,iy)*V.x(k,ix,iy)
  	- U.x(k,ix,iy)*U.x(k,ix,iy) - V.y(k,ix,iy)*V.y(k,ix,iy)
	+ pressure.forcing(& F_u_x, k, 0)
	+ pressure.forcing(& F_v_y, k, 0)
	;
    }
  setPressureBCs(RHS, pressure, U, V, &(*Kimera));
  pressurePoisson.solve(RHS, &(*Kimera), 0, allNeumann, noBC);
  //  pressure.initializeGridData();
  //  cout << "ALL neumann? " << allNeumann << endl;
  //  pressurePoisson.solve(RHS, &(*Kimera), 0);
  /*****************Time integration is finished here*********************/
  
  U.saveToFile("Uins");
  V.saveToFile("Vins");
  div.saveToFile("div");
  pressure.saveToFile("Pins");
  Kimera->saveCoordinatesToFile();
  
  V.createHDF5File(RESULTFILE);
  (*Kimera).saveCoordinatesToHDF5File(RESULTFILE);

  U.saveToHDF5File(RESULTFILE, "uAn");
  V.saveToHDF5File(RESULTFILE, "vAn");
  pressure.saveToHDF5File(RESULTFILE, "pAn");

  double compTime = MPI_Wtime() - startTime;
  double maxTime;
  MPI_Reduce(&compTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, myComm);
  
  if (myid == 0)
    {
      cout << "\nTime-stepping finished\n";
      cout << "  Simulated time is: " << U.getTime() << endl;
      cout << "  Wall clock execution time was " << maxTime << "s\n";
    }
  
  if (myid == 0)
    cout << "Checking U-solution :\n";

  U.checkErrors(& U_b);
  
  if (myid == 0)
    cout << "Checking V-solution :\n";
  
  V.checkErrors(& V_b);
  
  if (myid == 0)
    cout << "Checking P-solution :\n";

  pressure.checkErrors(& P_b);

  delete Kimera;

  // Exit PETSc
  PetscFinalize();

  Optimization_Manager::Exit_Virtual_Machine();
  
  return 0;
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
	      //	      cout << "Setting p-dir, (grid,side) = (" << k << "," << side << ")" << endl;
	      pressure.dirichlet(side, k, RHS);
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
    //    U.forcing(& p_x, k, 0,ix,iy);
    - U.forcing(& dudt, k, dt, ix, iy)
    - U(k)(ix,iy)*U.x(k,ix,iy) - V(k)(ix,iy)*U.y(k,ix,iy) 
    - viscosityG*(V.xy(k,ix,iy) - U.yy(k,ix,iy)) 
    + U.forcing(& zero, k, dt, ix, iy)
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
    //    V.forcing(& p_y, k, 0, ix,iy);
    - V.forcing(& dvdt, k, dt, ix, iy) 
    - U(k)(ix,iy)*V.x(k,ix,iy) - V(k)(ix,iy)*V.y(k,ix,iy) 
    - viscosityG*(U.xy(k,ix,iy) - V.xx(k,ix,iy)) 
    + V.forcing(& zero, k, dt, ix, iy)
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

