////////////////////////////////////////
//                                    //
// Miniature code that solves the     //
// laminar Navier-Stokes equations    //
// for incompressible flow on         //
// an Overlapping grid generated      //
// by Xcog                            //
//                                    //
////////////////////////////////////////

#include <algorithm>
#include <sstream>

#include <mpi.h>
#include "petscsystypes.h"

// OGEquation.hh contains declarations  
// for implicit equations Ax=b
#include "OGEquation.hh"
// GridFunc.hh contains
// declarations for gridFunctions
// and their methods
#include "GridFunc.hh"
// ins_sharp contains some functions
// and constants used in the simulation
#include "ins_sharp.hh"

// This is necessary because Aztec defines
// min|max macros
#undef min
#undef max

// File to save data to
#define RESULTS "INS_data"

using std::cout;
using std::endl;

////////////////////////////////////////
//                                    //
// Main program starts here           //
//                                    //
////////////////////////////////////////

int
main(int argc, char** argv)
{
  std::ios::sync_with_stdio();
  
  int Number_of_Processors = 0;
  
  // flg is used to pick out 
  // arguments from argv
  PetscBool flg;

  // Initialize P++ and MPI 
  Optimization_Manager::Initialize_Virtual_Machine("", Number_of_Processors, argc, argv);

  // initialize PETSc
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  // We use a ghostCellWidth of 1 in both 
  // dimensions to make parallel interpolation possible
  // between the component grids
  Partitioning_Type::SpecifyDefaultInternalGhostBoundaryWidths(1, 1);

  // myid is the id of this process in the 
  // parallel virtual machine (probably the
  // same as its MPI id)
  int myid = Communication_Manager::My_Process_Number;
  
  // Create a new MPI communicator for 
  // interprocess communication outside P++
  MPI_Comm myComm;
  MPI_Group myGroup;
  
  MPI_Comm_group(MPI_COMM_WORLD, &myGroup);
  MPI_Comm_create(MPI_COMM_WORLD, myGroup, &myComm);
  
  Index::setBoundsCheck(off);

  //
  // Kimera is a composite grid built from P++ arrays
  //
  // It can be read from an ASCII file...
  // ...or from an HDF5-file
  char filename[100];
  PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "--filename", filename, 100, &flg);
  if (flg != PETSC_TRUE)
  {
    std::cerr << "Hey!, we need a grid from somewhere (no filename given, exiting...)... \n\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  CompositeGrid Kimera(filename);

  //
  // Check (using PetscOptionsGet...) if some 
  // arguments are present
  // that will modify the simulation
  //
  double starttime = 0;

  // convergence check?
  bool twilightFlow = false;
  PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "--twilight", &flg);
  if (flg == PETSC_TRUE)
  {
    twilightFlow = true;
  }
  ins::drivers inflow = ins::constantInflow;
  char driverName[100];
  PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "--inflow", driverName, 100, &flg);
  if (flg == PETSC_TRUE)
  {
    if (!strcmp("parabolic", driverName))
      inflow = ins::parabolicInflow;
    else if (!strcmp("lid", driverName))
      inflow = ins::drivenLid;
    else if (!strcmp("constant", driverName))
      inflow = ins::constantInflow;
  }
  // restart decides if we are
  // restarting the simulation
  // from an old result
  bool restart = false;
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "--restart", &starttime, &flg);
  if (flg == PETSC_TRUE)
  {
    restart = true;
  }

  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "--viscosity", &ins::viscosityG, &flg);

  if (twilightFlow)
  {
    ins::viscosityG = 1;
  }

  // 
  //  Declare some gridfunctions living on the Kimera grid
  //

  // U and V are the velocity components
  gridFunction U(Kimera, starttime), V(Kimera, starttime);
  
  // old|temp[UV] are used during the numerical
  // time-stepping
  gridFunction oldU(Kimera, starttime), oldV(Kimera, starttime);
  gridFunction tempU(Kimera, starttime), tempV(Kimera, starttime);

  // (old)pressure and the divergence of (U,V)
  gridFunction pressure(Kimera, starttime), oldpressure(Kimera, starttime), div(Kimera, starttime);
  
  // gridFunction to be used as right-hand-side
  // in the pressure equation
  gridFunction RHS(Kimera, starttime);

  // reset basic mask in composite grid
  Kimera.simplifyFlagValues();

  // We want to use a semi-implicit
  // time-stepping algorithm so we
  // create Overlapping grid equations
  // to be solved for the velocity
  // components. An equation for the
  // pressure is also needed.
  //   In the declaration we also
  // decide which sparse solver should
  // be used (Aztec or PETSc)
  OGEquation pressurePoisson(&Kimera, PETSc);
  OGEquation U_momentum(&Kimera, Aztec), V_momentum(&Kimera, Aztec);

  // Set the type of boundary condition to be used
  // for the gridFunctions (to be used together
  // with information from the Xcog-file)
  U.setBC(xComponentOfVector); V.setBC(yComponentOfVector);
  oldU.setBC(xComponentOfVector); oldV.setBC(yComponentOfVector);
  tempU.setBC(xComponentOfVector); tempV.setBC(yComponentOfVector);

  RHS.setBC(scalar); 
  pressure.setBC(scalar); 
  div.setBC(DirichletAllOver);

  functionType pressureType = scalar;

  PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "--allDir", &flg);
  if (flg == PETSC_TRUE)
  {
    pressureType = DirichletAllOver;
  }
  RHS.setBC(scalar); 
  pressure.setBC(pressureType); 
  div.setBC(DirichletAllOver);
  // Set flag values for all grid points the
  // gridFunctions are
  // defined on
  pressure.setFlagValues(); 
  U.setFlagValues(); 
  V.setFlagValues();
  // Attach gridFunctions to the
  // Overlapping grid equations
  pressurePoisson.setUnknown(pressure); 
  U_momentum.setUnknown(U); 
  V_momentum.setUnknown(V);
  // allNeumann should be set to true
  // if Neumann conditions are
  // used for the pressure on all 
  // boundaries
  bool allNeumann = false;
  PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "--allNeumann", &flg);
  if (flg == PETSC_TRUE && !(pressureType == DirichletAllOver) )
  {
    allNeumann = true;
  }
  // Set type of equation system for the
  // Overlapping grid equations
  // Heat-eqations for U,V
  U_momentum.setOperator(Poisson_with_TD_diagonal);
  V_momentum.setOperator(Poisson_with_TD_diagonal);
  // Poisson equation for the pressure
  pressurePoisson.setOperator(Poisson, allNeumann);
  // Set boundary value functions
  // for the gridFunctions
  if (inflow == ins::constantInflow)
  {
    U.set_Dirichlet(&ins::inlet_U); U.set_bc_x(&ins::zero); U.set_bc_y(&ins::zero);
    U.set_Dirichlet_time_derivative(&ins::inlet_U_dt);
    V.set_Dirichlet(&ins::inlet_V); V.set_bc_x(&ins::zero); V.set_bc_y(&ins::zero);
    V.set_Dirichlet_time_derivative(&ins::inlet_V_dt);
    oldU.set_Dirichlet(&ins::inlet_U); oldU.set_bc_x(&ins::zero); oldU.set_bc_y(&ins::zero);
    oldV.set_Dirichlet(&ins::inlet_V); oldV.set_bc_x(&ins::zero); oldV.set_bc_y(&ins::zero);
    pressure.set_Dirichlet(&ins::zero); pressure.set_bc_x(&ins::zero); pressure.set_bc_y(&ins::zero);
  }
  else if (inflow == ins::drivenLid)
  {
    U.set_Dirichlet(&ins::lidDriver);     U.set_bc_x(&ins::zero);        U.set_bc_y(&ins::zero);
    V.set_Dirichlet(&ins::zero);     V.set_bc_x(&ins::zero);        V.set_bc_y(&ins::zero);
    oldU.set_Dirichlet(&ins::lidDriver); oldU.set_bc_x(&ins::zero); oldU.set_bc_y(&ins::zero);
    oldV.set_Dirichlet(&ins::zero); oldV.set_bc_x(&ins::zero); oldV.set_bc_y(&ins::zero);
    pressure.set_Dirichlet(&ins::zero); pressure.set_bc_x(&ins::zero); pressure.set_bc_y(&ins::zero);
  }
  if (twilightFlow)
  {
    U.set_Dirichlet(&ins::twilight_U); U.set_bc_x(&ins::zero); U.set_bc_y(&ins::zero);
    U.set_Dirichlet_time_derivative(&ins::twilight_U_dt);
    V.set_Dirichlet(&ins::twilight_V); V.set_bc_x(&ins::zero); V.set_bc_y(&ins::zero);
    V.set_Dirichlet_time_derivative(&ins::twilight_V_dt);
    oldU.set_Dirichlet(&ins::twilight_U); oldU.set_bc_x(&ins::zero); oldU.set_bc_y(&ins::zero);
    oldV.set_Dirichlet(&ins::twilight_V); oldV.set_bc_x(&ins::zero); oldV.set_bc_y(&ins::zero);
    pressure.set_Dirichlet(&ins::twilight_P); pressure.set_bc_x(&ins::zero); pressure.set_bc_y(&ins::zero);
  }
  // Initialize the gridFunctions .
  // How depends on if we restart or not
  if (restart)
  {
    U.readFromHDF5File(RESULTFILE, "U");
    V.readFromHDF5File(RESULTFILE, "V");
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
  div = 0;
  // Some variables used in the time-loop
  int maxNits = 1000;
  int printIntervall = 20;
  int saveIntervall = 20;
  int fileCounter = 0;
  double cfl = 0.5;
  std::stringstream fileString;
  
  double stoptime = 0.1, startTime = MPI_Wtime();
  
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "--stoptime", &stoptime, &flg);

  PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "--maxNits", &maxNits, &flg); 

  PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "--saveIntervall", &saveIntervall, &flg); 

  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "--cfl", &cfl, &flg);  

  PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "--savematrix", &flg);
  if (flg == PETSC_TRUE)
  // We can save the matrix operators
  // in matlab-files if we want to
  {
    U_momentum.writeMatrixOperatorToFile("U_OP");
    pressurePoisson.writeMatrixOperatorToFile("laplaceOp");
  }
  // Set options used for the solution of 
  // the velocity equations
  U_momentum.setOption(precond_method, Jacobi);
  V_momentum.setOption(precond_method, Jacobi);
  pressurePoisson.setOption(ksp_method, Gmres);

  double st = 0, dt = 0;
  bool noBC = false;
  int k;
  
  fileString << RESULTS << "_" << fileCounter << ".h5";
	  
  if (myid == 0)
  {
    cout << "\nSaving initial data to file " << fileString.str().c_str() << endl;
  }
  U.createHDF5File(fileString.str().c_str());
  Kimera.saveCoordinatesToHDF5File(fileString.str().c_str());
  
  U.saveToHDF5File(fileString.str().c_str(), "u-velocity");
  V.saveToHDF5File(fileString.str().c_str(), "v-velocity");
  pressure.saveToHDF5File(fileString.str().c_str(), "pressure");

  Kimera.saveCoordinatesToFile();
  //
  // Start time-stepping
  //
  // The equation simulated is 
  // U_t + u*U_x + v*U_y - nu*(U_xx + U_yy) = grad P
  // U = (u,v) (i.e. the INS-equation)

  //
  // Most operations are done on the
  // component grid level as I am not
  // confident that the operators for 
  // whole gridFunctions work as intended
  // when gridFunctions with different
  // boundary conditions are mixed in
  // an expression
  //
  if (myid == 0)
  {
    cout << "\nStarting time-stepping now\n";
  }
  for (int it = 0; it < maxNits; it++)
  {
    // Before the inlet has achieved full
    // speed we use another field to compute
    // the time step
    if (U.getTime() < ins::accelerationTime && !twilightFlow)
    {
      tempU = 0.5;
      tempV = 0.;
      dt = U.get_CFL_Dt(cfl, tempU, tempV, Kimera);
    }
    else
    {
      dt = U.get_CFL_Dt(cfl, U, V, Kimera);
    }
    if (myid == 0)
    {
      cout << "\n New time step (" << it << "), time is; " << U.getTime() << ", dt: " << dt << endl;
    }
    
    tempU = U;
    tempV = V;
    
    oldpressure = pressure;
    
    // Set div(U,V) to 0 at (most) boundaries
    U.setDivergenceToZeroAtBoundaries(V, xComponentOfVector);
    V.setDivergenceToZeroAtBoundaries(U, yComponentOfVector);
    
    // Compute div(U,V) at all points where
    // the pressure equation is solved
    for (k = 0; k < Kimera.nrGrids(); k++)
    {
      Index ix = pressure.getDiscrBounds(0, k);
      Index iy = pressure.getDiscrBounds(1, k);
      
      div(k)(ix, iy) = U.x(k, ix, iy)+V.y(k, ix, iy);
    }
    
    if (it % printIntervall == 0 && myid == 0)
    {
      cout << "|------------------------------------------------------------|"; 
      cout << "---PREDICTOR STEP---" << endl;
    }
    
    double l_hNormOfDivergence = div.l_h();
    if (myid == 0)
    {
      cout << " l_h norm of divergence is " << l_hNormOfDivergence << endl;
    }
      
    // Construct the right-hand-side
    // for the pressure equation ...
    for (k = 0; k < Kimera.nrGrids(); k++)
    {
      Index ix = pressure.getDiscrBounds(0, k);
      Index iy = pressure.getDiscrBounds(1, k);
      
      RHS(k)(ix, iy) = 
	+ ins::Cd(k, &Kimera, dt, ins::viscosityG, ix, iy)*div(k)(ix, iy)
	- 2.*U.y(k, ix, iy)*V.x(k, ix, iy)
	- U.x(k, ix, iy)*U.x(k, ix, iy) - V.y(k, ix, iy)*V.y(k, ix, iy);
    }
    // ... and assign correct boundary
    // conditions
    ins::setPressureBCs(RHS, pressure, U, V, &Kimera);
    
    if (it % printIntervall == 0 && myid == 0)
    {
      st = MPI_Wtime();
      cout << " Solving for pressure now";
    }
    // Solve the pressure equation
    pressurePoisson.solve(RHS, &Kimera, 0, allNeumann, noBC);
    if (it % printIntervall == 0 && myid == 0)
    {
      cout << ", solution time was: " << MPI_Wtime() - st << " [s]" << endl;
    }
    pressure.addTime(dt);
    if (it == 0)
    {
      oldpressure = pressure;
      oldU = U;
      oldV = V;
    }
    
    //
    // Do the predictor step
    // for the velocity equations.
    // 
    // A mixed Adams_Bashforth(2)|Crank-Nicolson
    // scheme is used
    for (k = 0; k < Kimera.nrGrids(); k++)
    {
      Index ui = U.getDiscrBounds(0, k);      Index vi = V.getDiscrBounds(0, k);
      Index uj = U.getDiscrBounds(1, k);      Index vj = V.getDiscrBounds(1, k);
      
      U(k)(ui, uj) = U(k)(ui, uj) +
	dt*( 
	    - 0.5*( + 3.*(U(k)(ui, uj)*U.x(k) + V(k)(ui, uj)*U.y(k))
		    - 1.*(oldU(k)(ui, uj)*oldU.x(k) + oldV(k)(ui, uj)*oldU.y(k)) )
	    + 0.5*( ins::viscosityG*U.laplacian(k) )
	    + U.forcing(&ins::zero, k, 0.5*dt)
	    - 0.5*( 3.*pressure.x(k, ui, uj) - 1.*oldpressure.x(k, ui, uj) )
	     );
      V(k)(vi, vj) = V(k)(vi, vj) + 
	dt*(
	    - 0.5*( + 3.*(tempU(k)(vi, vj)*V.x(k) + tempV(k)(vi, vj)*V.y(k))
		    - 1.*(oldU(k)(vi, vj)*oldV.x(k) + oldV(k)(vi, vj)*oldV.y(k)) )
	    + 0.5*( ins::viscosityG*V.laplacian(k) )
	    + V.forcing(&ins::zero, k, 0.5*dt)
	    - 0.5*( 3.*pressure.y(k, vi, vj) - 1.*oldpressure.y(k, vi, vj) )
	    );
      
    }
    // Advance the velocity one time level
    // and finish the Crank-Nicolson part
    // of the predictor step
    U.addTime(dt);
    V.addTime(dt);
    
    if (it % printIntervall == 0 && myid == 0)
    {
      st = MPI_Wtime();
      cout << " Solving for U now";
    }
    U_momentum.solve(U, &Kimera, 0.5*ins::viscosityG*dt);

    if (it % printIntervall == 0 && myid == 0)
    {
      cout << ", solution time was: " << MPI_Wtime() - st << " [s]" << endl;
    }
    if (it % printIntervall == 0 && myid == 0)
    {
      st = MPI_Wtime();
      cout << " Solving for V now";
    }
    V_momentum.solve(V, &Kimera, 0.5*ins::viscosityG*dt);
    
    if (it % printIntervall == 0 && myid == 0)
    {
      cout << ", solution time was: " << MPI_Wtime() - st << " [s]" << endl;
    }
    //
    // The predictor step is finished here
    //
    
    if (it % printIntervall == 0 && myid == 0)
    {
      cout << "---CORRECTOR STEP---" << endl;
    }
    oldpressure = pressure;
    
    oldU = tempU;
    oldV = tempV;
    
    tempU = U;
    tempV = V;
    
    U.setDivergenceToZeroAtBoundaries(V, xComponentOfVector);
    V.setDivergenceToZeroAtBoundaries(U, yComponentOfVector);

    // In the corrector step we first
    // use the predicted velocity field
    // to compute a corrected pressure field
    for (k = 0; k < Kimera.nrGrids(); k++)
    {
      Index ix = pressure.getDiscrBounds(0, k);
      Index iy = pressure.getDiscrBounds(1, k);
      
      div(k)(ix, iy) = U.x(k, ix, iy)+V.y(k, ix, iy);
    }
    div.finishBCs();
    for (k = 0; k < Kimera.nrGrids(); k++)
    {
      Index ix = pressure.getDiscrBounds(0, k);
      Index iy = pressure.getDiscrBounds(1, k);
      
      RHS(k)(ix, iy) = 
	+ ins::Cd(k, &Kimera, dt, ins::viscosityG, ix, iy)*div(k)(ix, iy)
	- 2.*U.y(k, ix, iy)*V.x(k, ix, iy)
	- U.x(k, ix, iy)*U.x(k, ix, iy) - V.y(k, ix, iy)*V.y(k, ix, iy);
    }
    ins::setPressureBCs(RHS, pressure, U, V, &Kimera);
    if (it % printIntervall == 0 && myid == 0)
    {
      st = MPI_Wtime();
      cout << " Solving for pressure now";
    }
    pressurePoisson.solve(RHS, &Kimera, 0, allNeumann, noBC);
    if (it % printIntervall == 0 && myid == 0)
    {
      cout << ", solution time was: " << MPI_Wtime() - st << " [s]" << endl;
    }
    // 
    // Then a corrected velocity field
    // is computed using the 
    // Adams-Moulton(2) (=Crank-Nicolson)
    //
    for (k = 0; k < Kimera.nrGrids(); k++)
    {
      Index ui = U.getDiscrBounds(0, k);      Index vi = V.getDiscrBounds(0, k);
      Index uj = U.getDiscrBounds(1, k);      Index vj = V.getDiscrBounds(1, k);
      
      U(k)(ui, uj) = oldU(k)(ui, uj) +
	dt*( 
	    - 0.5*( + 1.*(U(k)(ui, uj)*U.x(k) + V(k)(ui, uj)*U.y(k))
		    + 1.*(oldU(k)(ui, uj)*oldU.x(k) + oldV(k)(ui, uj)*oldU.y(k)) )
	    + 0.5*( ins::viscosityG*oldU.laplacian(k) )
	    + U.forcing(&ins::zero, k, 0.5*dt)
	    - 0.5*( 1.*pressure.x(k, ui, uj) + 1.*oldpressure.x(k, ui, uj) )
	     );
      V(k)(vi, vj) = oldV(k)(vi, vj) + 
	dt*( 
	    - 0.5*( + 1.*(tempU(k)(vi, vj)*V.x(k) + tempV(k)(vi, vj)*V.y(k))
		    + 1.*(oldU(k)(vi, vj)*oldV.x(k) + oldV(k)(vi, vj)*oldV.y(k)) )
	    + 0.5*( ins::viscosityG*oldV.laplacian(k) )
	    + V.forcing(&ins::zero, k, 0.5*dt)
	    - 0.5*( 1.*pressure.y(k, vi, vj) + 1.*oldpressure.y(k, vi, vj) )
	     );
      
    }
    if (it % printIntervall == 0 && myid == 0)
    {
      st = MPI_Wtime();
      cout << " Solving for U now";
    }
    U_momentum.solve(U, &Kimera, 0.5*ins::viscosityG*dt);
    
    if (it % printIntervall == 0 && myid == 0)
    {
      cout << ", solution time was: " << MPI_Wtime() - st << " [s]" << endl;
    }
    if (it % printIntervall == 0 && myid == 0)
    {
      st = MPI_Wtime();
      cout << " Solving for V now";
    }
    V_momentum.solve(V, &Kimera, 0.5*ins::viscosityG*dt);

    if (it % printIntervall == 0 && myid == 0)
    {
      cout << ", solution time was: " << MPI_Wtime() - st << " [s]" << endl;
    }
    // check if simulation can be stopped
    if (U.getTime() > stoptime)
    {
      break;
    }
      
    if (it % printIntervall == 0 && myid ==0)
    {
      cout << "--------DONE--------";
      cout << "|------------------------------------------------------------|" 
	   << endl;
    }
    if (it % saveIntervall == 0)
    {
      fileCounter++;
      fileString.str("");
      fileString << RESULTS << "_" << fileCounter << ".h5";
      
      if (myid == 0)
      {
	cout << "\nSaving result to file " << fileString.str().c_str() << endl;
      }
      U.createHDF5File(fileString.str().c_str());
      Kimera.saveCoordinatesToHDF5File(fileString.str().c_str());
      
      U.saveToHDF5File(fileString.str().c_str(), "u-velocity");
      V.saveToHDF5File(fileString.str().c_str(), "v-velocity");
      pressure.saveToHDF5File(fileString.str().c_str(), "pressure");
    }
  }

  // 
  // Compute the final pressure field
  //
  U.setDivergenceToZeroAtBoundaries(V, xComponentOfVector);
  V.setDivergenceToZeroAtBoundaries(U, yComponentOfVector);
  for (k=0; k<Kimera.nrGrids(); k++)
  {
    Index ix = pressure.getDiscrBounds(0, k);
    Index iy = pressure.getDiscrBounds(1, k);
    
    RHS(k)(ix, iy) = 
      - 2.*U.y(k, ix, iy)*V.x(k, ix, iy)
      - U.x(k, ix, iy)*U.x(k, ix, iy) - V.y(k, ix, iy)*V.y(k, ix, iy);
  }
  ins::setPressureBCs(RHS, pressure, U, V, &Kimera);
  pressurePoisson.solve(RHS, &Kimera, 0, allNeumann, noBC);
  //
  // Timestepping is finished here
  //
  
  double compTime = MPI_Wtime() - startTime;
  double maxTime;
  MPI_Reduce(&compTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, myComm);
  
  if (myid == 0)
  {
    cout << "\nTime-stepping finished\n";
    cout << "  Simulated time is: " << U.getTime() << endl;
    cout << "  Wall clock execution time was " << maxTime << "s\n";
  }
  if (twilightFlow)
  {
    if (myid == 0)
    {
      cout << "Checking U-solution :\n";
    }  
    U.checkErrors(&ins::twilight_U);
      
    if (myid == 0)
    {
      cout << "Checking V-solution :\n";
    }
      
    V.checkErrors(&ins::twilight_V);
      
    if (myid == 0)
    {
      cout << "Checking P-solution :\n";
    }  
    pressure.checkErrors(&ins::twilight_P);
  }
  // Exit PETSc
  PetscFinalize();
  
  // Exit P++ and MPI
  Optimization_Manager::Exit_Virtual_Machine();
  
  return 0;
}
