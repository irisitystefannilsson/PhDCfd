////////////////////////////////////////
//                                    //
// Miniature code that solves the     //
// heat equation on                   //
// an Overlapping grid generated      //
// by Xcog                            //
//                                    //
////////////////////////////////////////

// OGEquation.hh contains declarations  
// for implicit equations Ax=b
#include <mpi.h>

#include "OGEquation.hh"

// GridFunc.hh contains
// declarations for gridFunctions
// and their methods
#include "GridFunc.hh"
#include "P++.h"
#include "partitioning.h"

// This is necessary because Aztec defines
// min|max macros
#undef min
#undef max

#include "ins_sharp.hh"

// File to save data to
#define RESULTS "heat_data"

#include <algorithm>

#include <sstream>

////////////////////////////////////////
//                                    //
// Main program starts here           //
//                                    //
////////////////////////////////////////

int 
main(int argc, char** argv)
{
  ios::sync_with_stdio();
  
  int Number_of_Processors = 0;
  
  // flg is used to pick out 
  // arguments from argv
  PetscBool flg;

  // Initialize P++ and MPI 
  Optimization_Manager::Initialize_Virtual_Machine("", Number_of_Processors, argc, argv);

  // initialize PETSc
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);

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
  PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "--filename", filename, 100, &flg);
  if (flg != PETSC_TRUE) {
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

  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "--ingrid", &ins::grid, &flg);
  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "--inside", &ins::side, &flg);  

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "--viscosity", &ins::viscosityG, &flg);

  // 
  //  Declare some gridfunctions living on the Kimera grid
  //

  // U and V are the velocity components
  gridFunction U(Kimera, starttime);
  
  // reset basic mask in composite grid
  Kimera.simplifyFlagValues();

  // We want to use a semi-implicit
  // time-stepping algorithm so we
  // create Overlapping grid equations
  // to be solved for the velocity
  // components.
  //   In the declaration we also
  // decide which sparse solver should
  // be used (Aztec or PETSc)
  OGEquation U_momentum(&Kimera, Aztec);

  
  // Set the type of boundary condition to be used
  // for the gridFunctions (to be used together
  // with information from the Xcog-file)
  U.setBC(xComponentOfVector);

  U.setFlagValues(); 

  U_momentum.setUnknown(U); 
  
  // Set type of equation system for the
  // Overlapping grid equations
  // Heat-eqations for U,V
  U_momentum.setOperator(Poisson_with_TD_diagonal); 

  U.set_Dirichlet(& ins::inlet_U); U.set_bc_x(& ins::zero); U.set_bc_y(& ins::zero);
  U.set_Dirichlet_time_derivative(& ins::zero);

  U.initializeGridData();

  // Some variables used in the time-loop
  int maxNits = 1000;
  int printIntervall = 20;
  int saveIntervall = 20;
  int fileCounter = 0;
  double cfl = 0.5;
  stringstream fileString;
  
  double stoptime = 0.1, startTime = MPI_Wtime();
  
  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "--stoptime", &stoptime, &flg);

  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "--maxNits", &maxNits, &flg); 

  PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "--saveIntervall", &saveIntervall, &flg); 

  PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "--cfl", &cfl, &flg);  

  // Set options used for the solution of 
  // the velocity equations
  U_momentum.setOption(precond_method, Jacobi);

  double dt = 0;
  int k;
  
  fileString << RESULTS << "." << fileCounter << ".h5";
	  
  if (myid == 0) {
    std::cout << "\nSaving initial data to file " << fileString.str().c_str() << endl;
  }
  //  U.createHDF5File(fileString.str().c_str());
  //  Kimera.saveCoordinatesToHDF5File(fileString.str().c_str());
  
  U.saveToFile("temperature");

  Kimera.saveCoordinatesToFile();
  //
  // Start time-stepping
  //
  // The equation simulated is 
  // U_t - nu*(U_xx + U_yy) = 0 
  //
  // Most operations are done on the
  // component grid level as I am not
  // confident that the operators for 
  // whole gridFunctions work as intended
  // when gridFunctions with different
  // boundary conditions are mixed in
  // an expression
  //
  if (myid == 0) {
    std::cout << "\nStarting time-stepping now\n";
  }
  for (int it = 0; it < maxNits; it++) {
    // Before the inlet has achieved full
    // speed we use another field to compute
    // the time step
    dt = 0.001*cfl;
    
    if (myid == 0) {
      std::cout << "\nIteration nmbr " 
		<< it << ", time is " 
		<< U.getTime() << endl;
      std::cout << " : dt= " << dt << endl;
    }
    for (k = 0; k < Kimera.nrGrids(); k++) {
      Index ui = U.getDiscrBounds(0, k);
      Index uj = U.getDiscrBounds(1, k);
      
      U(k)(ui, uj) = U(k)(ui, uj) +
	dt*(
	    ins::viscosityG*( U.xx(k) + U.yy(k))
	    );
      
    }
    U_momentum.solve(U, &Kimera, ins::viscosityG*dt);
    U.interpolate();
    
    U.addTime(dt);
    
    U.updateBCs();
    // Advance the velocity one time level
    // and finish the Crank-Nicolson part
    // of the predictor step
    
    // check if simulation can be stopped
    if (U.getTime() > stoptime) {
      break;
    }
    if (it % printIntervall == 0 && myid == 0) {
      std::cout << "|----------------------------------------|" 
		<< endl << endl;
    }
    if (it % saveIntervall == 0) {
      fileCounter++;
      fileString.str("");
      fileString << RESULTS << "." << fileCounter << ".h5";
      
      if (myid == 0) {
	std::cout << "\nSaving result to file " << fileString.str().c_str() << endl;
      }
      U.saveToFile("temperature");
    }
  }
  //
  // Timestepping is finished here
  //
  double compTime = MPI_Wtime() - startTime;
  double maxTime;
  MPI_Reduce(&compTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, myComm);
  
  if (myid == 0) {
    std::cout << "\nTime-stepping finished\n";
    std::cout << "  Simulated time is: " << U.getTime() << endl;
    std::cout << "  Wall clock execution time was " << maxTime << "s\n";
  }
  // Exit PETSc
  PetscFinalize();

  // Exit P++ and MPI
  Optimization_Manager::Exit_Virtual_Machine();
  
  return 0;
}
