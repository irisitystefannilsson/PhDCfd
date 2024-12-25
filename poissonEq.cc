#include <mpi.h>
#include "petscsystypes.h"

// OGEquation.hh contains declarations  
// for implicit equations Ax=b
#include "OGEquation.hh"
// GridFunc.hh contains
// declarations for gridFunctions
// and their methods
#include "GridFunc.hh"

// This is necessary because Aztec defines
// min|max macros
#undef min
#undef max
using std::cout;
using std::endl;

doubleArray
twilight_Sol(const doubleArray& x, 
	     const doubleArray& y, 
	     const double t, 
	     int grid, 
	     int side)
{
  doubleArray resArray(x);
  resArray = sin(x)*sin(x)*cos(y)*cos(y) + cos(x)*sin(y);
  return resArray;
}

doubleArray
twilight_X(const doubleArray& x, 
	   const doubleArray& y, 
	   const double t, 
	   int grid, 
	   int side)
{
  doubleArray resArray(x);
  resArray = 2*cos(x)*sin(x)*cos(y)*cos(y) - sin(x)*sin(y);
  return resArray;
}

doubleArray
twilight_Y(const doubleArray& x, 
	   const doubleArray& y, 
	   const double t, 
	   int grid, 
	   int side)
{
  doubleArray resArray(x);
  resArray = -2*sin(x)*sin(x)*cos(y)*sin(y) + cos(x)*cos(y);
  return resArray;
}

doubleArray
twilight_forcing(const doubleArray& x, 
		 const doubleArray& y, 
		 const double t, 
		 int grid, 
		 int side)
{
  doubleArray resArray(x);
  resArray = 2*sin(x)*sin(x)*sin(y)*sin(y)
    - 2*cos(x)*sin(y)
    - 4*sin(x)*sin(x)*cos(y)*cos(y)
    + 2*cos(x)*cos(x)*cos(y)*cos(y);
  return resArray;
}


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
  if (flg != PETSC_TRUE) {
    std::cerr << "Hey!, we need a grid from somewhere (no filename given, exiting...)... \n\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  CompositeGrid chimera(filename);
  // convergence check?
  bool twilightFlow = false;
  PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "--twilight", &flg);
  if (flg == PETSC_TRUE) {
    twilightFlow = true;
  }
  bool allNeumann = false;
  PetscOptionsHasName(PETSC_NULL, PETSC_NULL, "--allNeumann", &flg);
  if (flg == PETSC_TRUE) {
    allNeumann = true;
  }
  // solution
  gridFunction sol(chimera);
  // gridFunction to be used as right-hand-side
  // in the poisson equation
  gridFunction rhs(chimera);

  // reset basic mask in composite grid
  chimera.simplifyFlagValues();

  OGEquation poissonEq(&chimera, PETSc);

  rhs.setBC(scalar); 
  if (allNeumann) {
    sol.setBC(scalar); 
  } else {
    sol.setBC(DirichletAllOver);
  }
  sol.setFlagValues();

  poissonEq.setUnknown(sol);

  poissonEq.setOperator(Poisson, allNeumann);

  if (twilightFlow) {
    sol.set_Dirichlet(&twilight_Sol); sol.set_bc_x(&twilight_X); sol.set_bc_y(&twilight_Y);
  }
  sol.initializeGridData();
  rhs.set_Dirichlet(&twilight_forcing);
  rhs.initializeGridData();

  poissonEq.setOption(ksp_method, Gmres);

  // Solve the Poisson equation
  poissonEq.solve(rhs, &chimera, 0, allNeumann);
  // Check errors
  sol.checkErrors(&twilight_Sol);

  chimera.saveCoordinatesToFile();
  sol.saveToFile("poissSol");
  //poissonEq.writeMatrixOperatorToFile("laplaceOp");
  // Exit PETSc
  PetscFinalize();
  // Exit P++ and MPI
  Optimization_Manager::Exit_Virtual_Machine();
  
  return 0;
}
