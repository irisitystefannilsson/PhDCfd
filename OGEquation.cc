#ifdef _OLD_STL_
#include <math.h>
#else
#include <cmath>
#endif

#include <mpi.h>
#include "petscsys.h"
#include "petscsystypes.h"

#include "OGEquation.hh"


DMSRMatrix::DMSRMatrix()
{
}

DMSRMatrix::~DMSRMatrix()
{
  if (matTypeM == Aztec)
    {
      delete []aM;
      delete []jaM;
    }
  else if (matTypeM == PETSc)
    {
      //      MatDestroy(Pmat);
      delete []localToGlobalM;
      delete []petscLocalRLowM;
      delete []petscLocalSLowM;
      delete []petscLocalRHiM;
      delete []petscLocalSHiM;
      delete []petscStartIndexM;
    }
  delete []rowScalingM;
}


void
DMSRMatrix::setupMatrix(CompositeGrid *Kimera, MatrixOperator operatorType, int N_update, int *update, intArray flagValues[], bool allNeumann, interpolationType iType, gridFunction *visc, bool reset, sparseSolver package, gridFunction *convFieldU, gridFunction *convFieldV)
{
  int me = Communication_Manager::My_Process_Number;
  int nrProcs = Communication_Manager::numberOfProcessors();
  
  MPI_Reduce(&N_update, &globalNrOfUnknownsM, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  NrOfNonZerosM = 0;

  matTypeM = package;

  int iW;
  if (iType == quadratic)
    iW = 3;
  else if (iType == cubic)
    iW = 4;
  else // (iType == linear)
    iW = 2;
  CoeffSendBufferM.redim(2*iW, Kimera->theSendBufferM.getLength(0));
  CoeffReceiveBufferM.redim(2*iW, Kimera->theReceiveBufferM.getLength(0));
  //  cout << "Buffer lengths (1) are =" << Kimera->theSendBuffer.getLength(0) << "," << Kimera->theReceiveBuffer.getLength(0) << endl;
  //  cout << "Buffer lengths (2) are =" << CoeffSendBuffer.getLength(0) << "," << CoeffReceiveBuffer.getLength(0) << endl;
  getInterpolationCoefficients(Kimera, iType);
  
  int discMolecule, intMolecule = iW*iW+1, grid;

  if (operatorType == Poisson || operatorType == Poisson_with_TD_diagonal || operatorType == Grad_Visc_Div || operatorType == Grad_Visc_Div_with_TD_diagonal || operatorType == PUM_with_TD_diagonal || operatorType == Conv_varDiff || operatorType == Conv_varDiff_with_TD_diagonal)
    {
      discMolecule = 9;
    }
  else // ( operatorType == PUM )
    {
      discMolecule = 21;
    }

  int sizeOfMatrix = N_update * discMolecule;

  if (operatorType != PUM)
    {
      for (grid=0; grid<Kimera->nrGrids(); grid++)
	sizeOfMatrix += Kimera->nrOfLocalIntPoints(grid) * (intMolecule - discMolecule);
    }

  int Total_N_update;

  MPI_Allreduce(&N_update, &Total_N_update, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (allNeumann) 
  {
    sizeOfMatrix += N_update;
    if (me == nrProcs-1)
    {
      sizeOfMatrix += Total_N_update;
      update[N_update-1] = Total_N_update-1;
    }
  }
  if (package == PETSc && !reset)
  {
    PetscInt* nnz = new PetscInt[N_update];
    PetscInt* o_nnz = new PetscInt[N_update];
    for (int i = 0; i < N_update; ++i)
    {
      nnz[i] = 11;
      o_nnz[i] = 10;
    }
    if (allNeumann && me == nrProcs-1)
    {
      nnz[N_update - 1] = N_update;
      o_nnz[N_update - 1] = Total_N_update - N_update;
    }
    MatCreateAIJ(PETSC_COMM_WORLD,N_update,N_update,
		 PETSC_DETERMINE,PETSC_DETERMINE,0,
		 nnz,0,o_nnz,&PmatM);
    MatGetOwnershipRange(PmatM, &startRow, &endRow); 
  }
  if (package == PETSc)
  {
    localToGlobalM = new intArray[Kimera->nrGrids()];
    getLocalToGlobal(Kimera);
  }

  if (!reset)
    {
      if (package == Aztec)
	{
	  aM = new double[sizeOfMatrix];
	  jaM = new int[sizeOfMatrix];
	}
      rowScalingM = new double[N_update];
    }

  if (package == Aztec)
    {
      jaM[0] = N_update + 1;
    }
  
  int localrow = 0;
  double r_step, s_step;
  
  for (grid=0; grid<Kimera->nrGrids(); grid++)
  {
    int gType = Kimera->gridTypeM[grid];
    if (gType == 1)
    {
      r_step = Kimera->xM[grid](1,0) - Kimera->xM[grid](0,0);
      s_step = Kimera->yM[grid](0,1) - Kimera->yM[grid](0,0);
    }
    else
    {
      r_step = Kimera->r_stepM[grid];
      s_step = Kimera->s_stepM[grid];
    }
    
    int rdim = Kimera->rdimM[grid];
    int sdim = Kimera->sdimM[grid];

    intSerialArray *Flag = flagValues[grid].getSerialArrayPointer();
    
    doubleSerialArray *X = Kimera->xM[grid].getSerialArrayPointer();
    doubleSerialArray *Y = Kimera->yM[grid].getSerialArrayPointer();

    doubleSerialArray *Xr = Kimera->xrM[grid].getSerialArrayPointer();
    doubleSerialArray *Xs = Kimera->xsM[grid].getSerialArrayPointer();
    doubleSerialArray *Yr = Kimera->yrM[grid].getSerialArrayPointer();
    doubleSerialArray *Ys = Kimera->ysM[grid].getSerialArrayPointer();
    doubleSerialArray *Xrr = Kimera->xrrM[grid].getSerialArrayPointer();
    doubleSerialArray *Xss = Kimera->xssM[grid].getSerialArrayPointer();
    doubleSerialArray *Yrr = Kimera->yrrM[grid].getSerialArrayPointer();
    doubleSerialArray *Yss = Kimera->yssM[grid].getSerialArrayPointer();
    doubleSerialArray *sqrtOfG = Kimera->sqrtOfGM[grid].getSerialArrayPointer();

    if (package == PETSc)
      {
	petscIndexM = localToGlobalM[grid].getSerialArrayPointer();
      }

    doubleSerialArray *viscosity = 0x0;
    if (operatorType == Grad_Visc_Div || operatorType == Grad_Visc_Div_with_TD_diagonal || operatorType == Conv_varDiff || operatorType == Conv_varDiff_with_TD_diagonal)
      {
	viscosity = (*visc)(grid).getSerialArrayPointer(); 
      }
    doubleSerialArray *cFU = 0x0;
    doubleSerialArray *cFV = 0x0;
    if (operatorType == Conv_varDiff || operatorType == Conv_varDiff_with_TD_diagonal)
      {
	cFU = (*convFieldU)(grid).getSerialArrayPointer(); 
	cFV = (*convFieldV)(grid).getSerialArrayPointer(); 
      }

    for (int j=Kimera->lb_1(grid); j<=Kimera->ub_1(grid); j++)
    {
      for (int i=Kimera->lb_0(grid); i<=Kimera->ub_0(grid); i++)
      {
	if ((*Flag)(i,j) == 0) // hole point
	{
	  if (package == Aztec)
	    {
	      aM[localrow] = 1.;
	      jaM[localrow+1] = jaM[localrow];
	      NrOfNonZerosM += 1;
	      
	      if (allNeumann)
		{
		  jaM[localrow+1]++;
		  jaM[jaM[localrow]] = Total_N_update - 1;
		  aM[jaM[localrow]] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  else if (package == PETSc)
	    {
	      vidxM = 0;
	      vM[vidxM] = 1.;
	      nM = 1;
	      //	      idxn[vidxM] = update[localrow];
	      idxnM[vidxM] = (*petscIndexM)(i,j);
	      mM = 1;
	      //	      idxm[vidxM] = update[localrow];
	      idxmM[vidxM] = (*petscIndexM)(i,j);
	      if (allNeumann)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	}
	else if ((*Flag)(i,j) < 0) // interpolation point
	  {
	    if (operatorType == PUM || operatorType == PUM_with_TD_diagonal)
	      {
		double R = 0, overlap = 0;
		double radius = (sqrt(pow((*X)(i,j),2.)+pow((*Y)(i,j),2.))-R)/overlap;
		if (grid == 0)
		  {
		    if (radius <= 0)
		      {
			if (package == Aztec)
			  {
			    aM[localrow] = 1.;
			    jaM[localrow+1] = jaM[localrow];
			    NrOfNonZerosM += 1;
			  }
			else if (package == PETSc)
			  {
			    NrOfNonZerosM += 1;
			    mM = 1;
			    idxmM[0] = (*petscIndexM)(i,j);
			    nM = 1;
			    idxnM[0] = (*petscIndexM)(i,j);
			    vM[0] = 1.;
			  }
		      }
		    else
		      {
			(*Flag)(i,j) = grid + 1;
			insertFD(update[localrow], localrow, rdim, grid, i, j, operatorType, gType, Xr, Xs, Yr, Ys, Xrr, Xss, Yrr, Yss, sqrtOfG, viscosity, r_step, s_step, package, cFU, cFV,X, Y,Kimera);
		      }
		  }
		else if (grid == 1)
		  {
		    if (radius >= 1)
		      {
			if (package == Aztec)
			  {
			    aM[localrow] = 1.;
			    jaM[localrow+1] = jaM[localrow];
			    NrOfNonZerosM += 1;
			  }
			else if (package == PETSc)
			  {
			    NrOfNonZerosM += 1;
			    mM = 1;
			    idxmM[0] = (*petscIndexM)(i,j);
			    nM = 1;
			    idxnM[0] = (*petscIndexM)(i,j);
			    vM[0] = 1.;
			  }
		      }
		    else
		      {
			(*Flag)(i,j) = grid + 1;
			insertFD(update[localrow], localrow, rdim, grid, i, j, operatorType, gType, Xr, Xs, Yr, Ys, Xrr, Xss, Yrr, Yss, sqrtOfG, viscosity, r_step, s_step, package,cFU, cFV,X,Y,Kimera);
		      }
		  }
	      }
	    else
	      {
		insertInterpolationInfo(localrow, update, grid, -(*Flag)(i,j)-1, 
					i, j, Kimera, iType, package);
	      }
	    if (allNeumann)
	      {
		if (package == Aztec)
		  {
		    if (operatorType == PUM || operatorType == PUM_with_TD_diagonal)
		      {
			double R = 0, overlap = 0;
			double radius = (sqrt(pow((*X)(i,j),2.)+pow((*Y)(i,j),2.))-R)/overlap;
			if (grid == 0)
			  {
			    if (radius < 0)
			      {
				jaM[localrow+1]++;
				jaM[jaM[localrow]] = Total_N_update - 1;
				aM[jaM[localrow]] = 1.;
				NrOfNonZerosM += 1;
			      }
			    else 
			      {
				jaM[localrow+1]++;
				if (gType == 1)
				  {
				    jaM[jaM[localrow]+4] = Total_N_update - 1;
				    aM[jaM[localrow]+4] = 1.;
				    NrOfNonZerosM += 1;
				  }
				else
				  {
				    jaM[jaM[localrow]+8] = Total_N_update - 1;
				    aM[jaM[localrow]+8] = 1.;
				    NrOfNonZerosM += 1;
				  }
			      }
			  }
			else if (grid == 1)
			  {
			    if (radius > 1)
			      {
				jaM[localrow+1]++;
				jaM[jaM[localrow]] = Total_N_update - 1;
				aM[jaM[localrow]] = 1.;
				NrOfNonZerosM += 1;
			      }
			    else 
			      {
				jaM[localrow+1]++;
				if (gType == 1)
				  {
				    jaM[jaM[localrow]+4] = Total_N_update - 1;
				    aM[jaM[localrow]+4] = 1.;
				    NrOfNonZerosM += 1;
				  }
				else
				  {
				    jaM[jaM[localrow]+8] = Total_N_update - 1;
				    aM[jaM[localrow]+8] = 1.;
				    NrOfNonZerosM += 1;
				  }
			      }
			  }
		      }
		    else
		      {
			jaM[localrow+1]++;
			jaM[jaM[localrow]+iW*iW] = Total_N_update - 1;
			aM[jaM[localrow]+iW*iW] = 1.;
			NrOfNonZerosM += 1;
		      }
		  }
		else if (package == PETSc)
		  {
		    nM++;
		    idxnM[nM-1] = Total_N_update - 1;
		    vM[nM-1] = 1.;
		    NrOfNonZerosM += 1;
		  }	 
	      }
	  }
	else // some sort of discretization point
	{
	  if ((*Flag)(i,j) == grid + 1) // inner point, apply finite-difference operator
	  {
	    if (operatorType == PUM)
	      {
		insertFD(update[localrow], localrow, rdim, grid, i, j, Poisson, gType, Xr, Xs, Yr, Ys, Xrr, Xss, Yrr, Yss, sqrtOfG, viscosity, r_step, s_step, package,cFU, cFV,X,Y);
	      }
	    else
	      {
		insertFD(update[localrow], localrow, rdim, grid, i, j, operatorType, gType, Xr, Xs, Yr, Ys, Xrr, Xss, Yrr, Yss, sqrtOfG, viscosity, r_step, s_step, package,cFU, cFV);
	      }
	    if (allNeumann)
	    {
	      if (package == Aztec)
		{
		  jaM[localrow+1]++;
		  if (gType == 1)
		    {
		      jaM[jaM[localrow]+4] = Total_N_update - 1;
		      aM[jaM[localrow]+4] = 1.;
		      NrOfNonZerosM += 1;
		    }
		  else
		    {
		      jaM[jaM[localrow]+8] = Total_N_update - 1;
		      aM[jaM[localrow]+8] = 1.;
		      NrOfNonZerosM += 1;
		    }
		}
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }

	  else if ((*Flag)(i,j) == grid + 1000) // Extrapolation, low_r
	  {
	    if (package == Aztec)
	      {
		aM[localrow] = 1.;
		jaM[localrow+1] = jaM[localrow] + 3;
		jaM[jaM[localrow]] = update[localrow] + 1;
		jaM[jaM[localrow]+1] = update[localrow] + 2;
		jaM[jaM[localrow]+2] = update[localrow] + 3;
		aM[jaM[localrow]] = -3.;
		aM[jaM[localrow]+1] = 3.;
		aM[jaM[localrow]+2] = -1.;
		NrOfNonZerosM += 4;
	      }
	    else if (package == PETSc)
	      {
		mM = 1;
		//		idxmM[0] = update[localrow];
		idxmM[0] = (*petscIndexM)(i,j);
		nM = 4;
		//		idxnM[0] = update[localrow];
		//		idxnM[1] = update[localrow] + 1;
		//		idxnM[2] = update[localrow] + 2;
		//		idxnM[3] = update[localrow] + 3;
		idxnM[0] = (*petscIndexM)(i,j);
		idxnM[1] = (*petscIndexM)(i+1,j);
		idxnM[2] = (*petscIndexM)(i+2,j);
		idxnM[3] = (*petscIndexM)(i+3,j);
		vM[0] = 1.;
		vM[1] = -3.;
		vM[2] = 3.;
		vM[3] = -1.;
		
	      }
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+3] = Total_N_update - 1;
		aM[jaM[localrow]+3] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 1 + 1000) // Dirichlet, low_r
	  {
	    if (package == Aztec)
	      {
		aM[localrow] = 1.;
		jaM[localrow+1] = jaM[localrow];
		NrOfNonZerosM += 1;
	      }
	    else if (package == PETSc)
	      {
		NrOfNonZerosM += 1;
		mM = 1;
		//		idxmM[0] = update[localrow];
		idxmM[0] = (*petscIndexM)(i,j);
		nM = 1;
		//		idxnM[0] = update[localrow];
		idxnM[0] = (*petscIndexM)(i,j);
		vM[0] = 1.;
	      }
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]] = Total_N_update - 1;
		aM[jaM[localrow]] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 2 + 1000) // Neumann, low_r
	  {
	    insertNeumann(update[localrow], localrow, 0, rdim, i+1, j, gType, Xr, Xs, Yr, Ys, sqrtOfG, r_step, s_step, package);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		if (gType == 1)
		  {
		    jaM[jaM[localrow]+1] = Total_N_update - 1;
		    aM[jaM[localrow]+1] = 1.;
		    NrOfNonZerosM += 1;
		  }
		else
		  {
		    jaM[jaM[localrow]+3] = Total_N_update - 1;
		    aM[jaM[localrow]+3] = 1.;
		    NrOfNonZerosM += 1;
		  }
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 4 + 1000) // Robin, low_r
	  {
	    insertRobin(update[localrow], localrow, 0, rdim, i+1, j, gType, Xr, Xs, Yr, Ys, sqrtOfG, r_step, s_step, package);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		if (gType == 1)
		  {
		    jaM[jaM[localrow]+2] = Total_N_update - 1;
		    aM[jaM[localrow]+2] = 1.;
		    NrOfNonZerosM += 1;
		  }
		else
		  {
		    jaM[jaM[localrow]+4] = Total_N_update - 1;
		    aM[jaM[localrow]+4] = 1.;
		    NrOfNonZerosM += 1;
		  }
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 3 + 1000) // periodic, low_r
	  {
	    insertPeriod(update[localrow], localrow, 0, rdim, sdim, package, grid,i,j);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+1] = Total_N_update - 1;
		aM[jaM[localrow]+1] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 2000) // Extrapolation, hi_r
	  {
	    if (package == Aztec)
	      {
		aM[localrow] = 1.;
		jaM[localrow+1] = jaM[localrow] + 3;
		jaM[jaM[localrow]] = update[localrow] - 1;
		jaM[jaM[localrow]+1] = update[localrow] - 2;
		jaM[jaM[localrow]+2] = update[localrow] - 3;
		aM[jaM[localrow]] = -3.;
		aM[jaM[localrow]+1] = 3.;
		aM[jaM[localrow]+2] = -1.;
		NrOfNonZerosM += 4;
	      }
	    else if (package == PETSc)
	      {
		NrOfNonZerosM += 4;
		mM = 1;
		//		idxmM[0] = update[localrow];
		idxmM[0] = (*petscIndexM)(i,j);
		nM = 4;
		//		idxnM[0] = update[localrow];
		//		idxnM[1] = update[localrow] - 1; 
		//		idxnM[2] = update[localrow] - 2;
		//		idxnM[3] = update[localrow] - 3;
		idxnM[0] = (*petscIndexM)(i,j);
		idxnM[1] = (*petscIndexM)(i-1,j);
		idxnM[2] = (*petscIndexM)(i-2,j);
		idxnM[3] = (*petscIndexM)(i-3,j);
		vM[0] = 1.;
		vM[1] = -3.;
		vM[2] = 3.;
		vM[3] = -1.;
	      }
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+3] = Total_N_update - 1;
		aM[jaM[localrow]+3] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 1 + 2000) // Dirichlet, hi_r
	  {
	    if (package == Aztec)
	      {
		aM[localrow] = 1.;
		jaM[localrow+1] = jaM[localrow];
		NrOfNonZerosM += 1;
	      }
	    else if (package == PETSc)
	      {
		NrOfNonZerosM += 1;
		mM = 1;
		//		idxmM[0] = update[localrow];
		idxmM[0] = (*petscIndexM)(i,j);
		nM = 1;
		//		idxnM[0] = update[localrow];
		idxnM[0] = (*petscIndexM)(i,j);
		vM[0] = 1.;
	      }
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]] = Total_N_update - 1;
		aM[jaM[localrow]] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 2 + 2000) // Neumann, hi_r
	  {
	    insertNeumann(update[localrow], localrow, 1, rdim, i-1, j, gType, Xr, Xs, Yr, Ys, sqrtOfG, r_step, s_step, package);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		if (gType == 1)
		  {
		    jaM[jaM[localrow]+1] = Total_N_update - 1;
		    aM[jaM[localrow]+1] = 1.;
		    NrOfNonZerosM += 1;
		  }
		else
		  {
		    jaM[jaM[localrow]+3] = Total_N_update - 1;
		    aM[jaM[localrow]+3] = 1.;
		    NrOfNonZerosM += 1;
		  }
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 4 + 2000) // Robin, hi_r
	  {
	    insertRobin(update[localrow], localrow, 1, rdim, i-1, j, gType, Xr, Xs, Yr, Ys, sqrtOfG, r_step, s_step, package);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		if (gType == 1)
		  {
		    jaM[jaM[localrow]+2] = Total_N_update - 1;
		    aM[jaM[localrow]+2] = 1.;
		    NrOfNonZerosM += 1;
		  }
		else
		  {
		    jaM[jaM[localrow]+4] = Total_N_update - 1;
		    aM[jaM[localrow]+4] = 1.;
		    NrOfNonZerosM += 1;
		  }
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 3 + 2000) // periodic, hi_r
	  {
	    insertPeriod(update[localrow], localrow, 1, rdim, sdim, package, grid,i,j);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+1] = Total_N_update - 1;
		aM[jaM[localrow]+1] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 3000) // Extrapolation, low_s
	  {
	    if (package == Aztec)
	      {
		aM[localrow] = 1.;
		jaM[localrow+1] = jaM[localrow] + 3;
		jaM[jaM[localrow]] = update[localrow] + rdim;
		jaM[jaM[localrow]+1] = update[localrow] + 2*rdim;
		jaM[jaM[localrow]+2] = update[localrow] + 3*rdim;
		aM[jaM[localrow]] = -3.;
		aM[jaM[localrow]+1] = 3.;
		aM[jaM[localrow]+2] = -1.;
		NrOfNonZerosM += 4;
	      }
	    else if (package == PETSc)
	      {
		NrOfNonZerosM += 4;
		mM = 1;
		//		idxmM[0] = update[localrow];
		idxmM[0] = (*petscIndexM)(i,j);
		nM = 4;
		//		idxnM[0] = update[localrow];
		//		idxnM[1] = update[localrow] + rdim;
		//		idxnM[2] = update[localrow] + 2*rdim;
		//		idxnM[3] = update[localrow] + 3*rdim;
		idxnM[0] = (*petscIndexM)(i,j);
		idxnM[1] = (*petscIndexM)(i,j+1);
		idxnM[2] = (*petscIndexM)(i,j+2);
		idxnM[3] = (*petscIndexM)(i,j+3);
		vM[0] = 1.;
		vM[1] = -3.;
		vM[2] = 3.;
		vM[3] = -1.;
	      }
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+3] = Total_N_update - 1;
		aM[jaM[localrow]+3] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 1 + 3000) // Dirichlet, low_s
	  {
	    if (package == Aztec)
	      {
		aM[localrow] = 1.;
		jaM[localrow+1] = jaM[localrow];
		NrOfNonZerosM += 1;
	      }
	    else if (package == PETSc)
	      {
		NrOfNonZerosM += 1;
		mM = 1;
		idxmM[0] = (*petscIndexM)(i,j);
		nM = 1;
		idxnM[0] = (*petscIndexM)(i,j);
		vM[0] = 1.;
	      }
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]] = Total_N_update - 1;
		aM[jaM[localrow]] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 2 + 3000) // Neumann, low_s
	  {
	    insertNeumann(update[localrow], localrow,2, rdim, i, j+1, gType, Xr, Xs, Yr, Ys, sqrtOfG, r_step, s_step, package);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		if (gType == 1)
		  {
		    jaM[jaM[localrow]+1] = Total_N_update - 1;
		    aM[jaM[localrow]+1] = 1.;
		    NrOfNonZerosM += 1;
		  }
		else
		  {
		    jaM[jaM[localrow]+3] = Total_N_update - 1;
		    aM[jaM[localrow]+3] = 1.;
		    NrOfNonZerosM += 1;
		  }
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 4 + 3000) // Robin, low_s
	  {
	    insertRobin(update[localrow], localrow,2, rdim, i, j+1, gType, Xr, Xs, Yr, Ys, sqrtOfG, r_step, s_step, package);
	    if (allNeumann)
	    {
	      if (package == Aztec)
		{
		  jaM[localrow+1]++;
		  if (gType == 1)
		    {
		      jaM[jaM[localrow]+2] = Total_N_update - 1;
		      aM[jaM[localrow]+2] = 1.;
		      NrOfNonZerosM += 1;
		    }
		  else
		    {
		      jaM[jaM[localrow]+4] = Total_N_update - 1;
		      aM[jaM[localrow]+4] = 1.;
		      NrOfNonZerosM += 1;
		    }
		}
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 3 + 3000) // periodic, low_s
	  {
	    insertPeriod(update[localrow], localrow,2, rdim, sdim, package, grid,i,j);
	    if (allNeumann)
	    {
	      if (package == Aztec)
		{
		  jaM[localrow+1]++;
		  jaM[jaM[localrow]+1] = Total_N_update - 1;
		  aM[jaM[localrow]+1] = 1.;
		  NrOfNonZerosM += 1;
		}
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 4000) // Extrapolation, hi_s
	  {
	    if (package == Aztec)
	      {
		aM[localrow] = 1.;
		jaM[localrow+1] = jaM[localrow] + 3;
		jaM[jaM[localrow]] = update[localrow] - rdim;
		jaM[jaM[localrow]+1] = update[localrow] - 2*rdim;
		jaM[jaM[localrow]+2] = update[localrow] - 3*rdim;
		aM[jaM[localrow]] = -3.;
		aM[jaM[localrow]+1] = 3.;
		aM[jaM[localrow]+2] = -1.;
		NrOfNonZerosM += 4;
	      }
	    else if (package == PETSc)
	      {
		NrOfNonZerosM += 4;
		mM = 1;
		//		idxmM[0] = update[localrow];
		idxmM[0] = (*petscIndexM)(i,j);
		nM = 4;
		//		idxnM[0] = update[localrow];
		//		idxnM[1] = update[localrow] - rdim;
		//		idxnM[2] = update[localrow] - 2*rdim;
		//		idxnM[3] = update[localrow] - 3*rdim;
		idxnM[0] = (*petscIndexM)(i,j);
		idxnM[1] = (*petscIndexM)(i,j-1);
		idxnM[2] = (*petscIndexM)(i,j-2);
		idxnM[3] = (*petscIndexM)(i,j-3);
		vM[0] = 1.;
		vM[1] = -3.;
		vM[2] = 3.;
		vM[3] = -1.;
	      }
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+3] = Total_N_update - 1;
		aM[jaM[localrow]+3] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 1 + 4000) // Dirichlet, hi_s
	  {
	    if (package == Aztec)
	      {
		aM[localrow] = 1.;
		jaM[localrow+1] = jaM[localrow];
		NrOfNonZerosM += 1;	
	      }
	    else if (package == PETSc)
	      {
		NrOfNonZerosM += 1;
		mM = 1;
		//		idxmM[0] = update[localrow];
		idxmM[0] = (*petscIndexM)(i,j);
		nM = 1;
		//		idxnM[0] = update[localrow];
		idxnM[0] = (*petscIndexM)(i,j);
		vM[0] = 1.;
	      }
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]] = Total_N_update - 1;
		aM[jaM[localrow]] = 1.;
		NrOfNonZerosM += 1;
	      }
	    else if (package == PETSc)
	      {
		nM++;
		idxnM[nM-1] = Total_N_update - 1;
		vM[nM-1] = 1.;
		NrOfNonZerosM += 1;
	      }
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 2 + 4000) // Neumann, hi_s
	  {
	    insertNeumann(update[localrow], localrow,3, rdim, i, j-1, gType, Xr, Xs, Yr, Ys, sqrtOfG, r_step, s_step, package);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		if (gType == 1)
		  {
		    jaM[jaM[localrow]+1] = Total_N_update - 1;
		    aM[jaM[localrow]+1] = 1.;
		    NrOfNonZerosM += 1;
		  }
		else
		  {
		    jaM[jaM[localrow]+3] = Total_N_update - 1;
		    aM[jaM[localrow]+3] = 1.;
		    NrOfNonZerosM += 1;
		  }
	      }
	      else if (package == PETSc)
	      {
		nM++;
		idxnM[nM-1] = Total_N_update - 1;
		vM[nM-1] = 1.;
		NrOfNonZerosM += 1;
	      }
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 4 + 4000) // Robin, hi_s
	  {
	    insertRobin(update[localrow], localrow,3, rdim, i, j-1, gType, Xr, Xs, Yr, Ys, sqrtOfG, r_step, s_step, package);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		if (gType == 1)
		  {
		    jaM[jaM[localrow]+2] = Total_N_update - 1;
		    aM[jaM[localrow]+2] = 1.;
		    NrOfNonZerosM += 1;
		  }
		else
		  {
		    jaM[jaM[localrow]+4] = Total_N_update - 1;
		    aM[jaM[localrow]+4] = 1.;
		    NrOfNonZerosM += 1;
		  }
	      }
	      else if (package == PETSc)
	      {
		nM++;
		idxnM[nM-1] = Total_N_update - 1;
		vM[nM-1] = 1.;
		NrOfNonZerosM += 1;
	      }
	    }
	  }
	  else if ((*Flag)(i,j) == grid + 3 + 4000) // periodic, hi_s
	  {
	    insertPeriod(update[localrow], localrow,3, rdim, sdim, package, grid,i,j);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+1] = Total_N_update - 1;
		aM[jaM[localrow]+1] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == 10000) // corner low_r, low_s
	    {
	    insertCorner(update[localrow], localrow,rdim, 0, package,i,j);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+3] = Total_N_update - 1;
		aM[jaM[localrow]+3] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == 20000) // corner low_r, hi_s
	  {
	    insertCorner(update[localrow], localrow,rdim, 1, package,i,j);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+3] = Total_N_update - 1;
		aM[jaM[localrow]+3] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == 30000) // corner hi_r, low_s
	  {
	    insertCorner(update[localrow], localrow,rdim, 2, package,i,j);
	    if (allNeumann)
	    {
	      if (package == Aztec)
	      {
		jaM[localrow+1]++;
		jaM[jaM[localrow]+3] = Total_N_update - 1;
		aM[jaM[localrow]+3] = 1.;
		NrOfNonZerosM += 1;
	      }
	      else if (package == PETSc)
		{
		  nM++;
		  idxnM[nM-1] = Total_N_update - 1;
		  vM[nM-1] = 1.;
		  NrOfNonZerosM += 1;
		}
	    }
	  }
	  else if ((*Flag)(i,j) == 40000) // corner hi_r, hi_s
	  {
	    insertCorner(update[localrow], localrow,rdim, 3, package,i,j);
	    if (allNeumann)
	      {
		if (package == Aztec)
		  {
		    jaM[localrow+1]++;
		    jaM[jaM[localrow]+3] = Total_N_update - 1;
		    aM[jaM[localrow]+3] = 1.;
		    NrOfNonZerosM += 1;
		  }
		else if (package == PETSc)
		  {
		    nM++;
		    idxnM[nM-1] = Total_N_update - 1;
		    vM[nM-1] = 1.;
		    NrOfNonZerosM += 1;
		  }
	      }
	  }
	}
	//	cout << "i= " << i << " j= " << j << " localrow= " << localrow << " jaM[localrow]= " << jaM[localrow] << " aM[localrow]= " << aM[localrow] << " (*Flag)(i,j) = " << (*Flag)(i,j) << endl;
	
	if (package == Aztec)
	  {
	    //----------------------------------------
	    // Scale each row in Matrix so that
	    // row sum is 1
	    //----------------------------------------
	    int col;
	    double sum;
	    sum = fabs(aM[localrow]);
	    for (col=jaM[localrow]; col<jaM[localrow+1]; col++)
	      {
		sum += fabs(aM[col]);
		//		if (jaM[col] == 842)
		//		  {
		//		    cout << "YELL!\n";
		//		    cout << (*Flag)(i,j) << " grid = " << grid << endl;
		//		    cout << " i,j = " << i << "," << j << endl;
		//		  }
	      }
	    //	cout << "Sum = " << sum << endl;
	    sum = sum*fabs(aM[localrow])/aM[localrow];
	    aM[localrow] /= sum;
	    for (col=jaM[localrow]; col<jaM[localrow+1]; col++)
	      aM[col] /= sum;
	    
	    rowScalingM[localrow] = 1./ sum;
	    
	    localrow++;
	    
	  }
	else if (package == PETSc)
	  {
	    int col;
	    double sum;
	    sum = 0.;
	    for (col=0; col<nM; col++)
	      {
		sum += fabs((double) vM[col]);
	      }
	    sum = sum*fabs(vM[0])/vM[0];
	    for (col=0; col<nM; col++)
	      {
		vM[col] /= sum;
		//		cout << idxnM[col] << " ";
	      }
	    //	    cout << (*Flag)(i,j) << endl;
	    rowScalingM[localrow] = 1./ sum;

	    localrow++;

	    //	    cout << "I am No. " << Communication_Manager::My_Process_Number << " , I now insert line no. " << idxmM[0] << endl;
	    //	    cout << " ----Column numbers are ";
	    //	    for (int c=0; c<n; c++) cout << idxnM[c] << " ";
	    //	    cout << endl;
	    MatSetValues(PmatM, mM, idxmM, nM, idxnM, vM, INSERT_VALUES);
	  }
      }
    }
  }
  if (allNeumann)
    {
      if (me == nrProcs-1)
	{
	  if (package == Aztec)
	    {
	      NrOfNonZerosM += Total_N_update;
	      
	      jaM[localrow+1] = jaM[localrow] + Total_N_update - 1;
	      aM[localrow] = 0.;
	      
	      for (int col=0; col<Total_N_update-1; col++)
		{
		  jaM[jaM[localrow]+col] = col;
		  aM[jaM[localrow]+col] = 1./(double) (Total_N_update-1);
		}
	    }
	  else if (package == PETSc)
	    {
	      int *lastRow = new int[Total_N_update];
	      PetscScalar *lastRowVals = new PetscScalar[Total_N_update];
	      for (int col=0; col<Total_N_update-1; col++)
		{
		  lastRow[col] = col;
		  lastRowVals[col] = 1/(double) (Total_N_update-1);
		}
	      lastRow[Total_N_update-1] = Total_N_update-1;
	      lastRowVals[Total_N_update-1] = 0.;
	      mM = 1;
	      idxmM[0] = Total_N_update - 1;
	      nM = Total_N_update;
	      MatSetValues(PmatM, mM, idxmM, nM, lastRow, lastRowVals, INSERT_VALUES);
	
	      delete []lastRowVals;
	      delete []lastRow;
	    }
	}
    }
  if (package == PETSc)
    {
      MatAssemblyBegin(PmatM, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(PmatM, MAT_FINAL_ASSEMBLY);
    }
}

void
DMSRMatrix::insertFD(int row, int localrow, int rdim, int grid, int i, int j, MatrixOperator operatorType, int gridType, doubleSerialArray *Xr, doubleSerialArray *Xs, doubleSerialArray *Yr, doubleSerialArray *Ys, doubleSerialArray *Xrr, doubleSerialArray *Xss, doubleSerialArray *Yrr, doubleSerialArray *Yss, doubleSerialArray *sqrtOfG, doubleSerialArray *viscosity, double r_step, double s_step, sparseSolver package, doubleSerialArray *convFieldU, doubleSerialArray *convFieldV, doubleSerialArray *X, doubleSerialArray *Y, CompositeGrid *Kimera)
{
  if (operatorType == PUM)
    {
      if (gridType == 1) // cartesian grid, also means the other grid is curvi-linear
	{
	  NrOfNonZerosM += 5; // laplacian part ...
	  NrOfNonZerosM += 12; // ... and discretization part from other grid

	  
	  if (package == Aztec)
	    {
	      jaM[localrow+1] = jaM[localrow] + 16;
	      jaM[jaM[localrow]] = row - rdim;
	      jaM[jaM[localrow]+1] = row - 1;
	      jaM[jaM[localrow]+2] = row + 1;
	      jaM[jaM[localrow]+3] = row + rdim;
	      
	      aM[localrow] = - 2./pow(s_step,2.) - 2./pow(r_step,2.);
	      aM[jaM[localrow]] = 1./pow(s_step,2.);
	      aM[jaM[localrow]+1] = 1./pow(r_step,2.);
	      aM[jaM[localrow]+2] = 1./pow(r_step,2.);
	      aM[jaM[localrow]+3] = 1./pow(s_step,2.);
#ifdef PUMP
	      double fi_xx = phi_xx((*X)(i,j),(*Y)(i,j));
	      double fi_yy = phi_yy((*X)(i,j),(*Y)(i,j));
	      double fi_x = phi_x((*X)(i,j),(*Y)(i,j));
	      double fi_y = phi_y((*X)(i,j),(*Y)(i,j));

	      aM[localrow] += fi_xx + fi_yy;
	      aM[jaM[localrow]] += fi_y*(-0.5/s_step);
	      aM[jaM[localrow]+1] += fi_x*(-0.5/r_step);
	      aM[jaM[localrow]+2] += fi_x*(0.5/r_step);
	      aM[jaM[localrow]+3] += fi_y*(0.5/s_step);

	      double r_step1 = Kimera->r_stepM[1];
	      double s_step1 = Kimera->s_stepM[1];
	      doubleSerialArray *Xr1 = Kimera->xrM[1].getSerialArrayPointer();
	      doubleSerialArray *Xs1 = Kimera->xsM[1].getSerialArrayPointer();
	      doubleSerialArray *Yr1 = Kimera->yrM[1].getSerialArrayPointer();
	      doubleSerialArray *Ys1 = Kimera->ysM[1].getSerialArrayPointer();
	      doubleSerialArray *sqrtOfG1 = Kimera->sqrtOfGM[1].getSerialArrayPointer();

	      int index = Kimera->getInterpIndex(i,j,grid);
	      
	      int i_loc = Kimera->My_i_L[grid](index);
	      int j_loc = Kimera->My_j_L[grid](index);
	      
	      int start = 0;
	      int offset = 0;
	      
	      Index all;
	      
	      for (int k=0; k<1; k++)
		{
		  start += Kimera->rdimM[k] * Kimera->sdimM[k];
		}
	      for (int g=0; g<grid; g++)
		offset += Kimera->myIntPoints[g].getLength(0);
	      

	      jaM[jaM[localrow]+4] = start + j_loc*(Kimera->rdimM[1]) + i_loc;
	      jaM[jaM[localrow]+5] = start + j_loc*(Kimera->rdimM[1]) + i_loc + 1;
	      jaM[jaM[localrow]+6] = start + j_loc*(Kimera->rdimM[1]) + i_loc + Kimera->rdimM[1];
	      jaM[jaM[localrow]+7] = start + j_loc*(Kimera->rdimM[1]) + i_loc + Kimera->rdimM[1] + 1;
	      jaM[jaM[localrow]+8] = start + j_loc*(Kimera->rdimM[1]) + i_loc - Kimera->rdimM[1];
	      jaM[jaM[localrow]+9] = start + j_loc*(Kimera->rdimM[1]) + i_loc - Kimera->rdimM[1] + 1;
	      jaM[jaM[localrow]+10] = start + j_loc*(Kimera->rdimM[1]) + i_loc -1;
	      jaM[jaM[localrow]+11] = start + j_loc*(Kimera->rdimM[1]) + i_loc + 2;
	      jaM[jaM[localrow]+12] = start + j_loc*(Kimera->rdimM[1]) + i_loc + Kimera->rdimM[1] - 1;
	      jaM[jaM[localrow]+13] = start + j_loc*(Kimera->rdimM[1]) + i_loc + Kimera->rdimM[1] + 2;
	      jaM[jaM[localrow]+14] = start + j_loc*(Kimera->rdimM[1]) + i_loc + 2*Kimera->rdimM[1];
	      jaM[jaM[localrow]+15] = start + j_loc*(Kimera->rdimM[1]) + i_loc + 2*Kimera->rdimM[1] + 1;

	      aM[jaM[localrow]+4] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(fi_xx+fi_yy);
	      aM[jaM[localrow]+5] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(fi_xx+fi_yy);
	      aM[jaM[localrow]+6] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*(fi_xx+fi_yy);
	      aM[jaM[localrow]+7] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*(fi_xx+fi_yy);
	      double coeff01[3][3] = {0,1.,0, 0,0,0, 0,0,0};
	      aM[jaM[localrow]+8] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01) + 2*fi_y*yCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step,s_step,coeff01));
	      double coeff10[3][3] = {0,0,0, 1.,0,0, 0,0,0};
	      aM[jaM[localrow]+10] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10) + 2*fi_y*yCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10));
	      double coeff12[3][3] = {0,0,0, 0,0,1, 0,0,0};
	      aM[jaM[localrow]+5] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12) + 2*fi_y*yCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12));
	      double coeff21[3][3] = {0,0,0, 0,0,0, 0,1,0};
	      aM[jaM[localrow]+6] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21) + 2*fi_y*yCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21));

	      aM[jaM[localrow]+9] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01) + 2*fi_y*yCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01));
	      aM[jaM[localrow]+4] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10) + 2*fi_y*yCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10));
	      aM[jaM[localrow]+11] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12) + 2*fi_y*yCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12));
	      aM[jaM[localrow]+7] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21) + 2*fi_y*yCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21));

	      aM[jaM[localrow]+4] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01) + 2*fi_y*yCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01));
	      aM[jaM[localrow]+12] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10) + 2*fi_y*yCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10));
	      aM[jaM[localrow]+7] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12) + 2*fi_y*yCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12));
	      aM[jaM[localrow]+14] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21) + 2*fi_y*yCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21));
	      
	      aM[jaM[localrow]+5] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01) + 2*fi_y*yCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01));
	      aM[jaM[localrow]+6] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10) + 2*fi_y*yCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10));
	      aM[jaM[localrow]+13] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12) + 2*fi_y*yCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12));
	      aM[jaM[localrow]+15] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21) + 2*fi_y*yCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21));
#endif
	    }
	  else if (package == PETSc)
	    {
	      mM = 1;
	      //	      idxmM[0] = row;
	      idxmM[0] = (*petscIndexM)(i,j);
	      nM = 17;
	      //	      idxnM[0] = row;
	      //	      idxnM[1] = row - rdim;
	      //	      idxnM[2] = row - 1;
	      //	      idxnM[3] = row + 1;
	      //	      idxnM[4] = row + rdim;
	      idxnM[0] = (*petscIndexM)(i,j);
	      idxnM[1] = (*petscIndexM)(i,j-1);
	      idxnM[2] = (*petscIndexM)(i-1,j);
	      idxnM[3] = (*petscIndexM)(i+1,j);
	      idxnM[4] = (*petscIndexM)(i,j+1);

	      vM[0] = - 2./pow(s_step,2.) - 2./pow(r_step,2.);
	      vM[1] = 1./pow(s_step,2.);
	      vM[2] = 1./pow(r_step,2.);
	      vM[3] = 1./pow(r_step,2.);
	      vM[4] = 1./pow(s_step,2.);

#ifdef PUMP
	      double fi_xx = phi_xx((*X)(i,j),(*Y)(i,j));
	      double fi_yy = phi_yy((*X)(i,j),(*Y)(i,j));
	      double fi_x = phi_x((*X)(i,j),(*Y)(i,j));
	      double fi_y = phi_y((*X)(i,j),(*Y)(i,j));

	      vM[0] += fi_xx + fi_yy;
	      vM[1] += fi_y*(-0.5/s_step);
	      vM[2] += fi_x*(-0.5/r_step);
	      vM[3] += fi_x*(0.5/r_step);
	      vM[4] += fi_y*(0.5/s_step);

	      int index = Kimera->getInterpIndex(i,j,grid);
	      int i_loc = Kimera->My_i_L[grid](index);
	      int j_loc = Kimera->My_j_L[grid](index);
	      

	      double r_step1 = Kimera->r_stepM[1];
	      double s_step1 = Kimera->s_stepM[1];
	      doubleSerialArray *Xr1 = Kimera->xrM[1].getSerialArrayPointer();
	      doubleSerialArray *Xs1 = Kimera->xsM[1].getSerialArrayPointer();
	      doubleSerialArray *Yr1 = Kimera->yrM[1].getSerialArrayPointer();
	      doubleSerialArray *Ys1 = Kimera->ysM[1].getSerialArrayPointer();
	      doubleSerialArray *sqrtOfG1 = Kimera->sqrtOfGM[1].getSerialArrayPointer();
	      int start = 0;
	      int offset = 0;
	      
	      Index all;
	      
	      for (int g=0; g<grid; g++)
		offset += Kimera->myIntPoints[g].getLength(0);
	      

	      int lRDim;
	      
	      lRDim =  petscLocalRHi[1](0) - petscLocalRLow[1](0) + 1;

	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc-petscLocalSLow[1](0))
		+ (i_loc-petscLocalRLow[1](0));
	      
	      idxnM[5] = start;
	      
	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc-petscLocalSLow[1](0))
		+ (i_loc+1-petscLocalRLow[1](0));
	      
	      idxnM[6] = start;
	      
	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc+1-petscLocalSLow[1](0))
		+ (i_loc-petscLocalRLow[1](0));
	      
	      idxnM[7] = start;
	      
	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc+1-petscLocalSLow[1](0))
		+ (i_loc+1-petscLocalRLow[1](0));
	      
	      idxnM[8] = start;
	      
	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc-1-petscLocalSLow[1](0))
		+ (i_loc-petscLocalRLow[1](0));
	      
	      idxnM[9] = start;
	      
	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc-1-petscLocalSLow[1](0))
		+ (i_loc+1-petscLocalRLow[1](0));
	      
	      idxnM[10] = start;
	      
	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc-petscLocalSLow[1](0))
		+ (i_loc-1-petscLocalRLow[1](0));
	      
	      idxnM[11] = start;	      

	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc-petscLocalSLow[1](0))
		+ (i_loc+2-petscLocalRLow[1](0));
	      
	      idxnM[12] = start;

	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc+1-petscLocalSLow[1](0))
		+ (i_loc-1-petscLocalRLow[1](0));
	      
	      idxnM[13] = start;
	      
	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc+1-petscLocalSLow[1](0))
		+ (i_loc+2-petscLocalRLow[1](0));
	      
	      idxnM[14] = start;

	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc+2-petscLocalSLow[1](0))
		+ (i_loc-petscLocalRLow[1](0));
	      
	      idxnM[15] = start;

	      start = petscStartIndex[1](0)
		+ lRDim*(j_loc+2-petscLocalSLow[1](0))
		+ (i_loc+1-petscLocalRLow[1](0));
	      
	      idxnM[16] = start;


	      vM[5] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(fi_xx+fi_yy);
	      vM[6] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(fi_xx+fi_yy);
	      vM[7] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*(fi_xx+fi_yy);
	      vM[8] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*(fi_xx+fi_yy);
	      double coeff01[3][3] = {0,1.,0, 0,0,0, 0,0,0};
	      vM[9] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01) + 2*fi_y*yCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step,s_step,coeff01));
	      double coeff10[3][3] = {0,0,0, 1.,0,0, 0,0,0};
	      vM[11] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10) + 2*fi_y*yCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10));
	      double coeff12[3][3] = {0,0,0, 0,0,1, 0,0,0};
	      vM[6] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12) + 2*fi_y*yCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12));
	      double coeff21[3][3] = {0,0,0, 0,0,0, 0,1,0};
	      vM[7] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21) + 2*fi_y*yCoeff(i_loc,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21));

	      vM[10] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01) + 2*fi_y*yCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01));
	      vM[5] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10) + 2*fi_y*yCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10));
	      vM[12] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12) + 2*fi_y*yCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12));
	      vM[8] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21) + 2*fi_y*yCoeff(i_loc+1,j_loc,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21));

	      vM[5] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01) + 2*fi_y*yCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01));
	      vM[13] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10) + 2*fi_y*yCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10));
	      vM[8] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12) + 2*fi_y*yCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12));
	      vM[15] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*(2*fi_x*xCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21) + 2*fi_y*yCoeff(i_loc,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21));
	      
	      vM[6] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01) + 2*fi_y*yCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff01));
	      vM[7] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10) + 2*fi_y*yCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff10));
	      vM[14] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12) + 2*fi_y*yCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff12));
	      vM[16] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*(2*fi_x*xCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21) + 2*fi_y*yCoeff(i_loc+1,j_loc+1,Xr1,Xs1,Yr1,Ys1,sqrtOfG1,r_step1,s_step1,coeff21));
#endif
	    }
	}
      else
	{
	  NrOfNonZerosM += 9;
	  NrOfNonZerosM += 12;
	  if (package == Aztec)
	    {
	      jaM[localrow+1] = jaM[localrow] + 20;
	      jaM[jaM[localrow]] = row - rdim - 1;
	      jaM[jaM[localrow]+1] = row - rdim;
	      jaM[jaM[localrow]+2] = row - rdim + 1;
	      jaM[jaM[localrow]+3] = row - 1;
	      jaM[jaM[localrow]+4] = row + 1;
	      jaM[jaM[localrow]+5] = row + rdim - 1;
	      jaM[jaM[localrow]+6] = row + rdim;
	      jaM[jaM[localrow]+7] = row + rdim + 1;
	      //      cout << " column = " <<   row + rdim + 1 << " i,j = " << i << " , " << j << " ";
	      
	      double coeff11[3][3] = {{0,0,0}, {0,1.,0}, {0,0,0}};
	      aM[localrow] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11);
	      double coeff00[3][3] = {{1.,0,0}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00);
	      double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]+1] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01);
	      double coeff02[3][3] = {{0,0,1.}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]+2] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02);
	      double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	      aM[jaM[localrow]+3] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10);
	      double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	      aM[jaM[localrow]+4]= xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12);
	      double coeff20[3][3] = {{0,0,0}, {0,0,0}, {1.,0,0}};
	      aM[jaM[localrow]+5] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20);
	      double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	      aM[jaM[localrow]+6] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21);
	      double coeff22[3][3] = {{0,0,0}, {0,0,0}, {0,0,1.}};
	      aM[jaM[localrow]+7] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22);

#ifdef PUMP
	      double fi_xx = phi_xx((*X)(i,j),(*Y)(i,j));
	      double fi_yy = phi_yy((*X)(i,j),(*Y)(i,j));
	      double fi_x = phi_x((*X)(i,j),(*Y)(i,j));
	      double fi_y = phi_y((*X)(i,j),(*Y)(i,j));
	      
	      aM[localrow] -= fi_xx + fi_yy;

	      aM[localrow] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff11)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff11);

	      aM[jaM[localrow]] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff00)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff00);

	      aM[jaM[localrow]+1] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);

	      aM[jaM[localrow]+2] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff02)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff02);

	      aM[jaM[localrow]+3] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);

	      aM[jaM[localrow]+4] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);

	      aM[jaM[localrow]+5] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff20)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff20);

	      aM[jaM[localrow]+6] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);

	      aM[jaM[localrow]+7] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff22)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff22);

	      double r_step0 = Kimera->xM[0](1,0) - Kimera->xM[0](0,0);
	      double s_step0 = Kimera->yM[0](0,1) - Kimera->yM[0](0,0);

	      int index = Kimera->getInterpIndex(i,j,grid);
	      
	      int i_loc = Kimera->My_i_L[grid](index);
	      int j_loc = Kimera->My_j_L[grid](index);
	      
	      int start = 0;
	      int offset = 0;
	      
	      Index all;
	      
	      for (int g=0; g<grid; g++)
		offset += Kimera->myIntPoints[g].getLength(0);
	      

	      jaM[jaM[localrow]+8] = start + j_loc*(Kimera->rdimM[0]) + i_loc;
	      jaM[jaM[localrow]+9] = start + j_loc*(Kimera->rdimM[0]) + i_loc + 1;
	      jaM[jaM[localrow]+10] = start + j_loc*(Kimera->rdimM[0]) + i_loc + Kimera->rdimM[0];
	      jaM[jaM[localrow]+11] = start + j_loc*(Kimera->rdimM[0]) + i_loc + Kimera->rdimM[0] + 1;
	      jaM[jaM[localrow]+12] = start + j_loc*(Kimera->rdimM[0]) + i_loc - Kimera->rdimM[0];
	      jaM[jaM[localrow]+13] = start + j_loc*(Kimera->rdimM[0]) + i_loc - Kimera->rdimM[0] + 1;
	      jaM[jaM[localrow]+14] = start + j_loc*(Kimera->rdimM[0]) + i_loc -1;
	      jaM[jaM[localrow]+15] = start + j_loc*(Kimera->rdimM[0]) + i_loc + 2;
	      jaM[jaM[localrow]+16] = start + j_loc*(Kimera->rdimM[0]) + i_loc + Kimera->rdimM[0] - 1;
	      jaM[jaM[localrow]+17] = start + j_loc*(Kimera->rdimM[0]) + i_loc + Kimera->rdimM[0] + 2;
	      jaM[jaM[localrow]+18] = start + j_loc*(Kimera->rdimM[0]) + i_loc + 2*Kimera->rdimM[0];
	      jaM[jaM[localrow]+19] = start + j_loc*(Kimera->rdimM[0]) + i_loc + 2*Kimera->rdimM[0] + 1;

	      
	      aM[jaM[localrow]+8] = -CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(fi_xx+fi_yy);
	      aM[jaM[localrow]+9] = -CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(fi_xx+fi_yy);
	      aM[jaM[localrow]+10] = -CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*(fi_xx+fi_yy);
	      aM[jaM[localrow]+11] = -CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*(fi_xx+fi_yy);

	      aM[jaM[localrow]+12] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*2*fi_y*(0.5/s_step0);
	      aM[jaM[localrow]+10] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*-2*fi_y*(0.5/s_step0);

	      aM[jaM[localrow]+14] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*2*fi_x*(0.5/r_step0);
	      aM[jaM[localrow]+9] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*-2*fi_x*(0.5/r_step0);

	      aM[jaM[localrow]+13] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*2*fi_y*(0.5/s_step0);
	      aM[jaM[localrow]+11] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*-2*fi_y*(0.5/s_step0);

	      aM[jaM[localrow]+8] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*2*fi_x*(0.5/r_step0);
	      aM[jaM[localrow]+15] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*-2*fi_x*(0.5/r_step0);

	      aM[jaM[localrow]+8] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*2*fi_y*(0.5/s_step0);
	      aM[jaM[localrow]+18] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*-2*fi_y*(0.5/s_step0);

	      aM[jaM[localrow]+16] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*2*fi_x*(0.5/r_step0);
	      aM[jaM[localrow]+11] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*-2*fi_x*(0.5/r_step0);

	      aM[jaM[localrow]+9] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*2*fi_y*(0.5/s_step0);
	      aM[jaM[localrow]+19] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*-2*fi_y*(0.5/s_step0);

	      aM[jaM[localrow]+10] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*2*fi_x*(0.5/r_step0);
	      aM[jaM[localrow]+17] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*-2*fi_x*(0.5/r_step0);
#endif
	    }
	  else if (package == PETSc)
	    {
	      mM = 1;

	      idxmM[0] = (*petscIndexM)(i,j);
	      nM = 21;

	      idxnM[0] = (*petscIndexM)(i,j);
	      idxnM[1] = (*petscIndexM)(i-1,j-1);
	      idxnM[2] = (*petscIndexM)(i,j-1);
	      idxnM[3] = (*petscIndexM)(i+1,j-1);
	      idxnM[4] = (*petscIndexM)(i-1,j);
	      idxnM[5] = (*petscIndexM)(i+1,j);
	      idxnM[6] = (*petscIndexM)(i-1,j+1);
	      idxnM[7] = (*petscIndexM)(i,j+1);
	      idxnM[8] = (*petscIndexM)(i+1,j+1);

	      double coeff11[3][3] = {{0,0,0}, {0,1.,0}, {0,0,0}};
	      vM[0] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11);
	      double coeff00[3][3] = {{1.,0,0}, {0,0,0}, {0,0,0}};
	      vM[1] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00);
	      double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	      vM[2] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01);
	      double coeff02[3][3] = {{0,0,1.}, {0,0,0}, {0,0,0}};
	      vM[3] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02);
	      double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	      vM[4] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10);
	      double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	      vM[5] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12);
	      double coeff20[3][3] = {{0,0,0}, {0,0,0}, {1.,0,0}};
	      vM[6] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20);
	      double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	      vM[7] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21);
	      double coeff22[3][3] = {{0,0,0}, {0,0,0}, {0,0,1.}};
	      vM[8] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22);
	      
#ifdef PUMP
	      double fi_xx = phi_xx((*X)(i,j),(*Y)(i,j));
	      double fi_yy = phi_yy((*X)(i,j),(*Y)(i,j));
	      double fi_x = phi_x((*X)(i,j),(*Y)(i,j));
	      double fi_y = phi_y((*X)(i,j),(*Y)(i,j));
	      
	      vM[0] -= fi_xx + fi_yy;

	      vM[0] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff11)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff11);

	      vM[1] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff00)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff00);

	      vM[2] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);

	      vM[3] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff02)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff02);

	      vM[4] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);

	      vM[5] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);

	      vM[6] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff20)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff20);

	      vM[7] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);

	      vM[8] -= 2*fi_x*xCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff22)
		+ 2*fi_y*yCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff22);

	      double r_step0 = Kimera->xM[0](1,0) - Kimera->xM[0](0,0);
	      double s_step0 = Kimera->yM[0](0,1) - Kimera->yM[0](0,0);
	      
	      int index = Kimera->getInterpIndex(i,j,grid);
	      int i_loc = Kimera->My_i_L[grid](index);
	      int j_loc = Kimera->My_j_L[grid](index);
	      
	      int start = 0;
	      int offset = 0;
	      
	      Index all;
	      
	      for (int g=0; g<grid; g++)
		offset += Kimera->myIntPoints[g].getLength(0);
	      

	      int iProc, lRDim;
	      
	      lRDim =  petscLocalRHi[0](0) - petscLocalRLow[0](0) + 1;

	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc-petscLocalSLow[0](0))
		+ (i_loc-petscLocalRLow[0](0));
	      
	      idxnM[9] = start;
	      
	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc-petscLocalSLow[0](0))
		+ (i_loc+1-petscLocalRLow[0](0));
	      
	      idxnM[10] = start;
	      
	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc+1-petscLocalSLow[0](0))
		+ (i_loc-petscLocalRLow[0](0));
	      
	      idxnM[11] = start;
	      
	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc+1-petscLocalSLow[0](0))
		+ (i_loc+1-petscLocalRLow[0](0));
	      
	      idxnM[12] = start;
	      
	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc-1-petscLocalSLow[0](0))
		+ (i_loc-petscLocalRLow[0](0));
	      
	      idxnM[13] = start;
	      
	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc-1-petscLocalSLow[0](0))
		+ (i_loc+1-petscLocalRLow[0](0));
	      
	      idxnM[14] = start;
	      
	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc-petscLocalSLow[0](0))
		+ (i_loc-1-petscLocalRLow[0](0));
	      
	      idxnM[15] = start;	      

	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc-petscLocalSLow[0](0))
		+ (i_loc+2-petscLocalRLow[0](0));
	      
	      idxnM[16] = start;

	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc+1-petscLocalSLow[0](0))
		+ (i_loc-1-petscLocalRLow[0](0));
	      
	      idxnM[17] = start;
	      
	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc+1-petscLocalSLow[0](0))
		+ (i_loc+2-petscLocalRLow[0](0));
	      
	      idxnM[18] = start;

	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc+2-petscLocalSLow[0](0))
		+ (i_loc-petscLocalRLow[0](0));
	      
	      idxnM[19] = start;

	      start = petscStartIndex[0](0)
		+ lRDim*(j_loc+2-petscLocalSLow[0](0))
		+ (i_loc+1-petscLocalRLow[0](0));
	      
	      idxnM[20] = start;

	      vM[9] = -CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*(fi_xx+fi_yy);
	      vM[10] = -CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*(fi_xx+fi_yy);
	      vM[11] = -CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*(fi_xx+fi_yy);
	      vM[12] = -CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*(fi_xx+fi_yy);

	      vM[13] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*2*fi_y*(0.5/s_step0);
	      vM[11] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*-2*fi_y*(0.5/s_step0);

	      vM[15] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*2*fi_x*(0.5/r_step0);
	      vM[10] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(2,offset+index)*-2*fi_x*(0.5/r_step0);

	      vM[14] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*2*fi_y*(0.5/s_step0);
	      vM[12] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*-2*fi_y*(0.5/s_step0);

	      vM[9] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*2*fi_x*(0.5/r_step0);
	      vM[16] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(2,offset+index)*-2*fi_x*(0.5/r_step0);

	      vM[9] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*2*fi_y*(0.5/s_step0);
	      vM[19] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*-2*fi_y*(0.5/s_step0);

	      vM[17] = CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*2*fi_x*(0.5/r_step0);
	      vM[12] += CoeffReceiveBuffer(0,offset+index) * CoeffReceiveBuffer(3,offset+index)*-2*fi_x*(0.5/r_step0);

	      vM[10] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*2*fi_y*(0.5/s_step0);
	      vM[20] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*-2*fi_y*(0.5/s_step0);

	      vM[11] += CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*2*fi_x*(0.5/r_step0);
	      vM[18] = CoeffReceiveBuffer(1,offset+index) * CoeffReceiveBuffer(3,offset+index)*-2*fi_x*(0.5/r_step0);
#endif
	    } 
	}
    }
  if (operatorType == Poisson || operatorType == Poisson_with_TD_diagonal || operatorType == PUM_with_TD_diagonal)
    {
      if (gridType == 1) // cartesian grid
	{
	  NrOfNonZerosM += 5;
	  
	  if (package == Aztec)
	    {
	      jaM[localrow+1] = jaM[localrow] + 4;
	      jaM[jaM[localrow]] = row - rdim;
	      jaM[jaM[localrow]+1] = row - 1;
	      jaM[jaM[localrow]+2] = row + 1;
	      jaM[jaM[localrow]+3] = row + rdim;
	      
	      aM[localrow] = - 2./pow(s_step,2.) - 2./pow(r_step,2.);
	      aM[jaM[localrow]] = 1./pow(s_step,2.);
	      aM[jaM[localrow]+1] = 1./pow(r_step,2.);
	      aM[jaM[localrow]+2] = 1./pow(r_step,2.);
	      aM[jaM[localrow]+3] = 1./pow(s_step,2.);
	    }
	  else if (package == PETSc)
	    {
	      mM = 1;
	      //	      idxmM[0] = row;
	      idxmM[0] = (*petscIndexM)(i,j);
	      nM = 5;
	      //	      idxnM[0] = row;
	      //	      idxnM[1] = row - rdim;
	      //	      idxnM[2] = row - 1;
	      //	      idxnM[3] = row + 1;
	      //	      idxnM[4] = row + rdim;
	      idxnM[0] = (*petscIndexM)(i,j);
	      idxnM[1] = (*petscIndexM)(i,j-1);
	      idxnM[2] = (*petscIndexM)(i-1,j);
	      idxnM[3] = (*petscIndexM)(i+1,j);
	      idxnM[4] = (*petscIndexM)(i,j+1);

	      vM[0] = - 2./pow(s_step,2.) - 2./pow(r_step,2.);
	      vM[1] = 1./pow(s_step,2.);
	      vM[2] = 1./pow(r_step,2.);
	      vM[3] = 1./pow(r_step,2.);
	      vM[4] = 1./pow(s_step,2.);

	    }
	}
      else
	{
	  NrOfNonZerosM += 9;
	  
	  if (package == Aztec)
	    {
	      jaM[localrow+1] = jaM[localrow] + 8;
	      jaM[jaM[localrow]] = row - rdim - 1;
	      jaM[jaM[localrow]+1] = row - rdim;
	      jaM[jaM[localrow]+2] = row - rdim + 1;
	      jaM[jaM[localrow]+3] = row - 1;
	      jaM[jaM[localrow]+4] = row + 1;
	      jaM[jaM[localrow]+5] = row + rdim - 1;
	      jaM[jaM[localrow]+6] = row + rdim;
	      jaM[jaM[localrow]+7] = row + rdim + 1;
	      //      cout << " column = " <<   row + rdim + 1 << " i,j = " << i << " , " << j << " ";
	      
	      double coeff11[3][3] = {{0,0,0}, {0,1.,0}, {0,0,0}};
	      aM[localrow] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11);
	      double coeff00[3][3] = {{1.,0,0}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00);
	      double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]+1] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01);
	      double coeff02[3][3] = {{0,0,1.}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]+2] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02);
	      double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	      aM[jaM[localrow]+3] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10);
	      double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	      aM[jaM[localrow]+4]= xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12);
	      double coeff20[3][3] = {{0,0,0}, {0,0,0}, {1.,0,0}};
	      aM[jaM[localrow]+5] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20);
	      double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	      aM[jaM[localrow]+6] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21);
	      double coeff22[3][3] = {{0,0,0}, {0,0,0}, {0,0,1.}};
	      aM[jaM[localrow]+7] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22);
	    }
	  else if (package == PETSc)
	    {
	      mM = 1;
	      //	      idxmM[0] = row;
	      idxmM[0] = (*petscIndexM)(i,j);
	      nM = 9;
	      //	      idxnM[0] = row;
	      //	      idxnM[1] = row - rdim - 1;
	      //	      idxnM[2] = row - rdim;
	      //	      idxnM[3] = row - rdim + 1;
	      //	      idxnM[4] = row - 1;
	      //	      idxnM[5] = row + 1;
	      //	      idxnM[6] = row + rdim - 1;
	      //	      idxnM[7] = row + rdim;
	      //	      idxnM[8] = row + rdim + 1;
	      idxnM[0] = (*petscIndexM)(i,j);
	      idxnM[1] = (*petscIndexM)(i-1,j-1);
	      idxnM[2] = (*petscIndexM)(i,j-1);
	      idxnM[3] = (*petscIndexM)(i+1,j-1);
	      idxnM[4] = (*petscIndexM)(i-1,j);
	      idxnM[5] = (*petscIndexM)(i+1,j);
	      idxnM[6] = (*petscIndexM)(i-1,j+1);
	      idxnM[7] = (*petscIndexM)(i,j+1);
	      idxnM[8] = (*petscIndexM)(i+1,j+1);

	      double coeff11[3][3] = {{0,0,0}, {0,1.,0}, {0,0,0}};
	      vM[0] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11);
	      double coeff00[3][3] = {{1.,0,0}, {0,0,0}, {0,0,0}};
	      vM[1] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00);
	      double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	      vM[2] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01);
	      double coeff02[3][3] = {{0,0,1.}, {0,0,0}, {0,0,0}};
	      vM[3] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02);
	      double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	      vM[4] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10);
	      double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	      vM[5] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12);
	      double coeff20[3][3] = {{0,0,0}, {0,0,0}, {1.,0,0}};
	      vM[6] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20);
	      double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	      vM[7] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21);
	      double coeff22[3][3] = {{0,0,0}, {0,0,0}, {0,0,1.}};
	      vM[8] = xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22)
		+ yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22);
	    }
	}
    }
  else if (operatorType == Grad_Visc_Div || operatorType == Grad_Visc_Div_with_TD_diagonal)
    {
      if (gridType == 1) // cartesian grid
	{
	  NrOfNonZerosM += 5;
	  
	  if (package == Aztec)
	    {
	      jaM[localrow+1] = jaM[localrow] + 4;
	      jaM[jaM[localrow]] = row - rdim;
	      jaM[jaM[localrow]+1] = row - 1;
	      jaM[jaM[localrow]+2] = row + 1;
	      jaM[jaM[localrow]+3] = row + rdim;
	      
	      aM[localrow] = (- 2./pow(s_step,2.) - 2./pow(r_step,2.))*(*viscosity)(i,j);
	      aM[jaM[localrow]] = (1./pow(s_step,2.))*(*viscosity)(i,j);
	      aM[jaM[localrow]+1] = (1./pow(r_step,2.))*(*viscosity)(i,j);
	      aM[jaM[localrow]+2] = (1./pow(r_step,2.))*(*viscosity)(i,j);
	      aM[jaM[localrow]+3] = (1./pow(s_step,2.))*(*viscosity)(i,j);
	    }
	  else if (package == PETSc)
	    {
	      mM = 1;
	      //	      idxmM[0] = row;
	      idxmM[0] = (*petscIndexM)(i,j);
	      nM = 5;
	      //	      idxnM[0] = row;
	      //	      idxnM[1] = row - rdim;
	      //	      idxnM[2] = row - 1;
	      //	      idxnM[3] = row + 1;
	      //	      idxnM[4] = row + rdim;
	      idxnM[0] = (*petscIndexM)(i,j);
	      idxnM[1] = (*petscIndexM)(i,j-1);
	      idxnM[2] = (*petscIndexM)(i-1,j);
	      idxnM[3] = (*petscIndexM)(i+1,j);
	      idxnM[4] = (*petscIndexM)(i,j+1);
	      
	      vM[0] = (- 2./pow(s_step,2.) - 2./pow(r_step,2.))*(*viscosity)(i,j);
	      vM[1] = (1./pow(s_step,2.))*(*viscosity)(i,j);
	      vM[2] = (1./pow(r_step,2.))*(*viscosity)(i,j);
	      vM[3] = (1./pow(r_step,2.))*(*viscosity)(i,j);
	      vM[4] = (1./pow(s_step,2.))*(*viscosity)(i,j);

	    }
	}	      
      else
	{
	  NrOfNonZerosM += 9;
	  
	  if (package == Aztec)
	    {
	      jaM[localrow+1] = jaM[localrow] + 8;
	      jaM[jaM[localrow]] = row - rdim - 1;
	      jaM[jaM[localrow]+1] = row - rdim;
	      jaM[jaM[localrow]+2] = row - rdim + 1;
	      jaM[jaM[localrow]+3] = row - 1;
	      jaM[jaM[localrow]+4] = row + 1;
	      jaM[jaM[localrow]+5] = row + rdim - 1;
	      jaM[jaM[localrow]+6] = row + rdim;
	      jaM[jaM[localrow]+7] = row + rdim + 1;
	      //      cout << " column = " <<   row + rdim + 1 << " i,j = " << i << " , " << j << " ";
	      
	      double coeff11[3][3] = {{0,0,0}, {0,1.,0}, {0,0,0}};
	      aM[localrow] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11)
			     + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11))
		*(*viscosity)(i,j);
	      double coeff00[3][3] = {{1.,0,0}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00)
				 + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00))
		*(*viscosity)(i,j);
	      double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]+1] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01))
		*(*viscosity)(i,j);
	      double coeff02[3][3] = {{0,0,1.}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]+2] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02))
		*(*viscosity)(i,j);
	      double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	      aM[jaM[localrow]+3] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10))
		*(*viscosity)(i,j);
	      double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	      aM[jaM[localrow]+4]= (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12)
				  + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12))
		*(*viscosity)(i,j);
	      double coeff20[3][3] = {{0,0,0}, {0,0,0}, {1.,0,0}};
	      aM[jaM[localrow]+5] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20))
		*(*viscosity)(i,j);
	      double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	      aM[jaM[localrow]+6] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21))
		*(*viscosity)(i,j);
	      double coeff22[3][3] = {{0,0,0}, {0,0,0}, {0,0,1.}};
	      aM[jaM[localrow]+7] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22))
		*(*viscosity)(i,j);
	    }
	  else if (package == PETSc)
	    {
	      mM = 1;
	      //	      idxmM[0] = row;
	      idxmM[0] = (*petscIndexM)(i,j);
	      nM = 9;
	      //	      idxnM[0] = row;
	      //	      idxnM[1] = row - rdim - 1;
	      //	      idxnM[2] = row - rdim;
	      //	      idxnM[3] = row - rdim + 1;
	      //	      idxnM[4] = row - 1;
	      //	      idxnM[5] = row + 1;
	      //	      idxnM[6] = row + rdim - 1;
	      //	      idxnM[7] = row + rdim;
	      //	      idxnM[8] = row + rdim + 1;
	      idxnM[0] = (*petscIndexM)(i,j);
	      idxnM[1] = (*petscIndexM)(i-1,j-1);
	      idxnM[2] = (*petscIndexM)(i,j-1);
	      idxnM[3] = (*petscIndexM)(i+1,j-1);
	      idxnM[4] = (*petscIndexM)(i-1,j);
	      idxnM[5] = (*petscIndexM)(i+1,j);
	      idxnM[6] = (*petscIndexM)(i-1,j+1);
	      idxnM[7] = (*petscIndexM)(i,j+1);
	      idxnM[8] = (*petscIndexM)(i+1,j+1);

	      double coeff11[3][3] = {{0,0,0}, {0,1.,0}, {0,0,0}};
	      vM[0] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11)
			     + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11))
		*(*viscosity)(i,j);
	      double coeff00[3][3] = {{1.,0,0}, {0,0,0}, {0,0,0}};
	      vM[1] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00)
				 + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00))
		*(*viscosity)(i,j);
	      double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	      vM[2] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01))
		*(*viscosity)(i,j);
	      double coeff02[3][3] = {{0,0,1.}, {0,0,0}, {0,0,0}};
	      vM[3] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02))
		*(*viscosity)(i,j);
	      double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	      vM[4] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10))
		*(*viscosity)(i,j);
	      double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	      vM[5] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12)
				  + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12))
		*(*viscosity)(i,j);
	      double coeff20[3][3] = {{0,0,0}, {0,0,0}, {1.,0,0}};
	      vM[6] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20))
		*(*viscosity)(i,j);
	      double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	      vM[7] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21))
		*(*viscosity)(i,j);
	      double coeff22[3][3] = {{0,0,0}, {0,0,0}, {0,0,1.}};
	      vM[8] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22))
		*(*viscosity)(i,j);
	    }
	}
    }
  else if (operatorType == Conv_varDiff || operatorType == Conv_varDiff_with_TD_diagonal)
    {
      if (gridType == 1) // cartesian grid
	{
	  NrOfNonZerosM += 5;
	  
	  if (package == Aztec)
	    {
	      jaM[localrow+1] = jaM[localrow] + 4;
	      jaM[jaM[localrow]] = row - rdim;
	      jaM[jaM[localrow]+1] = row - 1;
	      jaM[jaM[localrow]+2] = row + 1;
	      jaM[jaM[localrow]+3] = row + rdim;
	      
	      aM[localrow] = (- 2./pow(s_step,2.) - 2./pow(r_step,2.))*(*viscosity)(i,j);
	      aM[jaM[localrow]] = (1./pow(s_step,2.))*(*viscosity)(i,j)
		- 0.5*(*convFieldV)(i,j)/s_step;
	      aM[jaM[localrow]+1] = (1./pow(r_step,2.))*(*viscosity)(i,j)
		- 0.5*(*convFieldU)(i,j)/r_step;
	      aM[jaM[localrow]+2] = (1./pow(r_step,2.))*(*viscosity)(i,j)
		+ 0.5*(*convFieldV)(i,j)/r_step;
	      aM[jaM[localrow]+3] = (1./pow(s_step,2.))*(*viscosity)(i,j)
		+ 0.5*(*convFieldU)(i,j)/s_step;
	    }
	  else if (package == PETSc)
	    {
	      mM = 1;
	      //	      idxmM[0] = row;
	      idxmM[0] = (*petscIndexM)(i,j);
	      nM = 5;
	      //	      idxnM[0] = row;
	      //	      idxnM[1] = row - rdim;
	      //	      idxnM[2] = row - 1;
	      //	      idxnM[3] = row + 1;
	      //	      idxnM[4] = row + rdim;
	      idxnM[0] = (*petscIndexM)(i,j);
	      idxnM[1] = (*petscIndexM)(i,j-1);
	      idxnM[2] = (*petscIndexM)(i-1,j);
	      idxnM[3] = (*petscIndexM)(i+1,j);
	      idxnM[4] = (*petscIndexM)(i,j+1);
	      
	      vM[0] = (- 2./pow(s_step,2.) - 2./pow(r_step,2.))*(*viscosity)(i,j);
	      vM[1] = (1./pow(s_step,2.))*(*viscosity)(i,j)
		- 0.5*(*convFieldV)(i,j)/s_step;
	      vM[2] = (1./pow(r_step,2.))*(*viscosity)(i,j)
		- 0.5*(*convFieldU)(i,j)/r_step;
	      vM[3] = (1./pow(r_step,2.))*(*viscosity)(i,j)
		+ 0.5*(*convFieldV)(i,j)/r_step;
	      vM[4] = (1./pow(s_step,2.))*(*viscosity)(i,j)
		+ 0.5*(*convFieldU)(i,j)/s_step;

	    }
	}	      
      else
	{
	  NrOfNonZerosM += 9;
	  
	  if (package == Aztec)
	    {
	      jaM[localrow+1] = jaM[localrow] + 8;
	      jaM[jaM[localrow]] = row - rdim - 1;
	      jaM[jaM[localrow]+1] = row - rdim;
	      jaM[jaM[localrow]+2] = row - rdim + 1;
	      jaM[jaM[localrow]+3] = row - 1;
	      jaM[jaM[localrow]+4] = row + 1;
	      jaM[jaM[localrow]+5] = row + rdim - 1;
	      jaM[jaM[localrow]+6] = row + rdim;
	      jaM[jaM[localrow]+7] = row + rdim + 1;
	      //      cout << " column = " <<   row + rdim + 1 << " i,j = " << i << " , " << j << " ";
	      
	      double coeff11[3][3] = {{0,0,0}, {0,1.,0}, {0,0,0}};
	      aM[localrow] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11)
			     + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff11,convFieldU)*(*convFieldU)(i,j)
			     + upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff11,convFieldV)
		*(*convFieldV)(i,j);
	      
	      double coeff00[3][3] = {{1.,0,0}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00)
				 + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff00,convFieldU)*(*convFieldU)(i,j)
				 + upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff00,convFieldV)
		*(*convFieldV)(i,j);
	      double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]+1] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01,convFieldU)*(*convFieldU)(i,j)
				   + upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01,convFieldV)
		*(*convFieldV)(i,j);
	      double coeff02[3][3] = {{0,0,1.}, {0,0,0}, {0,0,0}};
	      aM[jaM[localrow]+2] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff02,convFieldU)*(*convFieldU)(i,j)
				   + upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff02,convFieldV)
		*(*convFieldV)(i,j);
	      double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	      aM[jaM[localrow]+3] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10,convFieldU)*(*convFieldU)(i,j)
				   + upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10,convFieldV)
		*(*convFieldV)(i,j);
	      double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	      aM[jaM[localrow]+4]= (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12)
				  + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12,convFieldU)*(*convFieldU)(i,j)
				  + upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12,convFieldV)
		*(*convFieldV)(i,j);
	      double coeff20[3][3] = {{0,0,0}, {0,0,0}, {1.,0,0}};
	      aM[jaM[localrow]+5] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff20,convFieldU)*(*convFieldU)(i,j)
				   + upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff20,convFieldV)
		*(*convFieldV)(i,j);
	      double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	      aM[jaM[localrow]+6] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21,convFieldU)*(*convFieldU)(i,j)
				   + upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21,convFieldV)
		*(*convFieldV)(i,j);
	      double coeff22[3][3] = {{0,0,0}, {0,0,0}, {0,0,1.}};
	      aM[jaM[localrow]+7] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff22,convFieldU)*(*convFieldU)(i,j)
				   + upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff22,convFieldV)
		*(*convFieldV)(i,j);
	    }
	  else if (package == PETSc)
	    {
	      mM = 1;
	      //	      idxmM[0] = row;
	      idxmM[0] = (*petscIndexM)(i,j);
	      nM = 9;
	      //	      idxnM[0] = row;
	      //	      idxnM[1] = row - rdim - 1;
	      //	      idxnM[2] = row - rdim;
	      //	      idxnM[3] = row - rdim + 1;
	      //	      idxnM[4] = row - 1;
	      //	      idxnM[5] = row + 1;
	      //	      idxnM[6] = row + rdim - 1;
	      //	      idxnM[7] = row + rdim;
	      //	      idxnM[8] = row + rdim + 1;
	      idxnM[0] = (*petscIndexM)(i,j);
	      idxnM[1] = (*petscIndexM)(i-1,j-1);
	      idxnM[2] = (*petscIndexM)(i,j-1);
	      idxnM[3] = (*petscIndexM)(i+1,j-1);
	      idxnM[4] = (*petscIndexM)(i-1,j);
	      idxnM[5] = (*petscIndexM)(i+1,j);
	      idxnM[6] = (*petscIndexM)(i-1,j+1);
	      idxnM[7] = (*petscIndexM)(i,j+1);
	      idxnM[8] = (*petscIndexM)(i+1,j+1);

	      double coeff11[3][3] = {{0,0,0}, {0,1.,0}, {0,0,0}};
	      vM[0] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11)
		      + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff11))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff11,convFieldU)*(*convFieldU)(i,j)
		+ upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff11,convFieldV)*(*convFieldV)(i,j);
	      double coeff00[3][3] = {{1.,0,0}, {0,0,0}, {0,0,0}};
	      vM[1] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00)
				 + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff00))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff00,convFieldU)*(*convFieldU)(i,j)
		+ upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff00,convFieldV)*(*convFieldV)(i,j);
	      double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	      vM[2] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff01))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01,convFieldU)*(*convFieldU)(i,j)
		+ upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01,convFieldV)*(*convFieldV)(i,j);
	      double coeff02[3][3] = {{0,0,1.}, {0,0,0}, {0,0,0}};
	      vM[3] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff02))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff02,convFieldU)*(*convFieldU)(i,j)
		+ upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff02,convFieldV)*(*convFieldV)(i,j);
	      double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	      vM[4] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff10))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10,convFieldU)*(*convFieldU)(i,j)
		+ upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10,convFieldV)*(*convFieldV)(i,j);
	      double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	      vM[5] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12)
				  + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff12))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12,convFieldU)*(*convFieldU)(i,j)
		+ upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12,convFieldV)*(*convFieldV)(i,j);
	      double coeff20[3][3] = {{0,0,0}, {0,0,0}, {1.,0,0}};
	      vM[6] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff20))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff20,convFieldU)*(*convFieldU)(i,j)
		+ upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff20,convFieldV)*(*convFieldV)(i,j);
	      double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	      vM[7] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff21))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21,convFieldU)*(*convFieldU)(i,j)
		+ upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21,convFieldV)*(*convFieldV)(i,j);
	      double coeff22[3][3] = {{0,0,0}, {0,0,0}, {0,0,1.}};
	      vM[8] = (xxCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22)
				   + yyCoeff(i,j,Xr,Xs,Yr,Ys,Xrr,Xss,Yrr,Yss,sqrtOfG,r_step,s_step,coeff22))
		*(*viscosity)(i,j)
		+ upWindxCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff22,convFieldU)*(*convFieldU)(i,j)
		+ upWindyCoeff(i,j,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff22,convFieldV)*(*convFieldV)(i,j);
	    }
	}
    }
}

void
DMSRMatrix::insertNeumann(int row, int localrow, int side, int rdim, int i, int j, int gridType, doubleSerialArray *Xr, doubleSerialArray *Xs, doubleSerialArray *Yr, doubleSerialArray *Ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, sparseSolver package)
{
  if (gridType == 1) // cartesian grid
  {
    if (package == Aztec)
      {
	jaM[localrow+1] = jaM[localrow] + 1;
	
	NrOfNonZerosM += 2;
	if (side==0) // low r
	  {
	    jaM[jaM[localrow]] = row + 2;
	    aM[localrow] = -1./(2*r_step);
	    aM[jaM[localrow]] = 1./(2*r_step);
	  }
	else if (side==1) // high r
	  {
	    jaM[jaM[localrow]] = row - 2;
	    aM[localrow] = 1./(2*r_step);
	    aM[jaM[localrow]] = -1./(2*r_step);
	  }
	else if (side==2) // low s
	  {
	    jaM[jaM[localrow]] = row + 2*rdim;
	    aM[localrow] = -1./(2*s_step);
	    aM[jaM[localrow]] = 1./(2*s_step);
	  }
	else if (side==3) // high s
	  {
	    jaM[jaM[localrow]] = row - 2*rdim;
	    aM[localrow] = 1./(2*s_step);
	    aM[jaM[localrow]] = -1./(2*s_step);
	  }
      }
    else if (package == PETSc)
      {
	NrOfNonZerosM += 2;

	mM = 1;
	//	idxmM[0] = row;
	nM = 2;
	//	idxnM[0] = row;
	if (side==0) // low r
	  {
	    idxmM[0] = (*petscIndexM)(i-1,j);
	    idxnM[0] = (*petscIndexM)(i-1,j);
	    //	    idxnM[1] = row + 2;
	    idxnM[1] = (*petscIndexM)(i+1,j);
	    vM[0] = -1./(2*r_step);
	    vM[1] = 1./(2*r_step);
	  }
	else if (side==1) // high r
	  {
	    idxmM[0] = (*petscIndexM)(i+1,j);
	    idxnM[0] = (*petscIndexM)(i+1,j);
	    //	    idxnM[1] = row - 2;
	    idxnM[1] = (*petscIndexM)(i-1,j);
	    vM[0] = 1./(2*r_step);
	    vM[1] = -1./(2*r_step);
	  }
	else if (side==2) // low s
	  {
	    idxmM[0] = (*petscIndexM)(i,j-1);
	    idxnM[0] = (*petscIndexM)(i,j-1);
	    //	    idxnM[1] = row + 2*rdim;
	    idxnM[1] = (*petscIndexM)(i,j+1);
	    vM[0] = -1./(2*s_step);
	    vM[1] = 1./(2*s_step);
	  }
	else if (side==3) // high s
	  {
	    idxmM[0] = (*petscIndexM)(i,j+1);
	    idxnM[0] = (*petscIndexM)(i,j+1);
	    //	    idxnM[1] = row - 2*rdim;
	    idxnM[1] = (*petscIndexM)(i,j-1);
	    vM[0] = 1./(2*s_step);
	    vM[1] = -1./(2*s_step);
	  }
      }
  }
  else
  {
    if (package == Aztec)
      {
	jaM[localrow+1] = jaM[localrow] + 3;
	
	NrOfNonZerosM += 4;
	if (side==0) // low r
	  {
	    jaM[jaM[localrow]] = row - rdim + 1;
	    jaM[jaM[localrow]+1] = row + 2;
	    jaM[jaM[localrow]+2] = row + rdim + 1;
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    aM[localrow] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    aM[jaM[localrow]] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    aM[jaM[localrow]+1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    aM[jaM[localrow]+2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	  }
	else if (side==1) // high r
	  {
	    jaM[jaM[localrow]] = row - rdim - 1;
	    jaM[jaM[localrow]+1] = row - 2;
	    jaM[jaM[localrow]+2] = row + rdim - 1;
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    aM[localrow] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    aM[jaM[localrow]] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    aM[jaM[localrow]+1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    aM[jaM[localrow]+2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	  }
	else if (side==2) // low s
	  {
	    jaM[jaM[localrow]] = row + rdim - 1;
	    jaM[jaM[localrow]+1] = row + rdim + 1;
	    jaM[jaM[localrow]+2] = row + 2*rdim;
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    aM[localrow] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    aM[jaM[localrow]] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    aM[jaM[localrow]+1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    aM[jaM[localrow]+2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	  }
	else if (side==3) // high s
	  {
	    jaM[jaM[localrow]] = row - 2*rdim;
	    jaM[jaM[localrow]+1] = row - rdim - 1;
	    jaM[jaM[localrow]+2] = row - rdim + 1;
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    aM[localrow] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    aM[jaM[localrow]] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    aM[jaM[localrow]+1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    aM[jaM[localrow]+2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	  }
      }
    else if (package == PETSc)
      {
	NrOfNonZerosM += 4;

	mM = 1;
	//	idxmM[0] = row;
	nM = 4;
	//	idxnM[0] = row;
	if (side==0) // low r
	  {
	    idxmM[0] = (*petscIndexM)(i-1,j);
	    idxnM[0] = (*petscIndexM)(i-1,j);
	    //	    idxnM[1] = row - rdim + 1;
	    //	    idxnM[2] = row + 2;
	    //	    idxnM[3] = row + rdim + 1;
	    idxnM[1] = (*petscIndexM)(i,j-1);
	    idxnM[2] = (*petscIndexM)(i+1,j);
	    idxnM[3] = (*petscIndexM)(i,j+1);

	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    vM[0] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    vM[1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    vM[2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    vM[3] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	  }
	else if (side==1) // high r
	  {
	    idxmM[0] = (*petscIndexM)(i+1,j);
	    idxnM[0] = (*petscIndexM)(i+1,j);
	    //	    idxnM[1] = row - rdim - 1;
	    //	    idxnM[2] = row - 2;
	    //	    idxnM[3] = row + rdim - 1;
	    idxnM[1] = (*petscIndexM)(i,j-1);
	    idxnM[2] = (*petscIndexM)(i-1,j);
	    idxnM[3] = (*petscIndexM)(i,j+1);

	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    vM[0] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    vM[1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    vM[2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    vM[3] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	  }
	else if (side==2) // low s
	  {
	    idxmM[0] = (*petscIndexM)(i,j-1);
	    idxnM[0] = (*petscIndexM)(i,j-1);
	    //	    idxnM[1] = row + rdim - 1;
	    //	    idxnM[2] = row + rdim + 1;
	    //	    idxnM[3] = row + 2*rdim;
	    idxnM[1] = (*petscIndexM)(i-1,j);
	    idxnM[2] = (*petscIndexM)(i+1,j);
	    idxnM[3] = (*petscIndexM)(i,j+1);

	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    vM[0] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    vM[1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    vM[2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    vM[3] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	  }
	else if (side==3) // high s
	  {
	    idxmM[0] = (*petscIndexM)(i,j+1);
	    idxnM[0] = (*petscIndexM)(i,j+1);
	    //	    idxnM[1] = row - 2*rdim;
	    //	    idxnM[2] = row - rdim - 1;
	    //	    idxnM[3] = row - rdim + 1;
	    idxnM[1] = (*petscIndexM)(i,j-1);
	    idxnM[2] = (*petscIndexM)(i-1,j);
	    idxnM[3] = (*petscIndexM)(i+1,j);

	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    vM[0] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    vM[1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    vM[2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    vM[3] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	  }
      }
  }
}

void
DMSRMatrix::insertRobin(int row, int localrow, int side, int rdim, int i, int j, int gridType, doubleSerialArray *Xr, doubleSerialArray *Xs, doubleSerialArray *Yr, doubleSerialArray *Ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, sparseSolver package)
{
  double RobConst = 0.1;
  if (gridType == 1) // cartesian grid
  {
    if (package == Aztec)
      {
	jaM[localrow+1] = jaM[localrow] + 2;
	
	NrOfNonZerosM += 3;
	if (side==0) // low r
	  {
	    jaM[jaM[localrow]] = row + 1;
	    jaM[jaM[localrow]+1] = row + 2;
	    aM[localrow] = -1./(2*r_step);
	    aM[jaM[localrow]] = RobConst;
	    aM[jaM[localrow]+1] = 1./(2*r_step);
	  }
	else if (side==1) // high r
	  {
	    jaM[jaM[localrow]] = row - 1;
	    jaM[jaM[localrow]+1] = row - 2;
	    aM[localrow] = 1./(2*r_step);
	    aM[jaM[localrow]] = RobConst;
	    aM[jaM[localrow]+1] = -1./(2*r_step);
	  }
	else if (side==2) // low s
	  {
	    jaM[jaM[localrow]] = row + rdim;
	    jaM[jaM[localrow]+1] = row + 2*rdim;
	    aM[localrow] = -1./(2*s_step);
	    aM[jaM[localrow]] = RobConst;
	    aM[jaM[localrow]+1] = 1./(2*s_step);
	  }
	else if (side==3) // high s
	  {
	    jaM[jaM[localrow]] = row - rdim;
	    jaM[jaM[localrow]+1] = row - 2*rdim;
	    aM[localrow] = 1./(2*s_step);
	    aM[jaM[localrow]] = RobConst;
	    aM[jaM[localrow]+1] = -1./(2*s_step);
	  }
      }
    else if (package == PETSc)
      {
	NrOfNonZerosM += 3;
	
	mM = 1;
	//	idxmM[0] = row;
	nM = 3;
	//	idxnM[0] = row;
	if (side==0) // low r
	  {
	    idxmM[0] = (*petscIndexM)(i-1,j);
	    idxnM[0] = (*petscIndexM)(i-1,j);
	    //	    idxnM[1] = row + 1;
	    //	    idxnM[2] = row + 2;
	    idxnM[1] = (*petscIndexM)(i,j);
	    idxnM[2] = (*petscIndexM)(i+1,j);

	    vM[0] = -1./(2*r_step);
	    vM[1] = RobConst;
	    vM[2] = 1./(2*r_step);
	  }
	else if (side==1) // high r
	  {
	    idxmM[0] = (*petscIndexM)(i+1,j);
	    idxnM[0] = (*petscIndexM)(i+1,j);
	    //	    idxnM[1] = row - 1;
	    //	    idxnM[2] = row - 2;
	    idxnM[1] = (*petscIndexM)(i,j);
	    idxnM[2] = (*petscIndexM)(i-1,j);

	    vM[0] = 1./(2*r_step);
	    vM[1] = RobConst;
	    vM[2] = -1./(2*r_step);
	  }
	else if (side==2) // low s
	  {
	    idxmM[0] = (*petscIndexM)(i,j-1);
	    idxnM[0] = (*petscIndexM)(i,j-1);
	    //	    idxnM[1] = row + rdim;
	    //	    idxnM[2] = row + 2*rdim;
	    idxnM[1] = (*petscIndexM)(i,j);
	    idxnM[2] = (*petscIndexM)(i,j+1);

	    vM[0] = -1./(2*s_step);
	    vM[1] = RobConst;
	    vM[2] = 1./(2*s_step);
	  }
	else if (side==3) // high s
	  {
	    idxmM[0] = (*petscIndexM)(i,j+1);
	    idxnM[0] = (*petscIndexM)(i,j+1);
	    //	    idxnM[1] = row - rdim;
	    //	    idxnM[2] = row - 2*rdim;
	    idxnM[1] = (*petscIndexM)(i,j);
	    idxnM[2] = (*petscIndexM)(i,j-1);

	    vM[0] = 1./(2*s_step);
	    vM[1] = RobConst;
	    vM[2] = -1./(2*s_step);
	  }
      }
  }
  else
  {
    if (package == Aztec)
      {
	jaM[localrow+1] = jaM[localrow] + 4;
	
	NrOfNonZerosM += 5;
	if (side==0) // low r
	  {
	    jaM[jaM[localrow]] = row - rdim + 1;
	    jaM[jaM[localrow]+1] = row + 2;
	    jaM[jaM[localrow]+2] = row + rdim + 1;
	    jaM[jaM[localrow]+3] = row + 1;
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    aM[localrow] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    aM[jaM[localrow]] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    aM[jaM[localrow]+1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    aM[jaM[localrow]+2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    aM[jaM[localrow]+3] = RobConst;
	  }
	else if (side==1) // high r
	  {
	    jaM[jaM[localrow]] = row - rdim - 1;
	    jaM[jaM[localrow]+1] = row - 2;
	    jaM[jaM[localrow]+2] = row + rdim - 1;
	    jaM[jaM[localrow]+3] = row - 1;
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    aM[localrow] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    aM[jaM[localrow]] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    aM[jaM[localrow]+1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    aM[jaM[localrow]+2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    aM[jaM[localrow]+3] = RobConst;
	  }
	else if (side==2) // low s
	  {
	    jaM[jaM[localrow]] = row + rdim - 1;
	    jaM[jaM[localrow]+1] = row + rdim + 1;
	    jaM[jaM[localrow]+2] = row + 2*rdim;
	    jaM[jaM[localrow]+3] = row + rdim;
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    aM[localrow] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    aM[jaM[localrow]] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    aM[jaM[localrow]+1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    aM[jaM[localrow]+2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    aM[jaM[localrow]+3] = RobConst; 
	  }
	else if (side==3) // high s
	  {
	    jaM[jaM[localrow]] = row - 2*rdim;
	    jaM[jaM[localrow]+1] = row - rdim - 1;
	    jaM[jaM[localrow]+2] = row - rdim + 1;
	    jaM[jaM[localrow]+3] = row - rdim;
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    aM[localrow] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    aM[jaM[localrow]] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    aM[jaM[localrow]+1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    aM[jaM[localrow]+2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    aM[jaM[localrow]+3] = RobConst;
	  }
      }
    else if (package == PETSc)
      {
	NrOfNonZerosM += 5;
	
	mM = 1;
	//	idxmM[0] = row;
	nM = 5;
	//	idxnM[0] = row;

	if (side==0) // low r
	  {
	    idxmM[0] = (*petscIndexM)(i-1,j);
	    idxnM[0] = (*petscIndexM)(i-1,j);
	    //	    idxnM[1] = row - rdim + 1;
	    //	    idxnM[2] = row + 2;
	    //	    idxnM[3] = row + rdim + 1;
	    //	    idxnM[4] = row + 1;
	    idxnM[1] = (*petscIndexM)(i,j-1);
	    idxnM[2] = (*petscIndexM)(i+1,j);
	    idxnM[3] = (*petscIndexM)(i,j+1);
	    idxnM[4] = (*petscIndexM)(i,j);

	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    vM[0] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    vM[1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    vM[2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    vM[3] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    vM[4] = RobConst;
	  }
	else if (side==1) // high r
	  {
	    idxmM[0] = (*petscIndexM)(i+1,j);
	    idxnM[0] = (*petscIndexM)(i+1,j);
	    //	    idxnM[1] = row - rdim - 1;
	    //	    idxnM[2] = row - 2;
	    //	    idxnM[3] = row + rdim - 1;
	    //	    idxnM[4] = row - 1;
	    idxnM[1] = (*petscIndexM)(i,j-1);
	    idxnM[2] = (*petscIndexM)(i-1,j);
	    idxnM[3] = (*petscIndexM)(i,j+1);
	    idxnM[4] = (*petscIndexM)(i,j);

	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    vM[0] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    vM[1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    vM[2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    vM[3] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    vM[4] = RobConst;
	  }
	else if (side==2) // low s
	  {
	    idxmM[0] = (*petscIndexM)(i,j-1);
	    idxnM[0] = (*petscIndexM)(i,j-1);
	    //	    idxnM[1] = row + rdim - 1;
	    //	    idxnM[2] = row + rdim + 1;
	    //	    idxnM[3] = row + 2*rdim;
	    //	    idxnM[4] = row + rdim;
	    idxnM[1] = (*petscIndexM)(i-1,j);
	    idxnM[2] = (*petscIndexM)(i+1,j);
	    idxnM[3] = (*petscIndexM)(i,j+1);
	    idxnM[4] = (*petscIndexM)(i,j);

	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    vM[0] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    vM[1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    vM[2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    vM[3] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    vM[4] = RobConst; 
	  }
	else if (side==3) // high s
	  {
	    idxmM[0] = (*petscIndexM)(i,j+1);
	    idxnM[0] = (*petscIndexM)(i,j+1);
	    //	    idxnM[1] = row - 2*rdim;
	    //	    idxnM[2] = row - rdim - 1;
	    //	    idxnM[3] = row - rdim + 1;
	    //	    idxnM[4] = row - rdim;
	    idxnM[1] = (*petscIndexM)(i,j-1);
	    idxnM[2] = (*petscIndexM)(i-1,j);
	    idxnM[3] = (*petscIndexM)(i-1,j);
	    idxnM[4] = (*petscIndexM)(i,j);

	    double coeff21[3][3] = {{0,0,0}, {0,0,0}, {0,1.,0}};
	    vM[0] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff21);
	    double coeff01[3][3] = {{0,1.,0}, {0,0,0}, {0,0,0}};
	    vM[1] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff01);
	    double coeff10[3][3] = {{0,0,0}, {1.,0,0}, {0,0,0}};
	    vM[2] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff10);
	    double coeff12[3][3] = {{0,0,0}, {0,0,1.}, {0,0,0}};
	    vM[3] = NeumannCoeff(i,j,side,Xr,Xs,Yr,Ys,sqrtOfG,r_step,s_step,coeff12);
	    vM[4] = RobConst;
	  }
      }
  }
}

void
DMSRMatrix::insertPeriod(int row, int localrow, int side, int rdim, int sdim, sparseSolver package, int grid, int i, int j)
{
  if (package == Aztec)
    {
      jaM[localrow+1] = jaM[localrow] + 1;
      
      NrOfNonZerosM += 2;
      if (side==0) // low r
	{
	  jaM[jaM[localrow]] = row + rdim - 2;
	  aM[localrow] = 1.;
	  aM[jaM[localrow]] = -1.;
	}
      else if (side==1) // high r
	{
	  jaM[jaM[localrow]] = row - rdim + 2;
	  aM[localrow] = 1.;
	  aM[jaM[localrow]] = -1.;
	}
      else if (side==2) // low s
	{
	  jaM[jaM[localrow]] = row + (sdim - 2)*rdim;
	  aM[localrow] = 1.;
	  aM[jaM[localrow]] = -1.;
	}
      else if (side==3) // high s
	{
	  jaM[jaM[localrow]] = row - (sdim - 2)*rdim;
	  aM[localrow] = 1.;
	  aM[jaM[localrow]] = -1.;
	}
    }
  else if (package == PETSc)
    {
      int iProc, start;
      NrOfNonZerosM += 2;

      mM = 1;
      //      idxmM[0] = row;
      idxmM[0] = (*petscIndexM)(i,j);
      nM = 2;
      //      idxnM[0] = row;
      idxnM[0] = (*petscIndexM)(i,j);
      if (side==0) // low r
	{
	  //	  idxnM[1] = row + rdim - 2;
	  iProc = ownerProcess(i+rdim-2,j,grid);
	  int lRDim =  petscLocalRHiM[grid](iProc) - petscLocalRLowM[grid](iProc) + 1;
	  start = petscStartIndexM[grid](iProc)
	    + lRDim*(j-petscLocalSLowM[grid](iProc))
	    + (i+rdim-2-petscLocalRLowM[grid](iProc));
	  idxnM[1] = start;
	  vM[0] = 1.;
	  vM[1] = -1.;
	}
      else if (side==1) // high r
	{
	  //	  idxnM[1] = row - rdim + 2;
	  iProc = ownerProcess(i-rdim+2,j,grid);
	  int lRDim =  petscLocalRHiM[grid](iProc) - petscLocalRLowM[grid](iProc) + 1;
	  start = petscStartIndexM[grid](iProc)
	    + lRDim*(j-petscLocalSLowM[grid](iProc))
	    + (i-rdim+2-petscLocalRLowM[grid](iProc));
	  idxnM[1] = start;
	  vM[0] = 1.;
	  vM[1] = -1.;
	}
      else if (side==2) // low s
	{
	  //	  idxnM[1] = row + (sdim - 2)*rdim;
	  iProc = ownerProcess(i,j+sdim-2,grid);
	  int lRDim =  petscLocalRHiM[grid](iProc) - petscLocalRLowM[grid](iProc) + 1;
	  start = petscStartIndexM[grid](iProc)
	    + lRDim*(j+sdim-2-petscLocalSLowM[grid](iProc))
	    + (i-petscLocalRLowM[grid](iProc));
	  idxnM[1] = start;
	  vM[0] = 1.;
	  vM[1] = -1.;
	}
      else if (side==3) // high s
	{
	  //	  idxnM[1] = row - (sdim - 2)*rdim;
	  iProc = ownerProcess(i,j-sdim+2,grid);
	  int lRDim =  petscLocalRHiM[grid](iProc) - petscLocalRLowM[grid](iProc) + 1;
	  start = petscStartIndexM[grid](iProc)
	    + lRDim*(j-sdim+2-petscLocalSLowM[grid](iProc))
	    + (i-petscLocalRLowM[grid](iProc));
	  idxnM[1] = start;
	  vM[0] = 1.;
	  vM[1] = -1.;
	}
    }
}

void
DMSRMatrix::insertCorner(int row, int localrow, int rdim, int corner, sparseSolver package, int i, int j)
{
  if (package == Aztec)
    {
      jaM[localrow+1] = jaM[localrow] + 3;
      
      NrOfNonZerosM += 4;
      if (corner==0) // low r , low s
	{
	  jaM[jaM[localrow]] = row + rdim + 1;
	  jaM[jaM[localrow]+1] = row + 2*rdim + 2;
	  jaM[jaM[localrow]+2] = row + 3*rdim + 3;
	  
	  aM[localrow] = 1.;
	  aM[jaM[localrow]] = -3;
	  aM[jaM[localrow]+1] = 3;
	  aM[jaM[localrow]+2] = -1.;
	}
      else if (corner==1) // low r , high s
	{
	  jaM[jaM[localrow]] = row - rdim + 1;
	  jaM[jaM[localrow]+1] = row - 2*rdim + 2;
	  jaM[jaM[localrow]+2] = row - 3*rdim + 3;
	  
	  aM[localrow] = 1.;
	  aM[jaM[localrow]] = -3;
	  aM[jaM[localrow]+1] = 3;
	  aM[jaM[localrow]+2] = -1.;
	}
      else if (corner==2) // high r , low s
	{
	  jaM[jaM[localrow]] = row + rdim - 1;
	  jaM[jaM[localrow]+1] = row + 2*rdim - 2;
	  jaM[jaM[localrow]+2] = row + 3*rdim - 3;
	  
	  aM[localrow] = 1.;
	  aM[jaM[localrow]] = -3;
	  aM[jaM[localrow]+1] = 3;
	  aM[jaM[localrow]+2] = -1.;
	}
      else if (corner==3) // high r , high s
	{
	  jaM[jaM[localrow]] = row - rdim - 1;
	  jaM[jaM[localrow]+1] = row - 2*rdim - 2;
	  jaM[jaM[localrow]+2] = row - 3*rdim - 3;
	  
	  aM[localrow] = 1.;
	  aM[jaM[localrow]] = -3;
	  aM[jaM[localrow]+1] = 3;
	  aM[jaM[localrow]+2] = -1.;
	}
    }
  else if (package == PETSc)
    {
      NrOfNonZerosM += 4;
     
      mM = 1;
      //      idxmM[0] = row;
      idxmM[0] = (*petscIndexM)(i,j);
      nM = 4;
      //      idxnM[0] = row;
      idxnM[0] = (*petscIndexM)(i,j);
      if (corner==0) // low r , low s
	{
	  //	  idxnM[1] = row + rdim + 1;
	  //	  idxnM[2] = row + 2*rdim + 2;
	  //	  idxnM[3] = row + 3*rdim + 3;
	  idxnM[1] = (*petscIndexM)(i+1,j+1);
	  idxnM[2] = (*petscIndexM)(i+2,j+2);
	  idxnM[3] = (*petscIndexM)(i+3,j+3);

	  vM[0] = 1.;
	  vM[1] = -3;
	  vM[2] = 3;
	  vM[3] = -1.;
	}
      else if (corner==1) // low r , high s
	{
	  //	  idxnM[1] = row - rdim + 1;
	  //	  idxnM[2] = row - 2*rdim + 2;
	  //	  idxnM[3] = row - 3*rdim + 3;
	  idxnM[1] = (*petscIndexM)(i+1,j-1);
	  idxnM[2] = (*petscIndexM)(i+2,j-2);
	  idxnM[3] = (*petscIndexM)(i+3,j-3);
	  
	  vM[0] = 1.;
	  vM[1] = -3;
	  vM[2] = 3;
	  vM[3] = -1.;
	}
      else if (corner==2) // high r , low s
	{
	  //	  idxnM[1] = row + rdim - 1;
	  //	  idxnM[2] = row + 2*rdim - 2;
	  //	  idxnM[3] = row + 3*rdim - 3;
	  idxnM[1] = (*petscIndexM)(i-1,j+1);
	  idxnM[2] = (*petscIndexM)(i-2,j+2);
	  idxnM[3] = (*petscIndexM)(i-3,j+3);
	  
	  vM[0] = 1.;
	  vM[1] = -3;
	  vM[2] = 3;
	  vM[3] = -1.;
	}
      else if (corner==3) // high r , high s
	{
	  //	  idxnM[1] = row - rdim - 1;
	  //	  idxnM[2] = row - 2*rdim - 2;
	  //	  idxnM[3] = row - 3*rdim - 3;
	  idxnM[1] = (*petscIndexM)(i-1,j-1);
	  idxnM[2] = (*petscIndexM)(i-2,j-2);
	  idxnM[3] = (*petscIndexM)(i-3,j-3);
	  
	  vM[0] = 1.;
	  vM[1] = -3;
	  vM[2] = 3;
	  vM[3]= -1.;
	}
    }
}

void
DMSRMatrix::insertInterpolationInfo(int localrow, int *update, int grid, const int donorGrid, int i, int j, CompositeGrid *Kimera, interpolationType iType, sparseSolver package)
{
  if (iType == cubic)
    {
      jaM[localrow+1] = jaM[localrow] + 16;
      aM[localrow] = -1.;
      
      NrOfNonZerosM += 17;
      
      int index = Kimera->getInterpIndex(i,j,grid);
      int i_loc = Kimera->My_i_LM[grid](index);
      int j_loc = Kimera->My_j_LM[grid](index);
      
      int start = 0;
      int offset = 0;
      
      Index all;
      
      for (int k=0; k<donorGrid; k++)
	{
	  start += Kimera->rdimM[k] * Kimera->sdimM[k];
	}
      for (int g=0; g<grid; g++)
	offset += Kimera->myIntPointsM[g].getLength(0);
      
      //  cout << " offset = " << offset << " index = " << index << endl;
      
      jaM[jaM[localrow]] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc;
      jaM[jaM[localrow]+1] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 1;
      jaM[jaM[localrow]+2] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2;
      jaM[jaM[localrow]+3] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 3;
      jaM[jaM[localrow]+4] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid];
      jaM[jaM[localrow]+5] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid] + 1;
      jaM[jaM[localrow]+6] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid] + 2;
      jaM[jaM[localrow]+7] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid] + 3;
      jaM[jaM[localrow]+8] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]);
      jaM[jaM[localrow]+9] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]) + 1;
      jaM[jaM[localrow]+10] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]) + 2;
      jaM[jaM[localrow]+11] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]) + 3;
      jaM[jaM[localrow]+12] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 3*(Kimera->rdimM[donorGrid]);
      jaM[jaM[localrow]+13] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 3*(Kimera->rdimM[donorGrid]) + 1;
      jaM[jaM[localrow]+14] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 3*(Kimera->rdimM[donorGrid]) + 2;
      jaM[jaM[localrow]+15] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 3*(Kimera->rdimM[donorGrid]) + 3;

      aM[jaM[localrow]] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(4,offset+index);
      aM[jaM[localrow]+1] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(4,offset+index);
      aM[jaM[localrow]+2] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(4,offset+index);
      aM[jaM[localrow]+3] = CoeffReceiveBufferM(3,offset+index) * CoeffReceiveBufferM(4,offset+index);
      aM[jaM[localrow]+4] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(5,offset+index);
      aM[jaM[localrow]+5] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(5,offset+index);
      aM[jaM[localrow]+6] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(5,offset+index);
      aM[jaM[localrow]+7] = CoeffReceiveBufferM(3,offset+index) * CoeffReceiveBufferM(5,offset+index);
      aM[jaM[localrow]+8] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(6,offset+index);
      aM[jaM[localrow]+9] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(6,offset+index);
      aM[jaM[localrow]+10] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(6,offset+index);
      aM[jaM[localrow]+11] = CoeffReceiveBufferM(3,offset+index) * CoeffReceiveBufferM(6,offset+index);
      aM[jaM[localrow]+12] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(7,offset+index);
      aM[jaM[localrow]+13] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(7,offset+index);
      aM[jaM[localrow]+14] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(7,offset+index);
      aM[jaM[localrow]+15] = CoeffReceiveBufferM(3,offset+index) * CoeffReceiveBufferM(7,offset+index);
    }
  else if (iType == quadratic)
    {
      if (package == Aztec)
	{
	  jaM[localrow+1] = jaM[localrow] + 9;
	  aM[localrow] = -1.;
	}
      
      NrOfNonZerosM += 10;
      
      int index = Kimera->getInterpIndex(i,j,grid);
      int i_loc = Kimera->My_i_LM[grid](index);
      int j_loc = Kimera->My_j_LM[grid](index);
      
      int start = 0;
      int offset = 0;
      
      Index all;
      
      for (int k=0; k<donorGrid; k++)
	{
	  start += Kimera->rdimM[k] * Kimera->sdimM[k];
	}
      for (int g=0; g<grid; g++)
	offset += Kimera->myIntPointsM[g].getLength(0);
      
      //  cout << " offset = " << offset << " index = " << index << endl;
      
      if (package == Aztec)
	{
	  jaM[jaM[localrow]] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc;
	  jaM[jaM[localrow]+1] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 1;
	  jaM[jaM[localrow]+2] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2;
	  jaM[jaM[localrow]+3] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid];
	  jaM[jaM[localrow]+4] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid] + 1;
	  jaM[jaM[localrow]+5] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid] + 2;
	  jaM[jaM[localrow]+6] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]);
	  jaM[jaM[localrow]+7] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]) + 1;
	  jaM[jaM[localrow]+8] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]) + 2;
	  
	  aM[jaM[localrow]] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(3,offset+index);
	  aM[jaM[localrow]+1] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(3,offset+index);
	  aM[jaM[localrow]+2] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(3,offset+index);
	  aM[jaM[localrow]+3] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(4,offset+index);
	  aM[jaM[localrow]+4] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(4,offset+index);
	  aM[jaM[localrow]+5] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(4,offset+index);
	  aM[jaM[localrow]+6] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(5,offset+index);
	  aM[jaM[localrow]+7] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(5,offset+index);
	  aM[jaM[localrow]+8] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(5,offset+index);
	}
      else if (package == PETSc)
	{
	  mM = 1;
	  //	  idxmM[0] = update[localrow];
	  idxmM[0] = (*petscIndexM)(i,j);
	  nM = 10;
	  //	  idxnM[0] = update[localrow];
	  //	  idxnM[1] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc;
	  //	  idxnM[2] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 1;
	  //	  idxnM[3] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2;
	  //	  idxnM[4] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid];
	  //	  idxnM[5] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid] + 1;
	  //	  idxnM[6] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid] + 2;
	  //	  idxnM[7] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]);
	  //	  idxnM[8] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]) + 1;
	  //	  idxnM[9] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 2*(Kimera->rdimM[donorGrid]) + 2;
	  idxnM[0] = (*petscIndexM)(i,j);

	  int iProc, lRDim;
	  
	  iProc = ownerProcess(i_loc, j_loc, donorGrid);
	  lRDim =  petscLocalRHiM[donorGrid](iProc) - petscLocalRLowM[donorGrid](iProc) + 1;
	  start = petscStartIndexM[donorGrid](iProc)
	    + lRDim*(j_loc-petscLocalSLowM[donorGrid](iProc))
	    + (i_loc-petscLocalRLowM[donorGrid](iProc));
	  
	  idxnM[1] = start;

	  iProc = ownerProcess(i_loc+1, j_loc, donorGrid);
	  lRDim =  petscLocalRHiM[donorGrid](iProc) - petscLocalRLowM[donorGrid](iProc) + 1;
	  start = petscStartIndexM[donorGrid](iProc)
	    + lRDim*(j_loc-petscLocalSLowM[donorGrid](iProc))
	    + (i_loc+1-petscLocalRLowM[donorGrid](iProc));

	  idxnM[2] = start;

	  iProc = ownerProcess(i_loc+2, j_loc, donorGrid);
	  lRDim =  petscLocalRHiM[donorGrid](iProc) - petscLocalRLowM[donorGrid](iProc) + 1;
	  start = petscStartIndexM[donorGrid](iProc)
	    + lRDim*(j_loc-petscLocalSLowM[donorGrid](iProc))
	    + (i_loc+2-petscLocalRLowM[donorGrid](iProc));

	  idxnM[3] = start;

	  iProc = ownerProcess(i_loc, j_loc+1, donorGrid);
	  lRDim =  petscLocalRHiM[donorGrid](iProc) - petscLocalRLowM[donorGrid](iProc) + 1;
	  start = petscStartIndexM[donorGrid](iProc)
	    + lRDim*(j_loc+1-petscLocalSLowM[donorGrid](iProc))
	    + (i_loc-petscLocalRLowM[donorGrid](iProc));

	  idxnM[4] = start;

	  iProc = ownerProcess(i_loc+1, j_loc+1, donorGrid);
	  lRDim =  petscLocalRHiM[donorGrid](iProc) - petscLocalRLowM[donorGrid](iProc) + 1;
	  start = petscStartIndexM[donorGrid](iProc)
	    + lRDim*(j_loc+1-petscLocalSLowM[donorGrid](iProc))
	    + (i_loc+1-petscLocalRLowM[donorGrid](iProc));

	  idxnM[5] = start;

	  iProc = ownerProcess(i_loc+2, j_loc+1, donorGrid);
	  lRDim =  petscLocalRHiM[donorGrid](iProc) - petscLocalRLowM[donorGrid](iProc) + 1;
	  start = petscStartIndexM[donorGrid](iProc)
	    + lRDim*(j_loc+1-petscLocalSLowM[donorGrid](iProc))
	    + (i_loc+2-petscLocalRLowM[donorGrid](iProc));

	  idxnM[6] = start;

	  iProc = ownerProcess(i_loc, j_loc+2, donorGrid);
	  lRDim =  petscLocalRHiM[donorGrid](iProc) - petscLocalRLowM[donorGrid](iProc) + 1;
	  start = petscStartIndexM[donorGrid](iProc)
	    + lRDim*(j_loc+2-petscLocalSLowM[donorGrid](iProc))
	    + (i_loc-petscLocalRLowM[donorGrid](iProc));

	  idxnM[7] = start;

	  iProc = ownerProcess(i_loc+1, j_loc+2, donorGrid);
	  lRDim =  petscLocalRHiM[donorGrid](iProc) - petscLocalRLowM[donorGrid](iProc) + 1;
	  start = petscStartIndexM[donorGrid](iProc)
	    + lRDim*(j_loc+2-petscLocalSLowM[donorGrid](iProc))
	    + (i_loc+1-petscLocalRLowM[donorGrid](iProc));

	  idxnM[8] = start;

	  iProc = ownerProcess(i_loc+2, j_loc+2, donorGrid);
	  lRDim =  petscLocalRHiM[donorGrid](iProc) - petscLocalRLowM[donorGrid](iProc) + 1;
	  start = petscStartIndexM[donorGrid](iProc)
	    + lRDim*(j_loc+2-petscLocalSLowM[donorGrid](iProc))
	    + (i_loc+2-petscLocalRLowM[donorGrid](iProc));

	  idxnM[9] = start;

	  vM[0] = -1.;
	  vM[1] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(3,offset+index);
	  vM[2] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(3,offset+index);
	  vM[3] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(3,offset+index);
	  vM[4] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(4,offset+index);
	  vM[5] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(4,offset+index);
	  vM[6] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(4,offset+index);
	  vM[7] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(5,offset+index);
	  vM[8] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(5,offset+index);
	  vM[9] = CoeffReceiveBufferM(2,offset+index) * CoeffReceiveBufferM(5,offset+index);
	}
  }
  else if (iType == linear)
    {
      int index = Kimera->getInterpIndex(i,j,grid);
      //
      // The following changes are for
      // the PUM 
      //
      if (index < 0)
	{
	  jaM[localrow+1] = jaM[localrow];
	  aM[localrow] = 1.;
	  NrOfNonZerosM += 1;
	  return;
	}
      int i_loc = Kimera->My_i_LM[grid](index);
      int j_loc = Kimera->My_j_LM[grid](index);


      jaM[localrow+1] = jaM[localrow] + 4;
      aM[localrow] = -1.;
      
      NrOfNonZerosM += 5;
      
      
      int start = 0;
      int offset = 0;
      
      Index all;
      
      for (int k=0; k<donorGrid; k++)
	{
	  start += Kimera->rdimM[k] * Kimera->sdimM[k];
	}
      for (int g=0; g<grid; g++)
	offset += Kimera->myIntPointsM[g].getLength(0);
      
      //  cout << " offset = " << offset << " index = " << index << endl;
      
      //      if (start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid] + 1 >= 842)
      //	{
      //	  cout << "xvdf" << i_loc << " " << j_loc << "\n";
      //	}
      jaM[jaM[localrow]] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc;
      jaM[jaM[localrow]+1] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + 1;
      jaM[jaM[localrow]+2] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid];
      jaM[jaM[localrow]+3] = start + j_loc*(Kimera->rdimM[donorGrid]) + i_loc + Kimera->rdimM[donorGrid] + 1;
      
      aM[jaM[localrow]] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(2,offset+index);
      aM[jaM[localrow]+1] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(2,offset+index);
      aM[jaM[localrow]+2] = CoeffReceiveBufferM(0,offset+index) * CoeffReceiveBufferM(3,offset+index);
      aM[jaM[localrow]+3] = CoeffReceiveBufferM(1,offset+index) * CoeffReceiveBufferM(3,offset+index);
    }
}



OGEquation::OGEquation(CompositeGrid *Kimera, sparseSolver spSolver, int itm, int prec, int over) : MatrixM()
{
  OGridM = Kimera;

  packageM = spSolver;

  reInitializeM = false;
  setFromOptionsM = true;

  if (packageM == Aztec)
    {
      AZ_defaults(optionsM, paramsM);

      optionsM[AZ_solver] = AZ_bicgstab;
      optionsM[AZ_precond] = AZ_dom_decomp;
      optionsM[AZ_overlap] = 0;
      
      optionsM[AZ_subdomain_solve] = AZ_ilu;
      
      optionsM[AZ_output] = AZ_none;
      optionsM[AZ_max_iter] = 500;
      
      optionsM[AZ_conv] = AZ_noscaled;
      
      paramsM[AZ_tol] = ConvConstM*Kimera->hSquare();
      
      
      AZ_set_proc_config(proc_configM, MPI_COMM_WORLD);
    }
      
  int grid;

  N_updateM = 0;
  for (grid = 0; grid < Kimera->nrGrids(); grid++)
    {
      N_updateM += ( Kimera->ub_1(grid) - Kimera->lb_1(grid) + 1 ) * ( Kimera->ub_0(grid) - Kimera->lb_0(grid) + 1 );
    }
  
  if (packageM == Aztec)
    {
      updateM = new int[N_updateM + 1]; // +1 in case we are solving a pure Neumann problem
      // and need to add an extra equation
    }
  else if (packageM == PETSc)
    {
      updateM = new int[N_updateM + 1]; // +1 in case we are solving a pure Neumann problem
      // and need to add an extra equation
      update_indexM = new int[1];
    }
  
  int globalCount = 0;
  int index = 0;
  
  for (grid=0; grid<Kimera->nrGrids(); grid++)
  {
    for (int j=Kimera->lb_1(grid); j<=Kimera->ub_1(grid); j++)
    {
      for (int i=Kimera->lb_0(grid); i<=Kimera->ub_0(grid); i++)
      {
	updateM[index] = globalCount + i + j*Kimera->flagValuesM[grid].getLength(0);
	index++;
      }
    }
    globalCount += Kimera->flagValuesM[grid].getLength(0) * Kimera->flagValuesM[grid].getLength(1);
  }
  
  int ierr;
  if (packageM == PETSc)
    {
      //----------------------------------------
      // Create a linear solver structure
      // for use by PETSc
      //----------------------------------------
      ierr = KSPCreate(PETSC_COMM_WORLD, &kspM);
      if (ierr != 0) cout << "ERROR in kspCreate(...)..." 
			  << endl << endl;
    }
}

OGEquation::~OGEquation()
{
  delete []updateM;

  if (packageM == Aztec)
    {
      delete []lhsM;
      delete []rhsM;
      free(data_orgM);
      free(externalM);
      free(update_indexM);
      free(extern_indexM);
      
      AZ_matrix_destroy(&AmatM);
      AZ_precond_destroy(&precM);
    }
  else if (packageM == PETSc)
    {
      //      VecDestroy(x);
      //      VecDestroy(b);
      delete []update_indexM;
    }
}

void
OGEquation::setOperator(MatrixOperator operatorType, bool allNeumann, gridFunction *visc, bool reset, gridFunction *convFieldU, gridFunction *convFieldV)
{
  int me = Communication_Manager::My_Process_Number;
  int nrProcs = Communication_Manager::numberOfProcessors();

  equationM = operatorType;
  
  if (allNeumann && me == nrProcs-1)     N_updateM++;

  if (operatorType == Grad_Visc_Div || operatorType == Grad_Visc_Div_with_TD_diagonal)
    MatrixM.setupMatrix(OGridM, operatorType, N_updateM, updateM, unKnownM->flagValuesM, allNeumann, unKnownM->iTypeM, visc, reset, packageM, convFieldU, convFieldV);
  else
    MatrixM.setupMatrix(OGridM, operatorType, N_updateM, updateM, unKnownM->flagValuesM, allNeumann, unKnownM->iTypeM, visc, reset, packageM, convFieldU, convFieldV);


  if (!reset)
    {
      if (packageM == PETSc)
	{
	  VecCreateMPI(PETSC_COMM_WORLD, N_updateM, PETSC_DETERMINE, &bM);
	  VecDuplicate(bM ,&xM);
	  PetscScalar zero = 0.;
	  VecSet(xM, zero);
	}
    }
  if (reset)
    {
      if (packageM == Aztec)
	{
	  AZ_matrix_destroy(&AmatM);
	  AZ_precond_destroy(&precM);	  
	}
      else if (packageM == PETSc)
	{
	}
      
    }

  
  if (packageM == Aztec)
    {
      AZ_transform(proc_configM, &externalM, MatrixM.bindx(), MatrixM.val(), updateM, &update_indexM,
		   &extern_indexM, &data_orgM, N_updateM, 0, 0, 0, 0, AZ_MSR_MATRIX);
      if (!reset)
	{
	  lhsM = new double[N_updateM + data_orgM[AZ_N_external]];
	  for (int row = 0; row < N_updateM + data_orgM[AZ_N_external]; row++)
	    lhsM[row] = 0.;
	}
	
    }
  else if (packageM == PETSc)
    {
      int ierr;
      KSPSetOperators(kspM, MatrixM.PmatM, MatrixM.PmatM);
      //      SLESGetKSP(slesM, &kspM);
      KSPGetPC(kspM, &pcM);
      double atol = ConvConstM*unKnownM->myGridM->hSquare();
      KSPSetPCSide(kspM, PC_RIGHT);
      KSPSetTolerances(kspM, rtolM, atol, PETSC_DEFAULT, PETSC_DEFAULT);
      ierr = KSPSetType(kspM, KSPBCGS);CHKERRCONTINUE(ierr);
      PCSetType(pcM, PCBJACOBI);
      //      KSPSetRhs(kspM, bM);
      //      KSPSetSolution(kspM, xM);
      KSPSetUp(kspM);  
      //SLESSetUp(slesM, xM, bM);

      PCBJacobiGetSubKSP(pcM, &nlocalM, &firstM, &subkspM);

      PC subpc;
      for (int i = 0; i < nlocalM; i++)
	{
          int ierr; 
          KSPGetPC(subkspM[i], &subpc);
          //SLESGetKSP(subslesM[i], &subksp);
	  PCSetType(subpc, PCILU);
	  ierr = KSPSetType(*subkspM, KSPPREONLY);CHKERRCONTINUE(ierr);
	}
    }
  
  if (packageM == Aztec)
    {
      //----------------------------------------
      // Use the Aztec2.x interface when solving.
      // We need ADT's for the matrix and its
      // preconditioner.
      //----------------------------------------
      AmatM = AZ_matrix_create(N_updateM);
      AZ_set_MSR(AmatM, MatrixM.bindx(), MatrixM.val(), data_orgM, 0, NULL, AZ_LOCAL);
      
      //----------------------------------------
      // If we change AZ_precondition, we can
      // use our own preconditioning routine
      //----------------------------------------
      precM = AZ_precond_create(AmatM, AZ_precondition, NULL);
    }
}


void
OGEquation::infNorm()
{
  if (packageM == PETSc)
    {
      double NORM;
      MatNorm(MatrixM.PmatM, NORM_INFINITY, &NORM);
      if (Communication_Manager::My_Process_Number == 0)
	cout << " Infinity_Norm of matrix is " << NORM << endl;
    }
  else
    {
      if (Communication_Manager::My_Process_Number == 0)
	cout << " Sorry Norm is not available for Aztec matrices \n";
    }
}

void
OGEquation::destroyOperator()
{
  if (packageM == Aztec)
    {
      free(data_orgM);
      free(externalM);
      free(update_indexM);
      free(extern_indexM);
    }
  else if (packageM == PETSc)
    {
      MatZeroEntries(MatrixM.PmatM);
    }
}

double
DMSRMatrix::xxCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *xrr, doubleSerialArray *xss, doubleSerialArray *yrr, doubleSerialArray *yss, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3])
{
  return ( pow((*ys)(i,j),2.)*(coeff[1][2] - 2*coeff[1][1] + coeff[1][0])/pow(r_step,2.) - 2*(*yr)(i,j)*(*ys)(i,j)*((coeff[2][2] - coeff[0][2]) - (coeff[2][0] - coeff[0][0]))/(4.*r_step*s_step) + pow((*yr)(i,j),2.)*(coeff[2][1] - 2.*coeff[1][1] + coeff[0][1])/pow(s_step,2.)) / pow((*sqrtOfG)(i,j),2.)
    + ( (pow((*ys)(i,j),2.)*(*yrr)(i,j) - 2.*(*yr)(i,j)*(*ys)(i,j)*((*yr)(i,j+1)-(*yr)(i,j-1))/(2.*s_step) + pow((*yr)(i,j),2.)*(*yss)(i,j))*((*xs)(i,j)*(coeff[1][2] - coeff[1][0])/(2.*r_step) - (*xr)(i,j)*(coeff[2][1] - coeff[0][1])/(2.*s_step)) + (pow((*ys)(i,j),2.)*(*xrr)(i,j) - 2.*(*yr)(i,j)*(*ys)(i,j)*((*xr)(i,j+1)-(*xr)(i,j-1))/(2.*s_step) + pow((*yr)(i,j),2.)*(*xss)(i,j))*((*yr)(i,j)*(coeff[2][1] - coeff[0][1])/(2.*s_step) - (*ys)(i,j)*(coeff[1][2] - coeff[1][0])/(2.*r_step)) ) / pow((*sqrtOfG)(i,j),3.);
}

double
DMSRMatrix::yyCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *xrr, doubleSerialArray *xss, doubleSerialArray *yrr, doubleSerialArray *yss, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3])
{
  return ( pow((*xs)(i,j),2.)*(coeff[1][2] - 2*coeff[1][1] + coeff[1][0])/pow(r_step,2.) - 2*(*xr)(i,j)*(*xs)(i,j)*((coeff[2][2] - coeff[0][2]) - (coeff[2][0] - coeff[0][0]))/(4.*r_step*s_step) + pow((*xr)(i,j),2.)*(coeff[2][1] - 2*coeff[1][1] + coeff[0][1])/pow(s_step,2.)) / pow((*sqrtOfG)(i,j),2.)
      + ( (pow((*xs)(i,j),2.)*(*yrr)(i,j) - 2.*(*xr)(i,j)*(*xs)(i,j)*((*yr)(i,j+1)-(*yr)(i,j-1))/(2.*s_step) + pow((*xr)(i,j),2.)*(*yss)(i,j))*((*xs)(i,j)*(coeff[1][2] - coeff[1][0])/(2.*r_step) - (*xr)(i,j)*(coeff[2][1] - coeff[0][1])/(2.*s_step)) + (pow((*xs)(i,j),2.)*(*xrr)(i,j) - 2.*(*xr)(i,j)*(*xs)(i,j)*((*xr)(i,j+1)-(*xr)(i,j-1))/(2.*s_step) + pow((*xr)(i,j),2.)*(*xss)(i,j))*((*yr)(i,j)*(coeff[2][1] - coeff[0][1])/(2.*s_step) - (*ys)(i,j)*(coeff[1][2] - coeff[1][0])/(2.*r_step)) ) / pow((*sqrtOfG)(i,j),3.);
}

double
DMSRMatrix::xCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3])
{
  return ( ( (*ys)(i,j)*(coeff[1][2] - coeff[1][0]) / (2*r_step) - (*yr)(i,j)*(coeff[2][1] - coeff[0][1]) / (2*s_step)) / (*sqrtOfG)(i,j) );
}

double
DMSRMatrix::upWindxCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3], doubleSerialArray *cF)
{
  if ((*cF)(i,j) < 0)
    {
      return ( ( (*ys)(i,j)*(coeff[1][1] - coeff[1][0]) / (r_step) - (*yr)(i,j)*(coeff[2][1] - coeff[0][1]) / (2*s_step)) / (*sqrtOfG)(i,j) );
    }
  else
    {
      return ( ( (*ys)(i,j)*(coeff[1][2] - coeff[1][1]) / (r_step) - (*yr)(i,j)*(coeff[2][1] - coeff[0][1]) / (2*s_step)) / (*sqrtOfG)(i,j) );
    }
}

double
DMSRMatrix::yCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3])
{
  return ( -(*xs)(i,j)*(coeff[1][2] - coeff[1][0]) / (2*r_step) + (*xr)(i,j)*(coeff[2][1] - coeff[0][1]) / (2*s_step)) / (*sqrtOfG)(i,j);
}

double
DMSRMatrix::upWindyCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3], doubleSerialArray *cF)
{
  if ((*cF)(i,j) < 0)
    {
      return ( -(*xs)(i,j)*(coeff[1][2] - coeff[1][0]) / (2*r_step) + (*xr)(i,j)*(coeff[1][1] - coeff[0][1]) / (s_step)) / (*sqrtOfG)(i,j);
    }
  else
    {
      return ( -(*xs)(i,j)*(coeff[1][2] - coeff[1][0]) / (2*r_step) + (*xr)(i,j)*(coeff[2][1] - coeff[1][1]) / (s_step)) / (*sqrtOfG)(i,j);
    }
}

double
DMSRMatrix::NeumannCoeff(int i, int j, int side, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3])
{
  if (side==0 || side==1) // low r or high r
  {
    return ( (pow((*xs)(i,j),2.) + pow((*ys)(i,j),2.))*(coeff[1][2]-coeff[1][0])/(2.*r_step) - ((*xr)(i,j)*(*xs)(i,j) + (*yr)(i,j)*(*ys)(i,j))*(coeff[2][1]-coeff[0][1])/(2.*s_step) ) / ( (*sqrtOfG)(i,j)*sqrt(pow((*xs)(i,j),2.)+pow((*ys)(i,j),2.)) );
  }
  else // low s or high s
  {
    return ( -((*xr)(i,j)*(*xs)(i,j) + (*yr)(i,j)*(*ys)(i,j))*(coeff[1][2]-coeff[1][0])/(2.*r_step) + (pow((*xr)(i,j),2.) + pow((*yr)(i,j),2.))*(coeff[2][1]-coeff[0][1])/(2.*s_step) ) / ( (*sqrtOfG)(i,j)*sqrt(pow((*xr)(i,j),2.)+pow((*yr)(i,j),2.)) );
  }
}

void
DMSRMatrix::getInterpolationCoefficients(CompositeGrid *Kimera, interpolationType iType)
{
  int k, p, grid, proc, nrProcs, myRank, currGrid;
  Index all;

  /*
    First check if we really need to interpolate between component grids,
    i.e. if there is more than one component grid
  */
  if (Kimera->nmbrOfGridsM == 1)
    return;

  MPI_Comm myComm;
  MPI_Group myGroup;
  
  MPI_Comm_group(MPI_COMM_WORLD, &myGroup);
  MPI_Comm_create(MPI_COMM_WORLD, myGroup, &myComm);
  MPI_Comm_size(myComm, &nrProcs);
  MPI_Comm_rank(myComm, &myRank);
 
  
  intSerialArray offsets(nrProcs,Kimera->nmbrOfGridsM,Kimera->nmbrOfGridsM);
  
  MPI_Datatype *receiveType = new MPI_Datatype[nrProcs];
  MPI_Datatype *sendType = new MPI_Datatype[nrProcs];


  double *theSendArray, *theReceiveArray;
  theSendArray = CoeffSendBufferM.getDataPointer();
  theReceiveArray = CoeffReceiveBufferM.getDataPointer();

  //  cout << "Pointers are = " << theSendArray << "," << theReceiveArray << endl;

  MPI_Request *Request1 = new MPI_Request[nrProcs];
  MPI_Request *Request2 = new MPI_Request[nrProcs];

  /*
    First post all receives needed
  */

  int recHandle;
  int currIndex=0;
  int *block_lengths = new int[(Kimera->nmbrOfGridsM)*(Kimera->nmbrOfGridsM)];
  int *displacements = new int[(Kimera->nmbrOfGridsM)*(Kimera->nmbrOfGridsM)];

  int iW;
  if (iType == quadratic)
    iW = 3;
  else if (iType == cubic)
    iW = 4;
  else // (iType == linear)
    iW = 2;

  for (p=0; p<nrProcs; p++)
  {
    for (k=0; k<Kimera->nmbrOfGridsM; k++)
    {
      displacements[k] = 0;
      for (grid=0; grid<k; grid++)
	displacements[k] += Kimera->myIntPointsM[grid].getLength(0) * 2*iW;
      
      for (proc=0; proc<p; proc++)
	displacements[k] += Kimera->myNrOfIntPointsM(k,proc) * 2*iW;
      
      block_lengths[k] = Kimera->myNrOfIntPointsM(k,p) * 2*iW;
    }
    
    MPI_Type_indexed(Kimera->nmbrOfGridsM, block_lengths, displacements, MPI_DOUBLE, &(receiveType[p]));
    MPI_Type_commit(&(receiveType[p]));
    
    recHandle = myRank*10+4;

    double *rpointer = 0x0;
    if (theReceiveArray == 0) // receiving message of len=0
      MPI_Irecv(rpointer, 1, receiveType[p], p, recHandle, myComm, &(Request2[p]));
    else
      MPI_Irecv(theReceiveArray, 1, receiveType[p], p, recHandle, myComm, &(Request2[p]));
  }
  
  currIndex = 0;
  offsets = 0;
  for (k=0; k<Kimera->nmbrOfGridsM; k++)
  {
    for (p=0; p<nrProcs; p++)
    {
      for (currGrid=0; currGrid<Kimera->nmbrOfGridsM; currGrid++)
      {
	for (grid=0; grid<k; grid++)
	  offsets(p,k,currGrid) += Kimera->intInterpolationLocationM[grid].getLength(0) * 2*iW;
	
	for (proc=0; proc<p; proc++)
	  for (grid=0; grid<Kimera->nmbrOfGridsM; grid++)
	    offsets(p,k,currGrid) += Kimera->nrOfPointsToComputeM(k,grid,proc) * 2*iW;
	
	for (grid=0; grid<currGrid; grid++)
	  offsets(p,k,currGrid) += Kimera->nrOfPointsToComputeM(k,grid,p) * 2*iW;
	
      }
      
    }
    Index sendIndex(currIndex, (Kimera->intInterpolationLocationM[k]).getLength(0));
    CoeffSendBufferM(all,sendIndex) = computeInterpolationWeights(k, Kimera, iType);
    currIndex += (Kimera->intInterpolationLocationM[k]).getLength(0);
  }
  
  /*
    Send all data to correct process
  */
  int sendHandle;
  for (p=0; p<nrProcs; p++)
  {
    int index = 0;
    for (int i=0; i<Kimera->nmbrOfGridsM; i++)
      for (int j=0; j<Kimera->nmbrOfGridsM; j++)
      {
	block_lengths[index] = Kimera->nrOfPointsToComputeM(j,i,p) * 2*iW;
	displacements[index] = offsets(p,j,i);
	index++;
      }

    MPI_Type_indexed((Kimera->nmbrOfGridsM)*(Kimera->nmbrOfGridsM),block_lengths,displacements,MPI_DOUBLE, &(sendType[p]));
    MPI_Type_commit(&(sendType[p]));
    
    sendHandle = p*10+4;
    double *spoint = 0x0;
    if (theSendArray == 0) // message of len=0
      MPI_Isend(spoint, 1, sendType[p], p, sendHandle, myComm, &(Request1[p]));
    else
      MPI_Isend(theSendArray, 1, sendType[p], p, sendHandle, myComm, &(Request1[p]));
  }
  
  delete[] block_lengths;
  delete[] displacements;

  MPI_Status *status1 = new MPI_Status[nrProcs];
  MPI_Status *status2 = new MPI_Status[nrProcs];
  
  MPI_Waitall(nrProcs, Request2, status2);
  MPI_Waitall(nrProcs, Request1, status1);

  for (p=0; p<nrProcs; p++)
    {
      MPI_Type_free(&(sendType[p]));
      MPI_Type_free(&(receiveType[p]));
    }
  delete[] Request1;
  delete[] Request2;
  delete[] status1;
  delete[] status2;
  delete[] sendType;
  delete[] receiveType;
  
}

doubleSerialArray
DMSRMatrix::computeInterpolationWeights(int k, CompositeGrid *Kimera, interpolationType iType)
#define i_L (Kimera->intInterpolationLocationM[k](all,0))
#define j_L (Kimera->intInterpolationLocationM[k](all,1))
{
  Index all;

  int iW;
  if (iType == quadratic)
    iW = 3;
  else if (iType == cubic)
    iW = 4;
  else // (iType == linear)
    iW = 2;

  int sumOther = i_L.getLength(0);
//  printf("   Sumother = %d\n",sumOther);

  Index sum(0,sumOther);
  
  doubleSerialArray result(2*iW,sumOther);
  
  if (sumOther > 0)
  {
    if (iType == cubic)
      {
	result(0,sum) = epsilon(-1,k,Kimera);
	result(1,sum) = epsilon(0,k,Kimera);
	result(2,sum) = epsilon(1,k,Kimera);
	result(3,sum) = epsilon(2,k,Kimera);
	result(4,sum) = zeta(-1,k,Kimera);
	result(5,sum) = zeta(0,k,Kimera);
	result(6,sum) = zeta(1,k,Kimera);
	result(7,sum) = zeta(2,k,Kimera);
      }
    else if (iType == quadratic)
      {
	result(0,sum) = alpha(-1,k,Kimera);
	result(1,sum) = alpha(0,k,Kimera);
	result(2,sum) = alpha(1,k,Kimera);
	result(3,sum) = beta(-1,k,Kimera);
	result(4,sum) = beta(0,k,Kimera);
	result(5,sum) = beta(1,k,Kimera);
      }
    else if (iType == linear)
      {
	result(0,sum) = gamma(-1,k,Kimera);
	result(1,sum) = gamma(0,k,Kimera);
	result(2,sum) = delta(-1,k,Kimera);
	result(3,sum) = delta(0,k,Kimera);
      }
    return result;
  }
  return result;
}
#undef i_L
#undef j_L

doubleSerialArray &
DMSRMatrix::alpha(int i,int k, CompositeGrid *Kimera)
#define rstep (Kimera->r_stepM[k])
#define r_L (Kimera->interpolationCoordinatesM[k](all,0))
#define r (Kimera->doubleInterpolationLocationM[k](all,0))
{
  Index all;

  int sumOther = r.getLength(0);

  if (i == -1)
  {
    doubleSerialArray & result1 = ( (r_L - (r+1-1)*rstep)*(r_L - (r+2-1)*rstep) / (2*rstep*rstep) );
    return result1.reshape(1,sumOther);
  }
  else if (i == 0)
  {
    doubleSerialArray & result2 = - (r_L - (r-1)*rstep)*(r_L - (r+2-1)*rstep) / (rstep*rstep);
    return result2.reshape(1,sumOther);
  }
  else 
  {
    doubleSerialArray & result3 = (r_L - (r-1)*rstep)*(r_L - (r+1-1)*rstep) / (2*rstep*rstep);
    return result3.reshape(1,sumOther);
  }
  
}
#undef rstep
#undef r_L
#undef r

doubleSerialArray &
DMSRMatrix::beta(int i,int k, CompositeGrid *Kimera)
#define sstep (Kimera->s_stepM[k])  
#define s_L (Kimera->interpolationCoordinatesM[k](all,1))
#define s (Kimera->doubleInterpolationLocationM[k](all,1))  
{
  Index all;

  int sumOther = s.getLength(0);

  if (i == -1) 
    {
      doubleSerialArray & result1 = (s_L - (s+1-1)*sstep)*(s_L - (s+2-1)*sstep) / (2*sstep*sstep);
      return result1.reshape(1,sumOther);
    }
  else if (i == 0)
    {
      doubleSerialArray & result2 = - (s_L - (s-1)*sstep)*(s_L - (s+2-1)*sstep) / (sstep*sstep);
      return result2.reshape(1,sumOther);
    }
  else 
    {
      doubleSerialArray & result3 = (s_L - (s-1)*sstep)*(s_L - (s+1-1)*sstep) / (2*sstep*sstep);
      return result3.reshape(1,sumOther);
    }

}
#undef sstep
#undef s_L
#undef s

doubleSerialArray &
DMSRMatrix::gamma(int i,int k, CompositeGrid *Kimera)
#define rstep (Kimera->r_stepM[k])
#define r_L (Kimera->interpolationCoordinatesM[k](all,0))
#define r (Kimera->doubleInterpolationLocationM[k](all,0))
{
  Index all;

  int sumOther = r.getLength(0);

  if (i == -1)
  {
    doubleSerialArray & result1 = (r_L - (r+1-1)*rstep) / (-rstep);
    return result1.reshape(1,sumOther);
  }
  else 
  {
    doubleSerialArray & result2 = (r_L - (r-1)*rstep) / (rstep);
    return result2.reshape(1,sumOther);
  }
}
#undef rstep
#undef r_L
#undef r

doubleSerialArray &
DMSRMatrix::delta(int i,int k, CompositeGrid *Kimera)
#define sstep (Kimera->s_stepM[k])  
#define s_L (Kimera->interpolationCoordinatesM[k](all,1))
#define s (Kimera->doubleInterpolationLocationM[k](all,1))  
{
  Index all;

  int sumOther = s.getLength(0);

  if (i == -1) 
    {
      doubleSerialArray & result1 = (s_L - (s+1-1)*sstep) / (-sstep);
      return result1.reshape(1,sumOther);
    }
  else 
    {
      doubleSerialArray & result2 = (s_L - (s-1)*sstep) / (sstep);
      return result2.reshape(1,sumOther);
    }
}
#undef sstep
#undef s_L
#undef s

doubleSerialArray &
DMSRMatrix::epsilon(int i,int k, CompositeGrid *Kimera)
#define rstep (Kimera->r_stepM[k])
#define r_L (Kimera->interpolationCoordinatesM[k](all,0))
#define r (Kimera->doubleInterpolationLocationM[k](all,0))
{
  Index all;

  int sumOther = r.getLength(0);

  if (i == -1)
    {
      doubleSerialArray & result1 = - (r_L - (r+1-1)*rstep)*(r_L - (r+2-1)*rstep)*(r_L - (r+3-1)*rstep) / (6*rstep*rstep*rstep);
      return result1.reshape(1,sumOther);
    }
  else if (i == 0)
    {
      doubleSerialArray & result2 = (r_L - (r-1)*rstep)*(r_L - (r+2-1)*rstep)*(r_L - (r+3-1)*rstep) / (2*rstep*rstep*rstep);
      return result2.reshape(1,sumOther);
    }
  else if (i == 1)
    {
      doubleSerialArray & result3 = - (r_L - (r-1)*rstep)*(r_L - (r+1-1)*rstep)*(r_L - (r+3-1)*rstep) / (2*rstep*rstep*rstep);
      return result3.reshape(1,sumOther);
    }
  else
    {
      doubleSerialArray & result4 = (r_L - (r-1)*rstep)*(r_L - (r+1-1)*rstep)*(r_L - (r+2-1)*rstep) / (6*rstep*rstep*rstep);
      return result4.reshape(1,sumOther);
    } 
}
#undef rstep
#undef r_L
#undef r

doubleSerialArray &
DMSRMatrix::zeta(int i,int k, CompositeGrid *Kimera)
#define sstep (Kimera->s_stepM[k])  
#define s_L (Kimera->interpolationCoordinatesM[k](all,1))
#define s (Kimera->doubleInterpolationLocationM[k](all,1))  
{
  Index all;

  int sumOther = s.getLength(0);

  if (i == -1) 
    {
      doubleSerialArray & result1 = (s_L - (s+1-1)*sstep)*(s_L - (s+2-1)*sstep)*(s_L - (s+3-1)*sstep) / (6*sstep*sstep*sstep);
      return result1.reshape(1,sumOther);
    }
  else if (i == 0)
    {
      doubleSerialArray & result2 = (s_L - (s-1)*sstep)*(s_L - (s+2-1)*sstep)*(s_L - (s+3-1)*sstep) / (2*sstep*sstep*sstep);
      return result2.reshape(1,sumOther);
    }
  else if (i == 1)
    {
      doubleSerialArray & result3 = - (s_L - (s-1)*sstep)*(s_L - (s+1-1)*sstep)*(s_L - (s+3-1)*sstep) / (2*sstep*sstep*sstep);
      return result3.reshape(1,sumOther);
    }
  else
    { 
      doubleSerialArray & result4 = (s_L - (s-1)*sstep)*(s_L - (s+1-1)*sstep)*(s_L - (s+2-1)*sstep) / (6*sstep*sstep*sstep);
      return result4.reshape(1,sumOther);
    }
}
#undef sstep
#undef s_L
#undef s

void
OGEquation::solve(gridFunction & RHS, CompositeGrid *Kimera, double dt, bool allNeumann, bool setBCs)
{
  int me = Communication_Manager::My_Process_Number;
  int nrProcs = Communication_Manager::numberOfProcessors();

  int row = 0;
  
  if (setBCs)
    {
      //      cout << "Setting BCs in OGEquation::solve " << endl;
      setBoundaryValues(RHS, Kimera);
    }

  if (equationM == Poisson_with_TD_diagonal || equationM == Grad_Visc_Div_with_TD_diagonal || equationM == PUM_with_TD_diagonal || equationM == Conv_varDiff_with_TD_diagonal)
    {
      //    cout << "Adding to diagonal in OGEquation::solve " << endl;
      MatrixM.addToDiagonal(Kimera, dt, unKnownM->flagValuesM, update_indexM, updateM, packageM);
    }
  PetscScalar value;
  
  if (packageM == Aztec)
    {
      for (int grid = 0; grid < OGridM->nrGrids(); grid++)
	{
	  intSerialArray *Flag = unKnownM->flagValuesM[grid].getSerialArrayPointer();
	  doubleSerialArray *field = RHS.fieldValuesM[grid].getSerialArrayPointer();
	  
	  for (int j = Kimera->lb_1(grid); j<=Kimera->ub_1(grid); j++)
	    {
	      for (int i = Kimera->lb_0(grid); i<=Kimera->ub_0(grid); i++)
		{
		  if ((*Flag)(i,j) == 0) // hole point
		    rhsM[update_indexM[row]] = 0.000001;
		  else if ((*Flag)(i,j) < 0) // interpolation point
		    rhsM[update_indexM[row]] = 0.;
		  else // some sort of discretization point
		    {
		      if ((*Flag)(i,j) == grid + 1) // inner point, rhs = forcing
			if (equationM == Poisson_with_TD_diagonal || equationM == Grad_Visc_Div_with_TD_diagonal || equationM == PUM_with_TD_diagonal || equationM == Conv_varDiff_with_TD_diagonal)
			  rhsM[update_indexM[row]] = - (*field)(i,j) / dt;
			else
			  rhsM[update_indexM[row]] = (*field)(i,j);
		      
		      else if ((*Flag)(i,j) == grid + 1000) // extrapolation, low_r
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 1 + 1000) // Dirichlet, low_r
			rhsM[update_indexM[row]] = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 2 + 1000) // Neumann, low_r
			rhsM[update_indexM[row]] = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 3 + 1000) // periodic, low_r
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 4 + 1000) // robin, low_r
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 2000) // extrapolation, hi_r
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 1 + 2000) // Dirichlet, hi_r
			rhsM[update_indexM[row]] = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 2 + 2000) // Neumann, hi_r
			rhsM[update_indexM[row]] = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 3 + 2000) // periodic, hi_r
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 4 + 2000) // robin, hi_r
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 3000) // extrapolation, low_s
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 1 + 3000) // Dirichlet, low_s
			rhsM[update_indexM[row]] = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 2 + 3000) // Neumann, low_s
			rhsM[update_indexM[row]] = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 3 + 3000) // periodic, low_s
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 4 + 3000) // robin, low_s
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 4000) // extrapolation, hi_s
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 1 + 4000) // Dirichlet, hi_s
			rhsM[update_indexM[row]] = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 2 + 4000) // Neumann, hi_s
			rhsM[update_indexM[row]] = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 3 + 4000) // periodic, hi_s
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == grid + 4 + 4000) // robin, hi_s
			rhsM[update_indexM[row]] = 0.;
		      
		      else if ((*Flag)(i,j) == 10000) // corner low_r, low_s
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == 20000) // corner low_r, hi_s
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == 30000) // corner hi_r, low_s
			rhsM[update_indexM[row]] = 0.;
		      else if ((*Flag)(i,j) == 40000) // corner hi_r, hi_s
			rhsM[update_indexM[row]] = 0.;
		      
		    }
		  //----------------------------------------
		  // Scale RHS with row scale factor
		  //----------------------------------------
		  rhsM[update_indexM[row]] *= MatrixM.rowScale(row);
		  
		  row++;
		}
	    }
	}
    }
  else if (packageM == PETSc)
    {
      int count = MatrixM.startRow;
      for (int grid = 0; grid < OGridM->nrGrids(); grid++)
	{
	  intSerialArray *Flag = unKnownM->flagValuesM[grid].getSerialArrayPointer();
	  doubleSerialArray *field = RHS.fieldValuesM[grid].getSerialArrayPointer();
	  
	  for (int j = Kimera->lb_1(grid); j <= Kimera->ub_1(grid); j++)
	    {
	      for (int i = Kimera->lb_0(grid); i <= Kimera->ub_0(grid); i++)
		{
		  if ((*Flag)(i,j) == 0) // hole point
		    value = 0.000001;
		  else if ((*Flag)(i,j) < 0) // interpolation point
		    value = 0.;
		  else // some sort of discretization point
		    {
		      if ((*Flag)(i,j) == grid + 1) // inner point, rhs = forcing
			if (equationM == Poisson_with_TD_diagonal || equationM == Grad_Visc_Div_with_TD_diagonal || equationM == PUM_with_TD_diagonal || equationM == Conv_varDiff_with_TD_diagonal)
			  value = - (*field)(i,j) / dt;
			else
			  value = (*field)(i,j);
		      
		      else if ((*Flag)(i,j) == grid + 1000) // extrapolation, low_r
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 1 + 1000) // Dirichlet, low_r
			value = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 2 + 1000) // Neumann, low_r
			value = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 3 + 1000) // periodic, low_r
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 4 + 1000) // robin, low_r
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 2000) // extrapolation, hi_r
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 1 + 2000) // Dirichlet, hi_r
			value = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 2 + 2000) // Neumann, hi_r
			value = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 3 + 2000) // periodic, hi_r
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 4 + 2000) // robin, hi_r
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 3000) // extrapolation, low_s
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 1 + 3000) // Dirichlet, low_s
			value = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 2 + 3000) // Neumann, low_s
			value = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 3 + 3000) // periodic, low_s
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 4 + 3000) // robin, low_s
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 4000) // extrapolation, hi_s
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 1 + 4000) // Dirichlet, hi_s
			value = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 2 + 4000) // Neumann, hi_s
			value = (*field)(i,j);
		      else if ((*Flag)(i,j) == grid + 3 + 4000) // periodic, hi_s
			value = 0.;
		      else if ((*Flag)(i,j) == grid + 4 + 4000) // robin, hi_s
			value = 0.;
			  
		      else if ((*Flag)(i,j) == 10000) // corner low_r, low_s
			value = 0.;
		      else if ((*Flag)(i,j) == 20000) // corner low_r, hi_s
			value = 0.;
		      else if ((*Flag)(i,j) == 30000) // corner hi_r, low_s
			value = 0.;
		      else if ((*Flag)(i,j) == 40000) // corner hi_r, hi_s
			value = 0.;
			  
		    }
		  //----------------------------------------
		  // Scale RHS with row scale factor
		  //----------------------------------------
		  value *= MatrixM.rowScale(row);
		  
		  VecSetValues(bM, 1, &(count), &value, INSERT_VALUES);
		  row++;
		  count++;
		}
	    }
	}
    }
  if (allNeumann)
    {
      if (packageM == Aztec)
	{
	  if (me == nrProcs - 1)
	    {
	      lhsM[update_indexM[row]] = 0.;
	      rhsM[update_indexM[row]] = 0.;
	    }
	}
      else if (packageM == PETSc)
	{
	  if (me == nrProcs - 1)
	    {
	      value = 0.;
	      VecSetValues(bM, 1, &(updateM[row]), &value, INSERT_VALUES);
	    }
	}
    }
  
  
  if (packageM == Aztec)
    {
      if (optionsM[AZ_precond] == AZ_dom_decomp) optionsM[AZ_keep_info] = 1;
      
      //
      // solve the equation system
      //
      AZ_iterate(lhsM, rhsM, optionsM, paramsM, statusM, proc_configM, AmatM, precM, NULL);
      
      if (optionsM[AZ_precond] == AZ_dom_decomp) optionsM[AZ_pre_calc] = AZ_reuse;
      
      if (statusM[AZ_why] != AZ_normal && me == 0)
	cout << "WARNING!, solution is not sufficiently converged after " 
	     << statusM[AZ_its] 
	     << " iterations" << endl;
      
//       if (allNeumann)
// 	{
// 	  if (me == nrProcs - 1)
// 	    {
// 	      cout << "  Compatability condition satisfied to " << lhs[update_index[row]] << endl;
// 	    }
// 	}
  
      row = 0;
      
      for (int k = 0; k < OGridM->nrGrids(); k++)
	{
	  doubleSerialArray *field = unKnownM->fieldValuesM[k].getSerialArrayPointer();
	  
	  for (int j = Kimera->lb_1(k); j <= Kimera->ub_1(k); j++)
	    {
	      for (int i = Kimera->lb_0(k); i <= Kimera->ub_0(k); i++)
		{
		  (*field)(i,j) = lhsM[update_indexM[row]];
		  row++;
		}
	    }
	}
    }

  else if (packageM == PETSc)
    {
      int its;

      if (setFromOptionsM)
	{
          KSPSetFromOptions(kspM);
	  setFromOptionsM = false;
	}

      VecAssemblyBegin(bM);
      VecAssemblyEnd(bM);
      
      if (!reInitializeM)  
	{
	  PetscBool flg = PETSC_TRUE;
	  KSPSetInitialGuessNonzero(kspM, flg);
	}

      //      ier = SLESSolve(slesM, bM, xM, &its);
      //      KSPSetRhs(kspM, bM);
      //      KSPSetSolution(kspM,xM);
      KSPSolve(kspM, bM, xM);
      KSPGetIterationNumber(kspM, &its);
      KSPConvergedReason reason;
      KSPGetConvergedReason(kspM, &reason);
      if (reason < 0)
      	{
	  if (me == 0)
	    {
	      cout << "WARNING!, solution is not sufficiently converged after " 
		   << its
		   << " iterations" << endl
		   << ", reason is = " << reason << endl;
	    }
	}

      PetscScalar *localPointer;
      VecGetArray(xM, &localPointer);

      row = 0;
      
      for (int k = 0; k < OGridM->nrGrids(); k++)
	{
	  doubleSerialArray *field = unKnownM->fieldValuesM[k].getSerialArrayPointer();
	  
	  for (int j = Kimera->lb_1(k); j <= Kimera->ub_1(k); j++)
	    {
	      for (int i = Kimera->lb_0(k); i <= Kimera->ub_0(k); i++)
		{
		  (*field)(i,j) = localPointer[row];
		  row++;
		}
	    }
	}
      VecRestoreArray(xM, &localPointer);
    } 
  unKnownM->updateGhostBoundaries();

  if (equationM == Poisson_with_TD_diagonal || equationM == Grad_Visc_Div_with_TD_diagonal || equationM == PUM_with_TD_diagonal || equationM == Conv_varDiff_with_TD_diagonal)
    {
      MatrixM.addToDiagonal(Kimera, -dt, unKnownM->flagValuesM, update_indexM, updateM, packageM);
    }
  
}

void
OGEquation::setBoundaryValues(gridFunction & RHS, CompositeGrid *Kimera)
{
#define ys(x,y) (Kimera->ysM[k](x,y))
#define yr(x,y) (Kimera->yrM[k](x,y))
#define xs(x,y) (Kimera->xsM[k](x,y))
#define xr(x,y) (Kimera->xrM[K](x,y))
#define X(a,b) (Kimera->xM[k](a,b))
#define Y(a,b) (Kimera->yM[k](a,b))
    
  int k;

  for (k = 0; k < Kimera->nrGrids(); k++)
  {
    Index ixb(0, Kimera->rdimM[k]);
    int low_r(1), hi_r((Kimera->rdimM[k]) - 2);
    Index iyb(0, Kimera->sdimM[k]);
    int low_s(1), hi_s((Kimera->sdimM[k]) - 2); 

    for (int side = Side(lowR); side <= Side(hiS); side++)
      {
	if (unKnownM->getBoundaryType(k, (Side) side) == INTERPOLATION) // interpolation
	  {
	  }
	else if (unKnownM->getBoundaryType(k, (Side) side) == DIRICHLET) // Dirichlet
	  {
	  }
	else if (unKnownM->getBoundaryType(k, (Side) side) == NEUMANN) // Neumann
	  {
	    if (side == 0)
	      {
		RHS(k)(low_r-1,iyb) = 
		  Kimera->normalVector_R_Component(k,side,ixb,low_r, hi_r,iyb,low_s,hi_s)*unKnownM->bc_x(X(low_r,iyb),Y(low_r,iyb),unKnownM->timeM, k, side) + Kimera->normalVector_S_Component(k,side,ixb,low_r, hi_r,iyb,low_s,hi_s)*unKnownM->bc_y(X(low_r,iyb),Y(low_r,iyb),unKnownM->timeM, k, side);
	      }
	    else if (side==1)
	      {
		RHS(k)(hi_r+1,iyb) = 
		  Kimera->normalVector_R_Component(k,side,ixb,low_r, hi_r,iyb,low_s,hi_s)*unKnownM->bc_x(X(hi_r,iyb),Y(hi_r,iyb),unKnownM->timeM, k, side) + Kimera->normalVector_S_Component(k,side,ixb,low_r, hi_r,iyb,low_s,hi_s)*unKnownM->bc_y(X(hi_r,iyb),Y(hi_r,iyb),unKnownM->timeM, k, side);
	      }
	    else if (side==2)
	      {
		RHS(k)(ixb,low_s-1) = 
		  Kimera->normalVector_R_Component(k,side,ixb,low_r, hi_r,iyb,low_s,hi_s)*unKnownM->bc_x(X(ixb,low_s),Y(ixb,low_s),unKnownM->timeM, k, side) + Kimera->normalVector_S_Component(k,side,ixb,low_r, hi_r,iyb,low_s,hi_s)*unKnownM->bc_y(X(ixb,low_s),Y(ixb,low_s),unKnownM->timeM, k, side);
	      }
	    else if (side==3)
	      {
		RHS(k)(ixb,hi_s+1) = 
		  Kimera->normalVector_R_Component(k,side,ixb,low_r, hi_r,iyb,low_s,hi_s)*unKnownM->bc_x(X(ixb,hi_s),Y(ixb,hi_s),unKnownM->timeM, k, side) + Kimera->normalVector_S_Component(k,side,ixb,low_r, hi_r,iyb,low_s,hi_s)*unKnownM->bc_y(X(ixb,hi_s),Y(ixb,hi_s),unKnownM->timeM, k, side);
	      }
	  }
	else if (unKnownM->getBoundaryType(k, (Side) side) == PERIODIC) // Periodic
	  {
	  }
      }

  }
  for (k = 0; k < Kimera->nrGrids(); k++)
  {
    for (int side = Side(lowR); side <= Side(hiS); side++)
      {
	if (unKnownM->getBoundaryType(k, (Side) side) == INTERPOLATION) // interpolation
	  {
	  }
	else if (unKnownM->getBoundaryType(k, (Side) side) == DIRICHLET) // Dirichlet
	  {
	    //	    cout  << "Setting p-dir, (grid,side) = (" << k << "," << side << ")" << endl;
            //            if ((Side) side == Side(lowR))
            unKnownM->dirichlet(side, k, RHS);
	  }
	else if (unKnownM->getBoundaryType(k, (Side) side) == NEUMANN) // Neumann
	  {
	  }
	else if (unKnownM->getBoundaryType(k, (Side) side) == PERIODIC) // Periodic
	  {
	  }
      }

  }
#undef ys
#undef yr
#undef xs
#undef xr
#undef X
#undef Y
}

void
DMSRMatrix::addToDiagonal(CompositeGrid *Kimera, const double dt, intArray flagValuesM[], int update_index[], int update[], sparseSolver package)
{
  int localrow = 0;
  
  if (package == Aztec)
    {
      for (int grid=0; grid<Kimera->nrGrids(); grid++)
	{
	  intSerialArray *Flag = flagValuesM[grid].getSerialArrayPointer();
	  
	  for (int j=Kimera->lb_1(grid); j<=Kimera->ub_1(grid); j++)
	    {
	      for (int i=Kimera->lb_0(grid); i<=Kimera->ub_0(grid); i++)
		{
		  if ((*Flag)(i,j) == grid + 1) // discretization point
		    {
		      aM[update_index[localrow]] -= 1./dt * rowScalingM[localrow];
		    }
		  localrow++;
		}
	    }
	}
    }
  else if (package == PETSc)
    {
      int count = startRow;
      mM = 1;
      nM = 1;
      for (int grid=0; grid<Kimera->nrGrids(); grid++)
	{
	  intSerialArray *Flag = flagValuesM[grid].getSerialArrayPointer();
	  
	  for (int j=Kimera->lb_1(grid); j<=Kimera->ub_1(grid); j++)
	    {
	      for (int i=Kimera->lb_0(grid); i<=Kimera->ub_0(grid); i++)
		{
		  if ((*Flag)(i,j) == grid + 1) // discretization point
		    {
		      idxmM[0] = count;
		      idxnM[0] = count;
		      vM[0] = -1./dt * rowScalingM[localrow];
		      MatSetValues(PmatM , mM, idxmM, nM, idxnM, vM, ADD_VALUES);
		    }
		  localrow++;
		  count++;
		}
	    }
	}
      MatAssemblyBegin(PmatM, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(PmatM, MAT_FINAL_ASSEMBLY);
    }
}

void
OGEquation::writeMatrixOperatorToFile(const char *name)
{
  int me = Communication_Manager::My_Process_Number;

  char filename1[100], filename2[100], line[100];

  sprintf(filename1,"matlabData/temp%d.m",me);

  if (packageM == Aztec)
    {
      FILE ofp;
      ofp = *stdout;
      FILE *fp = fopen(filename1, "w");
      *stdout = *fp;
      
      AZ_print_out(update_indexM, extern_indexM, updateM, externalM, MatrixM.val(), 
		   0, MatrixM.bindx(), 0, 0, 0, proc_configM, AZ_global_mat, AZ_MSR_MATRIX, 
		   N_updateM, data_orgM[AZ_N_external], 1);
      
      fflush(stdout);
      fclose(fp);
      *stdout = ofp;
      
      ifstream tempfile;
      tempfile.open(filename1,ios::in);
      
      sprintf(filename2,"matlabData/%s_matrix%d.m",name,me);
      
      fp = fopen(filename2, "w");
      
      tempfile.getline(line,100); // throw away first line 
      if (me == 0) fprintf(fp,"a = sparse(%d,%d);\n",MatrixM.getGlobalNrOfUnknowns(),MatrixM.getGlobalNrOfUnknowns());
      
      while (tempfile.getline(line,100)) fprintf(fp,"%s\n",line);

      fclose(fp);
      tempfile.close();

      remove(filename1);
    }
  else if (packageM == PETSc)
    {
      sprintf(filename2,"matlabData/%s_matrix%d.m",name,me);
      PetscViewer  v2;
      PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename2, &v2);
      PetscViewerSetFormat(v2,PETSC_VIEWER_ASCII_MATLAB);
      MatView(MatrixM.PmatM, v2);
      PetscViewerDestroy(&v2);
    }
}

void
OGEquation::setOption(OPTION option, int value)
{
  if (packageM == Aztec)
    {
      if (option == ksp_method)
	{
	  switch (value) {
	  case Cg:
	    optionsM[AZ_solver] = AZ_cg;
	    break;
	  case Gmres:
	    optionsM[AZ_solver] = AZ_gmres;
	    break;
	  case Cgs:
	    optionsM[AZ_solver] = AZ_cgs;
	    break;
	  case Tfqmr:
	    optionsM[AZ_solver] = AZ_tfqmr;
	    break;
	  case Bcgs:
	    optionsM[AZ_solver] = AZ_bicgstab;
	    break;
	  case Direct:
	    optionsM[AZ_solver] = AZ_lu;
	    break;
	  default:
	    if (Communication_Manager::My_Process_Number == 0)
	      cout << "setOption(ksp_method,value):  This value is not recognized for Aztec package " 
		   << endl;
	  }
	}
      else if (option == precond_method)
	switch (value) {
	case Jacobi:
	  optionsM[AZ_precond] = AZ_Jacobi;
	  break;
	case Symgs:
	  optionsM[AZ_precond] = AZ_sym_GS;
	  break;
	case Neumann:
	  optionsM[AZ_precond] = AZ_Neumann;
	  break;
	case Ls:
	  optionsM[AZ_precond] = AZ_ls;
	  break;
	case Icc:
	  optionsM[AZ_precond] = AZ_dom_decomp;
	  optionsM[AZ_subdomain_solve] = AZ_icc;
	  break;
	case Ilu:
	  optionsM[AZ_precond] = AZ_dom_decomp;
	  optionsM[AZ_subdomain_solve] = AZ_ilu;
	  break;
	case Lu:
	  optionsM[AZ_precond] = AZ_dom_decomp;
	  optionsM[AZ_subdomain_solve] = AZ_lu;
	  break;
	case Ilut:
	  optionsM[AZ_precond] = AZ_dom_decomp;
	  optionsM[AZ_subdomain_solve] = AZ_ilut;
	  break;
	case Rilu:
	  optionsM[AZ_precond] = AZ_dom_decomp;
	  optionsM[AZ_subdomain_solve] = AZ_rilu;
	  break;
	default:
	  if (Communication_Manager::My_Process_Number == 0)
	    cout << "setOption(precond_method,value):  This value is not recognized for Aztec package " 
		 << endl;
	}
      else if (option == sub_block_method)
	switch (value) {
	case Icc:
	  optionsM[AZ_subdomain_solve] = AZ_icc;
	  break;
	case Ilu:
	  optionsM[AZ_subdomain_solve] = AZ_ilu;
	  break;
	case Lu:
	  optionsM[AZ_subdomain_solve] = AZ_lu;
	  break;
	case Ilut:
	  optionsM[AZ_subdomain_solve] = AZ_ilut;
	  break;
	case Rilu:
	  optionsM[AZ_subdomain_solve] = AZ_rilu;
	  break;

	default:
	  if (Communication_Manager::My_Process_Number == 0)
	    cout << "setOption(sub_block_method,value):  This value is not recognized for Aztec package " 
		 << endl;
	}
      else 
	{
	  if (Communication_Manager::My_Process_Number == 0)
	    cout << "setOption(option,value): option is not recognized "
		 << endl;
	}
    }
  else if (packageM == PETSc)
    {
      if (option == ksp_method)
	{
          int ierr;
	  switch (value) {
	  case Cg:
	    ierr = KSPSetType(kspM, KSPCG);CHKERRCONTINUE(ierr);
	    break;
	  case Gmres:
	    ierr = KSPSetType(kspM, KSPGMRES);CHKERRCONTINUE(ierr);
	    break;
	  case Cgs:
	    ierr = KSPSetType(kspM, KSPCGS);CHKERRCONTINUE(ierr);
	    break;
	  case Tfqmr:
	    ierr = KSPSetType(kspM, KSPTFQMR);CHKERRCONTINUE(ierr);
	    break;
	  case Bcgs:
	    ierr = KSPSetType(kspM, KSPBCGS);CHKERRCONTINUE(ierr);
	    break;
	  case Direct:
	    ierr = KSPSetType(kspM, KSPPREONLY);CHKERRCONTINUE(ierr);
	    PCSetType(pcM, PCLU);
	    reInitializeM = true;
	    break;
	  case Richardson:
	    ierr = KSPSetType(kspM, KSPRICHARDSON);CHKERRCONTINUE(ierr);
	    break;
	  case Chebychev:
	    ierr = KSPSetType(kspM, KSPCHEBYSHEV);CHKERRCONTINUE(ierr);
	    break;
	  case Bicg:
	    ierr = KSPSetType(kspM, KSPBICG);CHKERRCONTINUE(ierr);
	    break;
	  case Tcqmr:
	    ierr = KSPSetType(kspM, KSPTCQMR);CHKERRCONTINUE(ierr);
	    break;
	  case Cr:
	    ierr = KSPSetType(kspM, KSPCR);CHKERRCONTINUE(ierr);
	    break;
	  case Lsqr:
	    ierr = KSPSetType(kspM, KSPLSQR);CHKERRCONTINUE(ierr);
	    break;
	  case Preonly:
	    ierr = KSPSetType(kspM, KSPPREONLY);CHKERRCONTINUE(ierr);
	    break;
	  default:
	    if (Communication_Manager::My_Process_Number == 0)
	      cout << "setOption(ksp_method,option):  This value is not recognized for PETSc package " 
		   << endl;
	  }
	}
      else if (option == precond_method)
	switch (value) {
	case Jacobi:
	  PCSetType(pcM, PCJACOBI);
	  break;
	case Icc:
	  PCSetType(pcM, PCICC);
	  break;
	case Ilu:
	  PCSetType(pcM, PCILU);
	  break;
	case Lu:
	  PCSetType(pcM, PCLU);
	  break;
	case Bjacobi:
	  PCSetType(pcM, PCBJACOBI);
	  break;
	case Sor:
	  PCSetType(pcM, PCSOR);
	  break;
	case Eisenstat:
	  PCSetType(pcM, PCEISENSTAT);
	  break;
	case Asm:
	  PCSetType(pcM, PCASM);
	  break;
	case Sles:
          PCSetType(pcM, PCKSP);
	  break;
	case Composite:
	  PCSetType(pcM, PCCOMPOSITE);
	  break;
	case Cholesky:
	  PCSetType(pcM, PCCHOLESKY);
	  break;
	case None:
	  PCSetType(pcM, PCNONE);
	  break;
	case Shell:
	  PCSetType(pcM, PCSHELL);
	  break;
	default:
	  if (Communication_Manager::My_Process_Number == 0)
	    cout << "setOption(precond_method,option):  This value is not recognized for PETSc package " 
		 << endl;
	}
      else if (option == sub_block_method)
	{
	  PCSetUp(pcM);

	  int nlocal, first;
          KSP *lsubksp;
	  PetscBool isBjacobi, isAsm;

	  PetscObjectTypeCompare((PetscObject)pcM, PCBJACOBI, &isBjacobi);
	  PetscObjectTypeCompare((PetscObject)pcM, PCASM, &isAsm);
	  
          if ((bool) isBjacobi)
            PCBJacobiGetSubKSP(pcM, &nlocal, &first, &lsubksp);
          else if ((bool) isAsm)
            PCASMGetSubKSP(pcM, &nlocal, &first, &lsubksp);

	  PC subpc;
	  for (int i=0; i<nlocal; i++)
	    {
              KSPGetPC(lsubksp[i],&subpc);
              //              KGetKSP(lsubsles[i],&subksp);	
	      switch (value) {
	      case Jacobi:
		PCSetType(subpc, PCJACOBI);
		break;
	      case Icc:
		PCSetType(subpc, PCICC);
		break;
	      case Ilu:
		PCSetType(subpc, PCILU);
		break;
	      case Lu:
		PCSetType(subpc, PCLU);
		break;
	      case Bjacobi:
		PCSetType(subpc, PCBJACOBI);
		break;
	      case Sor:
		PCSetType(subpc, PCSOR);
		break;
	      case Eisenstat:
		PCSetType(subpc, PCEISENSTAT);
		break;
	      case Asm:
		PCSetType(subpc, PCASM);
		break;
	      case Sles:
                PCSetType(subpc, PCKSP);
		break;
	      case Composite:
		PCSetType(subpc, PCCOMPOSITE);
		break;
	      case Cholesky:
		PCSetType(subpc, PCCHOLESKY);
		break;
	      case None:
		PCSetType(subpc, PCNONE);
		break;
	      default:
		if (Communication_Manager::My_Process_Number == 0)
		  cout << "setOption(sub_block_method,option):  This value is not recognized for PETSc package " 
		       << endl;
	      }
	      KSPSetType(lsubksp[i], KSPPREONLY);
	    }
	}
      else
	{
	  if (Communication_Manager::My_Process_Number == 0)
	    cout << "setOption(option,value): option is not recognized "
		 << endl;
	} 
    }
}


void
OGEquation::setParameter(OPTION parameter, double value)
{
  if (packageM == Aztec)
    {
      if (parameter == convergence)
	{
	  paramsM[AZ_conv] = value;
	}
      else 
	{
	  if (Communication_Manager::My_Process_Number == 0)
	    cout << "setParameter(parameter, value): parameter is not recognized "
		 << endl;
	}
    }
  else if (packageM == PETSc)
    {
      if (parameter == convergence)
	{
	  KSPSetTolerances(kspM, value, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	}
      else 
	{
	  if (Communication_Manager::My_Process_Number == 0)
	    cout << "setParameter(parameter, value): parameter is not recognized "
		 << endl;
	}
    }
}

void 
DMSRMatrix::getLocalToGlobal(CompositeGrid *Kimera)
{
  int count = startRow;

  for (int grid = 0; grid < Kimera->nrGrids(); grid++)
    {
      localToGlobalM[grid].redim(Kimera->rDim(grid), Kimera->sDim(grid));
      intSerialArray *myLocalToGlobal;
      myLocalToGlobal = localToGlobalM[grid].getSerialArrayPointer();
      for (int j = Kimera->lb_1(grid); j <= Kimera->ub_1(grid); j++)
	{
	  for (int i = Kimera->lb_0(grid); i <= Kimera->ub_0(grid); i++)
	    {
	      (*myLocalToGlobal)(i,j) = count;
	      count++;
	    }
	}
      localToGlobalM[grid].updateGhostBoundaries();
      //      localToGlobal[grid].display();
    }

  petscLocalRLowM = new intSerialArray[Kimera->nrGrids()];
  petscLocalSLowM = new intSerialArray[Kimera->nrGrids()];
  petscLocalRHiM = new intSerialArray[Kimera->nrGrids()];
  petscLocalSHiM = new intSerialArray[Kimera->nrGrids()];
  petscStartIndexM  = new intSerialArray[Kimera->nrGrids()];

  for (int grid=0; grid<Kimera->nrGrids(); grid++)
    {
      petscLocalRLowM[grid].redim(Communication_Manager::numberOfProcessors());
      petscLocalSLowM[grid].redim(Communication_Manager::numberOfProcessors());
      petscLocalRHiM[grid].redim(Communication_Manager::numberOfProcessors());
      petscLocalSHiM[grid].redim(Communication_Manager::numberOfProcessors());
      petscStartIndexM[grid].redim(Communication_Manager::numberOfProcessors());

      int myLowR = Kimera->lb_0(grid);
      MPI_Allgather(&myLowR, 1, MPI_INT, petscLocalRLowM[grid].getDataPointer(),
		    1, MPI_INT, MPI_COMM_WORLD);
      int myHiR = Kimera->ub_0(grid);
      MPI_Allgather(&myHiR, 1, MPI_INT, petscLocalRHiM[grid].getDataPointer(),
		    1, MPI_INT, MPI_COMM_WORLD);
      int myLowS = Kimera->lb_1(grid);
      MPI_Allgather(&myLowS, 1, MPI_INT, petscLocalSLowM[grid].getDataPointer(),
		    1, MPI_INT, MPI_COMM_WORLD);
      int myHiS = Kimera->ub_1(grid);
      MPI_Allgather(&myHiS, 1, MPI_INT, petscLocalSHiM[grid].getDataPointer(),
		    1, MPI_INT, MPI_COMM_WORLD);
      intSerialArray *myLocalToGlobal;
      myLocalToGlobal = localToGlobalM[grid].getSerialArrayPointer();
      int myStartIndex = (*myLocalToGlobal)(Kimera->lb_0(grid),Kimera->lb_1(grid));
      MPI_Allgather(&myStartIndex, 1, MPI_INT, petscStartIndexM[grid].getDataPointer(),
		    1, MPI_INT, MPI_COMM_WORLD);
    }
}

int 
DMSRMatrix::ownerProcess(int i, int j, int grid)
{
  for (int proc = 0; proc < Communication_Manager::numberOfProcessors(); proc++)
    {
      if ( i >= petscLocalRLowM[grid](proc) && i <= petscLocalRHiM[grid](proc) && 
	   j >= petscLocalSLowM[grid](proc) && j <= petscLocalSHiM[grid](proc) )
	{
	  return proc;
	}
    }
  cerr << " Something is wrong in DMSRMatrix::ownerProcess() , I am "
       << Communication_Manager::My_Process_Number 
       << ", i,j is " << i << " ," << j << endl;
  return -1;
}
    
