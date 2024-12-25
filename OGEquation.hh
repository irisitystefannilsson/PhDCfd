// -*-c++-*-
#ifndef OGEQUATION_H
#define OGEQUATION_H

#include "petscksp.h"
extern "C"
{
#include "az_aztec.h"
}
#undef min
#undef max

#include "A++.h"

#include "GridFunc.hh"
#include "CompGrid.hh"


enum MatrixOperator 
{
  Poisson, Poisson_with_TD_diagonal, Grad_Visc_Div, Grad_Visc_Div_with_TD_diagonal, PUM, PUM_with_TD_diagonal, Conv_varDiff , Conv_varDiff_with_TD_diagonal
};

enum sparseSolver
{
  Aztec, PETSc
};

enum OPTION
{
  ksp_method, precond_method, sub_block_method, convergence
};

enum SOLVER
{
  Richardson, Chebychev, Cg, Bicg, Gmres, Bcgs, Cgs, Tfqmr, Tcqmr, Cr, Lsqr, Preonly, Direct
};

enum PRECONDITIONER
{
  Jacobi, Bjacobi, Sor, Eisenstat, Icc, Ilu, Ilut, Rilu, Asm, Sles, Composite, Lu, Cholesky, None, Shell, Neumann, Ls, Symgs
};


const double TOLERANCEM = 1.e-12;
const double rtolM = 1.e-30;
const double ConvConstM = 0.001;

class DMSRMatrix
{
private :
  // Variables used to insert values in
  // PETSc Matrix
  PetscScalar vM[21];
  int mM, nM;
  int idxmM[1], idxnM[21];
  int vidxM;

  intArray *localToGlobalM;
  intSerialArray *petscIndexM;
  intSerialArray *petscLocalRLowM;
  intSerialArray *petscLocalSLowM;
  intSerialArray *petscLocalRHiM;
  intSerialArray *petscLocalSHiM;
  intSerialArray *petscStartIndexM;

  sparseSolver matTypeM;

  int NrOfNonZerosM;
  
  int globalNrOfUnknownsM;

  double *aM;
  double *rowScalingM;
  int *jaM;
  doubleSerialArray CoeffSendBufferM;
  doubleSerialArray CoeffReceiveBufferM;
  
  void getLocalToGlobal(CompositeGrid *Kimera);

  void insertInterpolationInfo(int localrow, int *update, int grid, const int fromGrid, int i, int j, CompositeGrid *Kimera, interpolationType iType, sparseSolver package);

  void insertFD(int row, int localrow, int rdim, int grid, int i, int j, MatrixOperator operatorType, int gridType, doubleSerialArray *Xr, doubleSerialArray *Xs, doubleSerialArray *Yr, doubleSerialArray *Ys, doubleSerialArray *Xrr, doubleSerialArray *Xss, doubleSerialArray *Yrr, doubleSerialArray *Yss, doubleSerialArray *sqrtOfG, doubleSerialArray *viscosity, double r_step, double s_step, sparseSolver package, doubleSerialArray *convFieldU, doubleSerialArray *convFieldV, doubleSerialArray *X=0, doubleSerialArray *Y=0, CompositeGrid *Kimera=0 );

  void insertNeumann(int row, int localrow, int side, int rdim, int i, int j, int gridType, doubleSerialArray *Xr, doubleSerialArray *Xs, doubleSerialArray *Yr, doubleSerialArray *Ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, sparseSolver package);

  void insertRobin(int row, int localrow, int side, int rdim, int i, int j, int gridType, doubleSerialArray *Xr, doubleSerialArray *Xs, doubleSerialArray *Yr, doubleSerialArray *Ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, sparseSolver package);

  void insertPeriod(int row, int localrow, int side, int rdim, int sdim, sparseSolver package, int grid, int i, int j);

  void insertCorner(int row, int localrow, int rdim, int corner, sparseSolver package, int i,int j);
  
  double xxCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *xrr, doubleSerialArray *xss, doubleSerialArray *yrr, doubleSerialArray *yss, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3]);
  
  double yyCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *xrr, doubleSerialArray *xss, doubleSerialArray *yrr, doubleSerialArray *yss, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3]);

  double xCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3]);
  
  double upWindxCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3], doubleSerialArray *cF);
  
  double yCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3]);

  double upWindyCoeff(int i, int j, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3], doubleSerialArray *cF);

  double NeumannCoeff(int i, int j, int side, doubleSerialArray *xr, doubleSerialArray *xs, doubleSerialArray *yr, doubleSerialArray *ys, doubleSerialArray *sqrtOfG, double r_step, double s_step, const double coeff[3][3]);

  void getInterpolationCoefficients(CompositeGrid *Kimera, interpolationType iType);
  
  doubleSerialArray computeInterpolationWeights(int k, CompositeGrid *Kimera, interpolationType iType);
  
  doubleSerialArray & alpha(int i,int k, CompositeGrid *Kimera);
  
  doubleSerialArray & beta(int i,int k, CompositeGrid *Kimera);
  
  doubleSerialArray & gamma(int i,int k, CompositeGrid *Kimera);
  
  doubleSerialArray & delta(int i,int k, CompositeGrid *Kimera);

  doubleSerialArray & epsilon(int i,int k, CompositeGrid *Kimera);
  
  doubleSerialArray & zeta(int i,int k, CompositeGrid *Kimera);

  int ownerProcess(int i, int j, int grid);

public :

  int startRow, endRow;

  void setupMatrix(CompositeGrid *Kimera, MatrixOperator operatorType, int N_update, int *update, intArray flagValues[], bool allNeumann, interpolationType iType, gridFunction *visc, bool reset, sparseSolver package, gridFunction *convFieldU, gridFunction *convFieldV);

  void addToDiagonal(CompositeGrid *Kimera, const double dt, intArray flagValues[], int update_index[], int update[], sparseSolver package);
  
  DMSRMatrix();
  
  ~DMSRMatrix();

  Mat PmatM;
  
  double rowScale(int row) { return rowScalingM[row]; };

  double *val() { return aM; };
  
  int *bindx() { return jaM; };
  
  int getGlobalNrOfUnknowns() { return globalNrOfUnknownsM; }; 
  
};

class OGEquation
{
private :
  
  //  SLES    slesM; 
  //  SLES    *subslesM;
  KSP kspM;
  KSP *subkspM;
  PC pcM;
  
  int nlocalM, firstM;

  bool reInitializeM;
  bool setFromOptionsM;

  sparseSolver packageM;

  AZ_MATRIX *AmatM;

  AZ_PRECOND *precM;

  MatrixOperator equationM;
  
  DMSRMatrix MatrixM;
  
  CompositeGrid *OGridM;
  
  gridFunction *unKnownM;
  double *lhsM, *rhsM;
  Vec xM, bM;

  int proc_configM[AZ_PROC_SIZE];

  int optionsM[AZ_OPTIONS_SIZE];

  double paramsM[AZ_PARAMS_SIZE];

  int *data_orgM;             

  double statusM[AZ_STATUS_SIZE];
                                
  int *updateM, *externalM;                  

  int *update_indexM;
         
  int *extern_indexM; 

  int N_updateM;

public :
  
  void setBoundaryValues(gridFunction & RHS, CompositeGrid *Kimera);
  
  OGEquation(CompositeGrid *Kimera, sparseSolver spSolver=Aztec, int itm=AZ_bicgstab, int prec=AZ_dom_decomp, int over=0);

  ~OGEquation(); 

  void destroyOperator();

  void solve(gridFunction & RHS, CompositeGrid *Kimera, double dt=0, bool allNeumann=false, bool setBCs=true);
  
  void setUnknown(gridFunction & LHS) { unKnownM = &LHS; if (packageM == Aztec) {rhsM = new double[N_updateM+1];} };//+1 if we need an extra equation for a pure Neumann problem
  
  void setOperator(MatrixOperator operatorType, bool allNeumann=false, gridFunction *visc=0, bool reset=false, gridFunction *convFieldU=0, gridFunction *convFieldV=0);

  void infNorm();

  void writeMatrixOperatorToFile(const char *name);
  
  void setOption(OPTION option, int value);

  void setParameter(OPTION parameter, double value);
};
#endif
