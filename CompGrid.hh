// -*-c++-*-
#ifndef COMPGRID_H

#define COMPGRID_H

#include <vector>
#include "A++.h"

using std::vector;

enum FileType
{
  ASCII, HDF5
};

enum
{
  I,J
};

enum Side
{
  lowR=0, hiR=1, lowS=2, hiS=3
};

enum boundaryType
{
  NOSLIP, SLIP, INFLOW, OUTFLOW, OUTFLOW_NEUM_P, SYMMETRY, DIRICHLET, INTERPOLATION, PERIODIC, NEUMANN, ROBIN, EXTRAPOLATION, ALL
};

// MAGIC NUMBERS IN XCOG-FILE
// 0 == interpolation
// 1 == noslip or inflow
// 3 == periodic
// 5 == slip
// 20 == inflow (not sure about this one)
// 21 == extrapolated outflow with mixed (Robin) pressure
// 22 == extrapolated outflow with Neumann pressure
// 50 == dirichlet (twilight zone)  
class BC
{
public:
  
  boundaryType lowR, hiR, lowS, hiS;
  
  int bP[4];
  BC() {};
  ~BC() {};
};

class CompositeGrid 
{
  //////////////////////////////
  friend class gridFunction;
  friend class OGEquation;
  friend class DMSRMatrix;
  
private :
  
  doubleArray  *xM, *yM, *xrM, *xsM, *yrM, *ysM, *xrrM, *xssM, *yrrM, *yssM, *sqrtOfGM, *maskM;
  //  std::vector<doubleArray> xM;

  vector<int> lowerbound_0M, lowerbound_1M, upperbound_0M, upperbound_1M;
  
  intSerialArray *myIntPointsM, *my_donor_pointM, *intInterpolationLocationM, myNrOfIntPointsM, nrOfPointsToComputeM, offsetsM, *My_i_LM, *My_j_LM;
  
  doubleSerialArray *my_donor_parameterM, *doubleInterpolationLocationM, *interpolationCoordinatesM, theSendBufferM, theReceiveBufferM;
  
  MPI_Datatype *receiveTypeM, *sendTypeM;
  
  double *r_stepM, *s_stepM;
  
  intArray *flagValuesM;
  
  vector<int> gridTypeM, rdimM, sdimM;

  int nmbrOfGridsM;

  int boundaryWidthM;

  BC *bcsM;

  Partitioning_Type *ddM;

  void getLocalInterp(int *nr_int_point, intSerialArray *i_point, 
                      intSerialArray *j_point, intSerialArray *i_loc, 
                      intSerialArray *j_loc, intSerialArray *gridLoc, 
                      doubleSerialArray *r_loc, doubleSerialArray *s_loc);

  void gridSort(intSerialArray *i_point, intSerialArray *j_point, 
                intSerialArray *i_loc, intSerialArray *j_loc, 
                intSerialArray *gridLoc, doubleSerialArray *r_loc, 
                doubleSerialArray *s_loc);

  void setBoundaryType();

  void checkIfGridIsRegular(int grid);

public :

  boundaryType getBoundaryType(int grid, Side side) const;

  int rDim(int grid) const;
  int sDim(int grid) const;
  int gType(int grid) const;
  double rStep(int grid) const;
  double sStep(int grid) const;

  void saveCoordinatesToFile() const;

  void saveCoordinatesToHDF5File(const char *name) const;

  void setGhostCellWidth();

  CompositeGrid(std::string xcogFileName, 
                Partitioning_Type gridDistribution[] = 0x0);

  CompositeGrid(int dim1, int dim2, 
                Partitioning_Type gridDistribution[] = 0x0);

  ~CompositeGrid();

  int nrGrids() const;

  Range getBounds(int direction, int k) const;

  void updateGhostBoundaries();

  void printIntInf(int k) const;

  int nrOfLocalIntPoints(int k) const { return (myIntPointsM[k]).getLength(0); };

  int getInterpIndex(int i, int j, int grid) const;

  int lb_0(int i) const { return lowerbound_0M[i]; };
  int lb_1(int i) const { return lowerbound_1M[i]; };
  int ub_0(int i) const { return upperbound_0M[i]; };
  int ub_1(int i) const { return upperbound_1M[i]; };

  doubleArray normalVector_R_Component(int grid, int side, 
                                       Index ixb=-1000, Index low_r=-1000, 
                                       Index hi_r=-1000, Index iyb=-1000, 
                                       Index low_s=-1000, Index hi_s=-1000) const;

  doubleArray normalVector_S_Component(int grid, int side, 
                                       Index ixb=-1000, Index low_r=-1000, 
                                       Index hi_r=-1000, Index iyb=-1000, 
                                       Index low_s=-1000, Index hi_s=-1000) const;

  doubleArray gridSize(int grid, int direction1, 
                       int direction2, Index ix=-1000, Index iy=-1000) const;

  double hSquare() const;

  void simplifyFlagValues();
};

#endif
