// -*-c++-*-
#ifndef GRIDFUNC_H
#define GRIDFUNC_H

#include "A++.h"

#include "CompGrid.hh"

enum functionType
{
  scalar, xComponentOfVector, yComponentOfVector, DirichletAllOver, ExtrapolationAllOver
};

enum interpolationType
{
  linear, quadratic, cubic
};

class gridFunction 
{
  
  friend class OGEquation;

private :
  
  bool conservativeDifferenceOperatorsM;

  int nmbrOfDiscrPointsM;

  interpolationType iTypeM;
  
  CompositeGrid* myGridM;

  void dirichlet(int i, int k);

  void neumann(int i, int k);

  void robin(int i, int k);

  void periodic(int i, int k);

  doubleSerialArray computeInterpolationPoints(int k);

  doubleSerialArray alpha(int i,int k);

  doubleSerialArray beta(int j,int k);  

  doubleSerialArray gamma(int i,int k);

  doubleSerialArray delta(int j,int k);  

  doubleSerialArray epsilon(int i,int k);

  doubleSerialArray zeta(int j,int k);

  void fillIntElements(doubleSerialArray& elementArray, gridFunction *target, int grid=-1);

  doubleArray (*boundaryValue) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side);
  doubleArray (*boundaryValue_d_dt) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side);

  doubleArray boundaryDiff(int i, int gType, const double t, const doubleArray& x, const doubleArray& y, const doubleArray& xr, const doubleArray& xs, const doubleArray& yr, const doubleArray& ys);

  doubleArray (*bc_x) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side);

  doubleArray (*bc_y) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side);

  doubleArray (*InitFunction) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side);

  double timeM;

  BC *bcsM;

  doubleArray *fieldValuesM;

public :

  void setBoundaryConditionType(int grid, Side side, boundaryType condition);
 
  void setUseConservativeDifferenceOperators(){conservativeDifferenceOperatorsM = true;};

  void dirichlet(int i, int k, gridFunction& F);

  boundaryType getBoundaryType(int grid, Side side);

  void setInterpolationType(interpolationType arg) {iTypeM = arg;};

  BC* getBC() { return bcsM;};

  void setBC(functionType fType);

  void setDivergenceToZeroAtBoundaries(gridFunction& otherVelocityComponent, 
				       functionType fType);

  doubleArray& operator()(int k) { 
    return fieldValuesM[k];
  }

  gridFunction& operator=(const double right) {
    for (int k = 0; k < myGridM->nrGrids(); k++) {
      Index ix = fieldValuesM[k].getFullRange(0);
      Index iy = fieldValuesM[k].getFullRange(1);
      fieldValuesM[k](ix, iy) = right;
    }
    return *this;
  }

  gridFunction& operator=(const gridFunction& right) {
    if (this == &right) return *this;
    for (int k = 0; k < myGridM->nrGrids(); k++) {
      Index ix = right.fieldValuesM[k].getFullRange(0);
      Index iy = right.fieldValuesM[k].getFullRange(1);
      fieldValuesM[k](ix, iy) = right.fieldValuesM[k];
    }
    return *this;
  }

  gridFunction& operator+=(const gridFunction& right) {
    for (int k = 0; k < myGridM->nrGrids(); k++) {
      Index ix = right.fieldValuesM[k].getFullRange(0);
      if (ix.getBase() < fieldValuesM[k].getBase(0) || ix.getBound() > fieldValuesM[k].getBound(0)) {
	ix = fieldValuesM[k].getFullRange(0);
      }
      Index iy = right.fieldValuesM[k].getFullRange(1);
      if (iy.getBase() < fieldValuesM[k].getBase(1) || iy.getBound() > fieldValuesM[k].getBound(1)) {
	iy = fieldValuesM[k].getFullRange(1);
      }
      fieldValuesM[k](ix, iy) += right.fieldValuesM[k](ix, iy);
    }
    return *this;
  }
  
  gridFunction& operator+=(const double right) {
    for (int k = 0; k < myGridM->nrGrids(); k++) {
      Index ix = fieldValuesM[k].getFullRange(0);
      Index iy = fieldValuesM[k].getFullRange(1);
      
      fieldValuesM[k](ix, iy) += right;
    }
    return *this;
  }

  gridFunction& operator-=(const gridFunction& right) {
    for (int k = 0; k < myGridM->nrGrids(); k++) {
      Index ix = right.fieldValuesM[k].getFullRange(0);
      if (ix.getBase() < fieldValuesM[k].getBase(0) || ix.getBound() > fieldValuesM[k].getBound(0)) {
	ix = fieldValuesM[k].getFullRange(0);
      }
      Index iy = right.fieldValuesM[k].getFullRange(1);
      if (iy.getBase() < fieldValuesM[k].getBase(1) || iy.getBound() > fieldValuesM[k].getBound(1)) {
	iy = fieldValuesM[k].getFullRange(1);
      }
      fieldValuesM[k](ix, iy) -= right.fieldValuesM[k](ix, iy);
    }
    return *this;
  }

  gridFunction& operator*=(const gridFunction& right) {
    for (int k = 0; k < myGridM->nrGrids(); k++) {
      Index ix = right.fieldValuesM[k].getFullRange(0);
      if (ix.getBase() < fieldValuesM[k].getBase(0) || ix.getBound() > fieldValuesM[k].getBound(0)) {
	ix = fieldValuesM[k].getFullRange(0);
      }
      Index iy = right.fieldValuesM[k].getFullRange(1);
      if (iy.getBase() < fieldValuesM[k].getBase(1) || iy.getBound() > fieldValuesM[k].getBound(1)) {
	iy = fieldValuesM[k].getFullRange(1);
      }
      fieldValuesM[k](ix, iy) *= right.fieldValuesM[k](ix, iy);
    }
    return *this;
  }

  gridFunction& operator/=(const gridFunction& right) {
    for (int k = 0; k < myGridM->nrGrids(); k++) {
      Index ix = right.fieldValuesM[k].getFullRange(0);
      if (ix.getBase() < fieldValuesM[k].getBase(0) || ix.getBound() > fieldValuesM[k].getBound(0)) {
	ix = fieldValuesM[k].getFullRange(0);
      }
      Index iy = right.fieldValuesM[k].getFullRange(1);
      if (iy.getBase() < fieldValuesM[k].getBase(1) || iy.getBound() > fieldValuesM[k].getBound(1)) {
	iy = fieldValuesM[k].getFullRange(1);
      }
      fieldValuesM[k](ix, iy) /= right.fieldValuesM[k](ix, iy);
    }
    return *this;
  }

  const gridFunction& operator+() const {
    return *this;
  }

  const gridFunction operator-() const {
    gridFunction result = *this;
    for (int k = 0; k < this->myGridM->nrGrids(); k++) {
      result.fieldValuesM[k] = - (this->fieldValuesM[k]);
    }
    return result;
  }

  intArray* flagValuesM;

  gridFunction(CompositeGrid& Kimera, double time = 0);

  gridFunction(const gridFunction& fun);

  gridFunction(int i);

  ~gridFunction();

  void printValues(char* info);

  void saveToFile(const char* name);

  void readFromFile(const char* name);

  void readFromHDF5File(const char* fileName, const char* varName);

  void createHDF5File(const char* name);

  void saveToHDF5File(const char* fileName, const char* varName);

  void initializeGridData();

  doubleArray x(int k, Index ix = -1000, Index iy = -1000) const; 

  doubleArray y(int k, Index ix = -1000, Index iy = -1000) const; 

  doubleArray xx(int k, Index ix = -1000, Index iy = -1000) const;

  doubleArray yy(int k, Index ix = -1000, Index iy = -1000) const;

  doubleArray xy(int k, Index ix = -1000, Index iy = -1000) const;

  doubleArray laplacian(int k, Index ix = -1000, Index iy = -1000) const;

  gridFunction x() const; 

  gridFunction y() const; 

  gridFunction xx() const;

  gridFunction yy() const;

  gridFunction xy() const;

  gridFunction laplacian() const;

  void updateBCs(boundaryType rv = ALL);

  void finishBCs();

  void interpolate(int grid = -1, gridFunction* target = 0x0);

  void updateGhostBoundaries();

  double get_CFL_Dt(const double& cfl, const gridFunction& a, const gridFunction& b, CompositeGrid& Kimera);
  
  double getDt(const double& cfl, const gridFunction& a, const gridFunction& b, const double visc, CompositeGrid& Kimera, const double alpha0 = -2, const double beta0 = 1);

  double getDt(const double& cfl, const gridFunction& a, const gridFunction& b, const gridFunction& visc, CompositeGrid& Kimera, const double alpha0 = -2, const double beta0 = 1);

  void applyMask(int k);

  void addTime(double time);

  double getTime() const;

  void checkErrors(doubleArray (*infunc) (const doubleArray& x, const doubleArray& y, const double t, int, int));

  doubleArray forcing(doubleArray (*infunc) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side), int grid, double tadd = 0., Index ix = -1000, Index iy = -1000);

  doubleArray boundary_d_dt(int grid, double tadd = 0., Index ix = -1000, Index iy = -1000) { return forcing(boundaryValue_d_dt, grid, tadd, ix, iy); }

  gridFunction forcing(doubleArray (*infunc) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side), double tadd = 0.);
  
  void set_Dirichlet(doubleArray (*infunc) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side)) {boundaryValue = infunc;}

  void set_Dirichlet_time_derivative(doubleArray (*infunc) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side)) {boundaryValue_d_dt = infunc;}

  void set_bc_x(doubleArray (*infunc) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side)) {bc_x = infunc;}

  void set_bc_y(doubleArray (*infunc) (const doubleArray& x, const doubleArray& y, const double t, int grid, int side)) {bc_y = infunc;}

  void setFlagValues();
  
  Range getDiscrBounds(int direction, int k) const;

  void extrapolateToBoundaries();
  
  int nrGrids() const {return myGridM->nrGrids();};

  double l_inf();

  double l_2();

  double l_h();
};

namespace GridFunction {
  doubleArray zero(const doubleArray& x, 
                   const doubleArray& y, 
                   const double t, 
                   int grid, 
                   int side);
};

const gridFunction operator+(const gridFunction& left, const gridFunction& right);

const gridFunction operator-(const gridFunction& left, const gridFunction& right);

const gridFunction operator*(const gridFunction& left, const gridFunction& right);

const gridFunction operator/(const gridFunction& left, const gridFunction& right);

const gridFunction operator*(const double left, const gridFunction& right);

const gridFunction operator*(const gridFunction& left, const double right);

const gridFunction operator/(const gridFunction& left, const double right);

#endif
