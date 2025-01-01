#include <mpi.h>
#include "hdf5.h"
#include <assert.h>

#include <algorithm>

#include "CompGrid.hh"

#define FAIL -1

using std::vector;

CompositeGrid::CompositeGrid(std::string xcogFileName, Partitioning_Type gridDistribution[]) {
  fromFileM = true;
  int *nr_int_point = 0x0, nrProcs, lineSkip, dim1, dim2;
  intSerialArray *i_point = 0x0, *j_point = 0x0, *i_loc = 0x0, *j_loc = 0x0, *gridLoc = 0x0;
  doubleSerialArray *r_loc = 0x0, *s_loc = 0x0;
  char word[200];
  std::ifstream xcogFile;
  MPI_Comm_size(MPI_COMM_WORLD, &nrProcs);
  
  FileType typ;
  if (xcogFileName.find(".acg") != std::string::npos) typ = ASCII;
  else if (xcogFileName.find(".h5") != std::string::npos) typ = HDF5;
  else throw "xcogFile is apparently neither of ASCII nor of hdf5 type!\n";
      
  if (typ == ASCII) {
    xcogFile.open(xcogFileName.c_str(), ios::in);
    xcogFile >> word;
    while (strcmp(word,"@begin_component_grid_information")) {
      if (!strcmp(word,"@n_grids"))
	xcogFile >> nmbrOfGridsM;
      xcogFile >> word;
    }
    ddM = new Partitioning_Type[nmbrOfGridsM];
    if (gridDistribution != 0) {
      for (int dist=0; dist<nmbrOfGridsM; dist++) {
	ddM[dist] = gridDistribution[dist];
      }
    } else { 
      for (int dist=0; dist<nmbrOfGridsM; dist++) {
	ddM[dist] = Partitioning_Type();
      }
    }
    /*
      Local arrays used to setup interpolation
      arrays for all processors
    */
    nr_int_point = new int[nmbrOfGridsM];
    i_point = new intSerialArray[nmbrOfGridsM];
    j_point = new intSerialArray[nmbrOfGridsM];
    i_loc = new intSerialArray[nmbrOfGridsM];
    j_loc = new intSerialArray[nmbrOfGridsM];
    gridLoc = new intSerialArray[nmbrOfGridsM];
    r_loc = new doubleSerialArray[nmbrOfGridsM];
    s_loc = new doubleSerialArray[nmbrOfGridsM];
    
    xM = new doubleArray[nmbrOfGridsM];
    yM = new doubleArray[nmbrOfGridsM];
    flagValuesM = new intArray[nmbrOfGridsM];
    maskM = new doubleArray[nmbrOfGridsM];
    xrM = new doubleArray[nmbrOfGridsM];
    xsM = new doubleArray[nmbrOfGridsM];
    yrM = new doubleArray[nmbrOfGridsM];
    ysM = new doubleArray[nmbrOfGridsM];
    xrrM = new doubleArray[nmbrOfGridsM];
    xssM = new doubleArray[nmbrOfGridsM];
    yrrM = new doubleArray[nmbrOfGridsM];
    yssM = new doubleArray[nmbrOfGridsM];
    sqrtOfGM = new doubleArray[nmbrOfGridsM];
    
    myNrOfIntPointsM.redim(nmbrOfGridsM, nrProcs);
    nrOfPointsToComputeM.redim(nmbrOfGridsM,nmbrOfGridsM,nrProcs);
    
    offsetsM.redim(nrProcs,nmbrOfGridsM,nmbrOfGridsM);
    
    receiveTypeM = new MPI_Datatype[nrProcs];
    sendTypeM = new MPI_Datatype[nrProcs];
    
    myIntPointsM = new intSerialArray[nmbrOfGridsM];
    My_i_LM = new intSerialArray[nmbrOfGridsM];
    My_j_LM = new intSerialArray[nmbrOfGridsM];
    intInterpolationLocationM = new intSerialArray[nmbrOfGridsM];
    doubleInterpolationLocationM = new doubleSerialArray[nmbrOfGridsM];
    interpolationCoordinatesM = new doubleSerialArray[nmbrOfGridsM];
    r_stepM = new double[nmbrOfGridsM];
    s_stepM = new double[nmbrOfGridsM];
    gridTypeM.resize(nmbrOfGridsM);
    rdimM.resize(nmbrOfGridsM);
    sdimM.resize(nmbrOfGridsM);
    bcsM = new BC[nmbrOfGridsM];
      
    for (int k = nmbrOfGridsM - 1; k >= 0; k--) {
      xcogFile >> word;
      
      while (strcmp(word, "@end_component_grid")) {
	if (!strcmp(word,"@grid_type"))
	  xcogFile >> gridTypeM[k];
	else if(!strcmp(word,"@n_interp"))
	  xcogFile >> nr_int_point[k];
	else if(!strcmp(word,"@r_dim"))
	  xcogFile >> rdimM[k];
	else if(!strcmp(word,"@s_dim"))
	  xcogFile >> sdimM[k];
	else if(!strcmp(word,"@r_period"))
	  xcogFile >> word;
	else if(!strcmp(word,"@s_period"))
	  xcogFile >> word;
	else if(!strcmp(word,"@first_point"))
	  xcogFile >> word;
	else if(!strcmp(word,"@last_point"))
	  xcogFile >> word;
	else if(!strcmp(word,"@x_min"))
	  xcogFile >> word;
	else if(!strcmp(word,"@x_max"))
	  xcogFile >> word;
	else if(!strcmp(word,"@y_min"))
	  xcogFile >> word;
	else if(!strcmp(word,"@y_max"))
	  xcogFile >> word;
	else if(!strcmp(word,"@r_step"))
	  xcogFile >> r_stepM[k];
	else if(!strcmp(word,"@s_step"))
	  xcogFile >> s_stepM[k];
	else if(!strcmp(word,"@int_array_1d")) {
	  xcogFile >> word;	
	  if (!strcmp(word,"i_point")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    if (dim1>0) {
	      i_point[k].redim(dim1);
	      for (lineSkip=0; lineSkip<8; lineSkip++)
		xcogFile >> word;
	      for (int i=0; i<dim1; i++)
		xcogFile >> i_point[k](i);
	    }
	  } else if (!strcmp(word,"j_point")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    if (dim1>0) {
	      j_point[k].redim(dim1);
	      for (lineSkip=0; lineSkip<8; lineSkip++)
		xcogFile >> word;
	      for (int i=0; i<dim1; i++)
		xcogFile >> j_point[k](i);
	    }
	  } else if (!strcmp(word,"i_loc")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    if (dim1>0) {
	      i_loc[k].redim(dim1);
	      for (lineSkip=0; lineSkip<8; lineSkip++)
		xcogFile >> word;
	      for (int i=0; i<dim1; i++)
		xcogFile >> i_loc[k](i);
	    }
	  } else if (!strcmp(word,"j_loc")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    if (dim1>0) {
	      j_loc[k].redim(dim1);
	      for (lineSkip=0; lineSkip<8; lineSkip++)
		xcogFile >> word;
	      for (int i=0; i<dim1; i++)
		xcogFile >> j_loc[k](i);
	    }
	  } else if (!strcmp(word,"grid_loc->priority")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    if (dim1>0) {
	      gridLoc[k].redim(dim1);
	      for (lineSkip=0; lineSkip<8; lineSkip++)
		xcogFile >> word;
	      for (int i=0; i<dim1; i++)
		xcogFile >> gridLoc[k](i);
	    }
	  }
	} else if(!strcmp(word,"@int_array_2d")) {
	  xcogFile >> word;
	  if (!strcmp(word,"range")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    xcogFile >> word;
	    xcogFile >> dim2;
	    for (lineSkip=0; lineSkip<12; lineSkip++)
	      xcogFile >> word;
	  } else if(!strcmp(word,"bc")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    xcogFile >> word;
	    xcogFile >> dim2;
	    for (lineSkip=0; lineSkip<8; lineSkip++)
	      xcogFile >> word;
	    xcogFile >> bcsM[k].bP[0];
	    xcogFile >> bcsM[k].bP[1];
	    xcogFile >> bcsM[k].bP[2];
	    xcogFile >> bcsM[k].bP[3];
	  } else if(!strcmp(word,"flag")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    xcogFile >> word;
	    xcogFile >> dim2;
	    for (lineSkip=0; lineSkip<8; lineSkip++)
	      xcogFile >> word;
	    
	    flagValuesM[k].partition(ddM[k]);
	    
	    flagValuesM[k].redim(dim1,dim2);
	    
	    intSerialArray *localFlag = flagValuesM[k].getSerialArrayPointer();
	    int i_low_corner = flagValuesM[k].getLocalBase(0);
	    int i_hi_corner =  flagValuesM[k].getLocalBound(0);
	    int j_low_corner = flagValuesM[k].getLocalBase(1);
	    int j_hi_corner =  flagValuesM[k].getLocalBound(1);
	    for (int j = 0; j < j_low_corner; j++)
	      for (int i = 0; i < dim1; i++)
		xcogFile >> word;
	    for (int j1 = j_low_corner; j1 <= j_hi_corner; j1++) {
	      for (int i = 0; i < i_low_corner; i++) 
		xcogFile >> word;
	      for (int i1 = i_low_corner; i1 <= i_hi_corner; i1++) 
		xcogFile >> (*localFlag)(i1, j1);
	      for (int i2 = i_hi_corner + 1; i2 < dim1; i2++) 
		xcogFile >> word;
	    }
	  }
	} else if (!strcmp(word, "@real_array_1d")) {
	  xcogFile >> word;
	  if (!strcmp(word, "r_loc")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    if (dim1>0) {
	      r_loc[k].redim(dim1);
	      for (lineSkip = 0; lineSkip < 8; lineSkip++)
		xcogFile >> word;
	      for (int i = 0; i < dim1; i++)
		xcogFile >> r_loc[k](i);
	    }
	  } else if(!strcmp(word, "s_loc")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    if (dim1 > 0) {
	      s_loc[k].redim(dim1);
	      for (lineSkip = 0; lineSkip < 8; lineSkip++)
		xcogFile >> word;
	      for (int i = 0; i < dim1; i++)
		xcogFile >> s_loc[k](i);
	    }
	  }
	} else if (!strcmp(word, "@real_array_2d")) {
	  xcogFile >> word;
	  if(!strcmp(word, "x")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    xcogFile >> word;
	    xcogFile >> dim2;
	    for (lineSkip = 0; lineSkip < 8; lineSkip++)
	      xcogFile >> word;
	    
	    xM[k].partition(ddM[k]);
	    
	    xM[k].redim(dim1, dim2);
	    
	    doubleSerialArray *localX;
	    localX = xM[k].getSerialArrayPointer();
	    int i_low_corner = xM[k].getLocalBase(0);
	    int i_hi_corner =  xM[k].getLocalBound(0);
	    int j_low_corner = xM[k].getLocalBase(1);
	    int j_hi_corner =  xM[k].getLocalBound(1);
	    for (int j=0; j < j_low_corner; j++)
	      for (int i=0; i < dim1; i++)
		xcogFile >> word;
	    for (int j1=j_low_corner; j1 <= j_hi_corner; j1++) {
	      for (int i = 0; i < i_low_corner; i++) 
		xcogFile >> word;
	      for (int i1 = i_low_corner; i1 <= i_hi_corner; i1++) 
		xcogFile >> (*localX)(i1, j1);
	      for (int i2 = i_hi_corner + 1; i2 < dim1; i2++) 
		xcogFile >> word;
	    }
	  } else if(!strcmp(word, "y")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    xcogFile >> word;
	    xcogFile >> dim2;
	    for (lineSkip = 0; lineSkip < 8; lineSkip++)
	      xcogFile >> word;
	    
	    yM[k].partition(ddM[k]);
	    
	    yM[k].redim(dim1, dim2);
	    
	    doubleSerialArray *localY;
	    localY = yM[k].getSerialArrayPointer();
	    int i_low_corner = yM[k].getLocalBase(0);
	    int i_hi_corner =  yM[k].getLocalBound(0);
	    int j_low_corner = yM[k].getLocalBase(1);
	    int j_hi_corner =  yM[k].getLocalBound(1);
	    for (int j = 0; j < j_low_corner; j++)
	      for (int i = 0; i < dim1; i++)
		xcogFile >> word;
	    for (int j1 = j_low_corner; j1 <= j_hi_corner; j1++) {
	      for (int i = 0; i < i_low_corner; i++) 
		xcogFile >> word;
	      for (int i1 = i_low_corner; i1 <= i_hi_corner; i1++) 
		xcogFile >> (*localY)(i1,j1);
	      for (int i2 = i_hi_corner + 1; i2 < dim1; i2++) 
		xcogFile >> word;
	    }
	    
	    if (gridTypeM[k] == 1) {
	      //----------------------------------------
	      // This means grid is at least Cartesian
	      // , and maybe even regular
	      //----------------------------------------
	      //----------------------------------------
	      // Control whether the grid really is regular
	      // or not. Both x and y are known now (I hope)
	      //----------------------------------------
	      
	      checkIfGridIsRegular(k);
	    } else {
	      //----------------------------------------
	      // Otherwise it is Curvilinear
	      //----------------------------------------
	      gridTypeM[k] = 3;
	    }
	  } else if (!strcmp(word, "xr")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    xcogFile >> word;
	    xcogFile >> dim2;
	    for (lineSkip = 0; lineSkip < 8; lineSkip++)
	      xcogFile >> word;
	    
	    if (gridTypeM[k] != 1) {  // curvilinear 	  
		
	      xrM[k].partition(ddM[k]);
	      
	      xrM[k].redim(dim1, dim2);
	      
	      doubleSerialArray *localXr;
	      localXr = xrM[k].getSerialArrayPointer();
	      int i_low_corner = xrM[k].getLocalBase(0);
	      int i_hi_corner =  xrM[k].getLocalBound(0);
	      int j_low_corner = xrM[k].getLocalBase(1);
	      int j_hi_corner =  xrM[k].getLocalBound(1);
	      for (int j = 0; j < j_low_corner; j++)
		for (int i = 0; i < dim1; i++)
		  xcogFile >> word;
	      for (int j1 = j_low_corner; j1 <= j_hi_corner; j1++) {
		for (int i = 0; i < i_low_corner; i++) 
		  xcogFile >> word;
		for (int i1 = i_low_corner; i1 <= i_hi_corner; i1++) 
		  xcogFile >> (*localXr)(i1,j1);
		for (int i2 = i_hi_corner + 1; i2 < dim1; i2++) 
		  xcogFile >> word;
	      }
	    }
	  } else if(!strcmp(word, "xs")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    xcogFile >> word;
	    xcogFile >> dim2;
	    for (lineSkip = 0; lineSkip < 8; lineSkip++)
	      xcogFile >> word;
	    
	    if (gridTypeM[k] != 1) {  // curvilinear 
	      xsM[k].partition(ddM[k]);
	      
	      xsM[k].redim(dim1, dim2);
	      
	      doubleSerialArray *localXs;
	      localXs = xsM[k].getSerialArrayPointer();
	      int i_low_corner = xsM[k].getLocalBase(0);
	      int i_hi_corner =  xsM[k].getLocalBound(0);
	      int j_low_corner = xsM[k].getLocalBase(1);
	      int j_hi_corner =  xsM[k].getLocalBound(1);
	      for (int j = 0; j < j_low_corner; j++)
		for (int i = 0; i < dim1; i++)
		  xcogFile >> word;
	      for (int j1 = j_low_corner; j1 <= j_hi_corner; j1++) {
		for (int i = 0; i < i_low_corner; i++) 
		  xcogFile >> word;
		for (int i1 = i_low_corner; i1 <= i_hi_corner; i1++) 
		  xcogFile >> (*localXs)(i1, j1);
		for (int i2 = i_hi_corner + 1; i2 < dim1; i2++) 
		  xcogFile >> word;
	      }
	    }
	  } else if(!strcmp(word, "yr")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    xcogFile >> word;
	    xcogFile >> dim2;
	    for (lineSkip = 0; lineSkip < 8; lineSkip++)
	      xcogFile >> word;
	    
	    if (gridTypeM[k] != 1) {  // curvilinear
			  
	      yrM[k].partition(ddM[k]);
	      
	      yrM[k].redim(dim1, dim2);
	      
	      doubleSerialArray *localYr;
	      localYr = yrM[k].getSerialArrayPointer();
	      int i_low_corner = yrM[k].getLocalBase(0);
	      int i_hi_corner =  yrM[k].getLocalBound(0);
	      int j_low_corner = yrM[k].getLocalBase(1);
	      int j_hi_corner =  yrM[k].getLocalBound(1);
	      for (int j = 0; j < j_low_corner; j++)
		for (int i = 0; i < dim1; i++)
		  xcogFile >> word;
	      for (int j1 = j_low_corner; j1 <= j_hi_corner; j1++) {
		for (int i = 0; i < i_low_corner; i++) 
		  xcogFile >> word;
		for (int i1 = i_low_corner; i1 <= i_hi_corner; i1++) 
		  xcogFile >> (*localYr)(i1,j1);
		for (int i2 = i_hi_corner+1; i2 < dim1; i2++) 
		  xcogFile >> word;
	      }
	    }
	  } else if (!strcmp(word, "ys")) {
	    xcogFile >> word;
	    xcogFile >> dim1;
	    xcogFile >> word;
	    xcogFile >> dim2;
	    for (lineSkip = 0; lineSkip < 8; lineSkip++)
	      xcogFile >> word;
	    
	    if (gridTypeM[k] != 1) {   //curvilinear grid 

	      ysM[k].partition(ddM[k]);
	      
	      ysM[k].redim(dim1, dim2);
	      
	      doubleSerialArray *localYs;
	      localYs = ysM[k].getSerialArrayPointer();
	      int i_low_corner = ysM[k].getLocalBase(0);
	      int i_hi_corner =  ysM[k].getLocalBound(0);
	      int j_low_corner = ysM[k].getLocalBase(1);
	      int j_hi_corner =  ysM[k].getLocalBound(1);
	      for (int j = 0; j < j_low_corner; j++)
		for (int i = 0; i < dim1; i++)
		  xcogFile >> word;
	      for (int j1 = j_low_corner; j1 <= j_hi_corner; j1++) {
		for (int i = 0; i < i_low_corner; i++) 
		  xcogFile >> word;
		for (int i1 = i_low_corner; i1 <= i_hi_corner; i1++) 
		  xcogFile >> (*localYs)(i1, j1);
		for (int i2 = i_hi_corner + 1; i2 < dim1; i2++) 
		  xcogFile >> word;
	      }
	    }
	  }
	  
	}
	xcogFile >> word;
      } //while
      
    } // for ... nmbrOfGridsM
      
    xcogFile.close();
  } else if (typ == HDF5) {
    hid_t       file_id, file_props, dataset_id, int_tid, char_tid;
    herr_t      status;
    
    file_props = H5Pcreate (H5P_FILE_ACCESS);
    status = H5Pset_fapl_mpio(file_props, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(status != FAIL);
    
    //----------------------------------------
    // Open an existing Xcog-file.
    //----------------------------------------
    file_id = H5Fopen(xcogFileName.c_str(), H5F_ACC_RDONLY, file_props);
    assert(file_id != FAIL);
    
    status = H5Pclose(file_props);
    assert(status != FAIL);
    //----------------------------------------
    // Open and read number of component grids,
    // a compound hdf5 data type
    //----------------------------------------
    dataset_id = H5Dopen(file_id, "/root/overlapping grid/n_components", H5P_DEFAULT);
    
    int_tid = H5Tcreate(H5T_COMPOUND, sizeof(int));
    status = H5Tinsert(int_tid, "n_components", 0, H5T_NATIVE_INT);
    assert(status != FAIL);
    status = H5Dread(dataset_id, int_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nmbrOfGridsM);
    assert(status != FAIL);
    status = H5Dclose(dataset_id);
    assert(status != FAIL);

    //----------------------------------------
    // Allocate the correct number of components
    // for all members of the composite grid
    //----------------------------------------
    nr_int_point = new int[nmbrOfGridsM];
    rdimM.resize(nmbrOfGridsM);
    sdimM.resize(nmbrOfGridsM);
    r_stepM = new double[nmbrOfGridsM];
    s_stepM = new double[nmbrOfGridsM];
    i_point = new intSerialArray[nmbrOfGridsM];
    j_point = new intSerialArray[nmbrOfGridsM];
    i_loc = new intSerialArray[nmbrOfGridsM];
    j_loc = new intSerialArray[nmbrOfGridsM];
    gridTypeM.resize(nmbrOfGridsM);
    gridLoc = new intSerialArray[nmbrOfGridsM];
    r_loc = new doubleSerialArray[nmbrOfGridsM];
    s_loc = new doubleSerialArray[nmbrOfGridsM];
    bcsM = new BC[nmbrOfGridsM];
    
    xM = new doubleArray[nmbrOfGridsM];
    yM = new doubleArray[nmbrOfGridsM];
    flagValuesM = new intArray[nmbrOfGridsM];
    maskM = new doubleArray[nmbrOfGridsM];
    xrM = new doubleArray[nmbrOfGridsM];
    xsM = new doubleArray[nmbrOfGridsM];
    yrM = new doubleArray[nmbrOfGridsM];
    ysM = new doubleArray[nmbrOfGridsM];
    xrrM = new doubleArray[nmbrOfGridsM];
    xssM = new doubleArray[nmbrOfGridsM];
    yrrM = new doubleArray[nmbrOfGridsM];
    yssM = new doubleArray[nmbrOfGridsM];
    sqrtOfGM = new doubleArray[nmbrOfGridsM];
    
    myNrOfIntPointsM.redim(nmbrOfGridsM, nrProcs);
    nrOfPointsToComputeM.redim(nmbrOfGridsM,nmbrOfGridsM,nrProcs);
    
    offsetsM.redim(nrProcs,nmbrOfGridsM,nmbrOfGridsM);
    
    receiveTypeM = new MPI_Datatype[nrProcs];
    sendTypeM = new MPI_Datatype[nrProcs];
    
    myIntPointsM = new intSerialArray[nmbrOfGridsM];
    My_i_LM = new intSerialArray[nmbrOfGridsM];
    My_j_LM = new intSerialArray[nmbrOfGridsM];
    intInterpolationLocationM = new intSerialArray[nmbrOfGridsM];
    doubleInterpolationLocationM = new doubleSerialArray[nmbrOfGridsM];
    interpolationCoordinatesM = new doubleSerialArray[nmbrOfGridsM];
    
    
    ddM = new Partitioning_Type[nmbrOfGridsM];
    if (gridDistribution != 0) {
      for (int dist=0; dist<nmbrOfGridsM; dist++) {
	ddM[dist] = gridDistribution[dist];
      }
    }      
    hid_t group_id, dset_id;
    size_t size;
    
    for (int k = nmbrOfGridsM; k >= 1; k--) {
      //----------------------------------------
      // Open the group for component grid [k-1]
      // and use as offset for all datasets 
      //----------------------------------------
      sprintf(word, "/root/overlapping grid/component grid %d",k);
      
      group_id = H5Gopen(file_id, word, H5P_DEFAULT);
      
      //----------------------------------------
      // Open and read type of this grid
      //----------------------------------------
      dset_id = H5Dopen(group_id, "grid type", H5P_DEFAULT);
      
      size = H5Dget_storage_size(dset_id);
      
      hid_t string_tid;
      string_tid = H5Tcopy(H5T_C_S1); 
      status = H5Tset_size(string_tid, size);
      assert(status != FAIL);
      status = H5Tset_strpad(string_tid, H5T_STR_SPACEPAD);
      assert(status != FAIL);
      status = H5Tset_cset(string_tid, H5T_CSET_ASCII );
      assert(status != FAIL);
      char_tid = H5Tcreate(H5T_COMPOUND, size);
      
      status = H5Tinsert(char_tid, "grid type", 0, string_tid);
      assert(status != FAIL);
      
      status = H5Dread(dset_id, char_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT,  word);
      assert(status != FAIL);
      
      if (!strcmp(word,"c         a         r         t         e         s         i         a         n         "))
	gridTypeM[k-1] = 1;
      else if (!strcmp(word,"c           u           r           v           i           l           i           n           e           a           r           "))
	gridTypeM[k-1] = 3;
      else
	gridTypeM[k-1] = 3;
      
      /* Close the dataset. */
      status = H5Dclose(dset_id);
      assert(status != FAIL);
      
      //----------------------------------------
      // Open dataset and read the number of 
      // interpolation points in this component
      // grid
      //----------------------------------------
      dset_id = H5Dopen(group_id, "n_interp", H5P_DEFAULT);
      int_tid = H5Tcreate(H5T_COMPOUND, sizeof(int));
      status = H5Tinsert(int_tid, "n_interp", 0, H5T_NATIVE_INT);
      assert(status != FAIL);
      status = H5Dread(dset_id, int_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(nr_int_point[k-1]));
      assert(status != FAIL);
      status = H5Dclose(dset_id);
      assert(status != FAIL);
      
      //----------------------------------------
      // If the number of interpolation points >0
      // , read rest of interpolation information
      //----------------------------------------
      if (nr_int_point[k-1] > 0) {
	Index even(0, nr_int_point[k-1], 2), odd(1, nr_int_point[k-1], 2);
	Index tri0(0, nr_int_point[k-1], 3), tri1(1, nr_int_point[k-1], 3), tri2(2, nr_int_point[k-1], 3);
	
	i_point[k-1].redim(nr_int_point[k-1]);
	j_point[k-1].redim(nr_int_point[k-1]);
	
	i_loc[k-1].redim(nr_int_point[k-1]);
	j_loc[k-1].redim(nr_int_point[k-1]);
	gridLoc[k-1].redim(nr_int_point[k-1]);
	
	r_loc[k-1].redim(nr_int_point[k-1]);
	s_loc[k-1].redim(nr_int_point[k-1]);
	
	intSerialArray readintArray(3*nr_int_point[k-1]);
	
	dset_id = H5Dopen(group_id, "interpolation point", H5P_DEFAULT);
	status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, readintArray.getDataPointer());
	assert(status != FAIL);
	status = H5Dclose(dset_id);
	assert(status != FAIL);
	i_point[k-1] = readintArray(even);
	j_point[k-1] = readintArray(odd);
	
	dset_id = H5Dopen(group_id, "donor point", H5P_DEFAULT);
	status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, readintArray.getDataPointer());
	assert(status != FAIL);
	status = H5Dclose(dset_id);
	assert(status != FAIL);
	i_loc[k-1] = readintArray(tri0);
	j_loc[k-1] = readintArray(tri1);
	gridLoc[k-1] = readintArray(tri2);
	
	doubleSerialArray readdoubleArray(2*nr_int_point[k-1]);
	dset_id = H5Dopen(group_id, "donor parameter", H5P_DEFAULT);
	status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, readdoubleArray.getDataPointer());
	assert(status != FAIL);
	status = H5Dclose(dset_id);
	assert(status != FAIL);
	r_loc[k-1] = readdoubleArray(even);
	s_loc[k-1] = readdoubleArray(odd);
      }
      
      //----------------------------------------
      // Read dimensions, gridsizes, and 
      // boundary conditions
      //----------------------------------------
      int dim[2];
      dset_id = H5Dopen(group_id, "dimension", H5P_DEFAULT);
      status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dim);
      status = H5Dclose(dset_id);
      assert(status != FAIL);
      rdimM[k-1] = dim[0]; sdimM[k-1] = dim[1];
      
      double step[2];
      dset_id = H5Dopen(group_id, "step", H5P_DEFAULT);
      status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &step);
      assert(status != FAIL);
      status = H5Dclose(dset_id);
      assert(status != FAIL);
      r_stepM[k-1] = step[0]; s_stepM[k-1] = step[1];
      
      int bc[4];
      dset_id = H5Dopen(group_id, "boundary condition", H5P_DEFAULT);
      status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bc);
      assert(status != FAIL);
      status = H5Dclose(dset_id);
      assert(status != FAIL);
      bcsM[k-1].bP[0] = bc[0]; bcsM[k-1].bP[1] = bc[1];
      bcsM[k-1].bP[2] = bc[2]; bcsM[k-1].bP[3] = bc[3];
      
      //----------------------------------------
      // Setup information needed to do
      // collective I/O, and store data
      // in P++-arrays
      //----------------------------------------
      
      hid_t xfer_plist;         
      xfer_plist = H5Pcreate (H5P_DATASET_XFER);
      assert(xfer_plist != FAIL);
      status = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
      assert(status != FAIL);
      
      hid_t file_dataspace, mem_dataspace;
      hsize_t count[2], in_count[2];
      hsize_t fstart[2], sstart[2];
      
      flagValuesM[k-1].partition(ddM[k-1]);
      flagValuesM[k-1].redim(rdimM[k-1], sdimM[k-1]);
      
      //----------------------------------------
      // flagValuesM is used as a template 
      // for distribution of data over processes
      //----------------------------------------
      int bottom0, bottom1, top0, top1;
      
      bottom0 = ( flagValuesM[k-1].getLocalBase(0) == flagValuesM[k-1].getBase(0) )?0:1;
      bottom1 = ( flagValuesM[k-1].getLocalBase(1) == flagValuesM[k-1].getBase(1) )?0:1;
      top0    = ( flagValuesM[k-1].getLocalBound(0) == flagValuesM[k-1].getBound(0) )?0:1;
      top1    = ( flagValuesM[k-1].getLocalBound(1) == flagValuesM[k-1].getBound(1) )?0:1;
      
      sstart[1] = flagValuesM[k-1].getLocalBase(0);
      sstart[0] = flagValuesM[k-1].getLocalBase(1);
      count[1] = flagValuesM[k-1].getLocalBound(0) - sstart[1] + 1;
      count[0] = flagValuesM[k-1].getLocalBound(1) - sstart[0] + 1;
      
      in_count[0] = flagValuesM[k-1].getLocalLength(1) - bottom1 - top1;
      in_count[1] = flagValuesM[k-1].getLocalLength(0) - bottom0 - top0;
      
      //----------------------------------------
      // Open and read flag dataset
      //----------------------------------------
      dset_id = H5Dopen(group_id, "flag", H5P_DEFAULT);
      
      file_dataspace = H5Dget_space (dset_id);
      assert(file_dataspace != FAIL);
      
      sstart[0] += bottom1; sstart[1] += bottom0;
      status=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, sstart, NULL,
				 in_count, NULL);
      assert(status != FAIL);
      
      if (bottom0 == 0) count[1]++; if (bottom1 == 0) count[0]++;
      if (top0 == 0)    count[1]++; if (top1 == 0)    count[0]++;
      
      mem_dataspace = H5Screate_simple (2, count, NULL);
      assert(mem_dataspace != FAIL);
      
      fstart[0] = 1; fstart[1] = 1;
      status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, fstart, NULL, 
				   in_count, NULL);
      assert(status != FAIL);
      
      status = H5Dread(dset_id, H5T_NATIVE_INT, mem_dataspace, file_dataspace,
		       xfer_plist, flagValuesM[k-1].getDataPointer());
      assert(status != FAIL);
      
      status = H5Sclose(mem_dataspace);
      status = H5Sclose(file_dataspace);
      status = H5Dclose(dset_id);
      
      xM[k - 1].partition(ddM[k - 1]);
      xM[k - 1].redim(rdimM[k - 1], sdimM[k - 1]);
      
      //----------------------------------------
      // Open and read x dataset
      //----------------------------------------
      dset_id = H5Dopen(group_id, "x", H5P_DEFAULT);
      
      file_dataspace = H5Dget_space (dset_id);                               
      assert(file_dataspace != FAIL);
      
      status=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, sstart, NULL,
				 in_count, NULL);
      assert(status != FAIL);
      
      mem_dataspace = H5Screate_simple (2, count, NULL);
      assert(mem_dataspace != FAIL);
      
      status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, fstart, NULL, 
				   in_count, NULL);
      assert(status != FAIL);
      
      status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace,
		       xfer_plist, xM[k-1].getDataPointer());
      assert(status != FAIL);
      
      status = H5Sclose(mem_dataspace);
      status = H5Sclose(file_dataspace);
      status = H5Dclose(dset_id);
      
      yM[k-1].partition(ddM[k-1]);
      yM[k-1].redim(rdimM[k-1], sdimM[k-1]);
      
      //----------------------------------------
      // Open and read y dataset
      //----------------------------------------
      dset_id = H5Dopen(group_id, "y", H5P_DEFAULT);
      
      file_dataspace = H5Dget_space (dset_id);                               
      assert(file_dataspace != FAIL);
      
      status=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, sstart, NULL,
				 in_count, NULL);
      assert(status != FAIL);
      
      mem_dataspace = H5Screate_simple (2, count, NULL);
      assert(mem_dataspace != FAIL);
      
      status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, fstart, NULL, 
				   in_count, NULL);
      assert(status != FAIL);
      
      status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace,
		       xfer_plist, yM[k-1].getDataPointer());
      assert(status != FAIL);
      
      status = H5Sclose(mem_dataspace);
      status = H5Sclose(file_dataspace);
      status = H5Dclose(dset_id);
      
      
      if (gridTypeM[k - 1] == 1) {
	//----------------------------------------
	// This means grid is at least Cartesian
	// , and maybe even regular
	//----------------------------------------
	//----------------------------------------
	// Control whether the grid really is regular
	// or not
	//----------------------------------------
	checkIfGridIsRegular(k - 1);
      } else {
	//----------------------------------------
	// Otherwise it is Curvilinear
	//----------------------------------------
	gridTypeM[k] = 3;
      }
      if (gridTypeM[k - 1] != 1) {
	xrM[k - 1].partition(ddM[k-1]);
	xrM[k - 1].redim(rdimM[k-1], sdimM[k-1]);
	
	//----------------------------------------
	// Open and read xr dataset
	//----------------------------------------
	dset_id = H5Dopen(group_id, "xr", H5P_DEFAULT);
	
	file_dataspace = H5Dget_space (dset_id);                               
	assert(file_dataspace != FAIL);
	
	status=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, sstart, NULL,
				   in_count, NULL);
	assert(status != FAIL);
	
	mem_dataspace = H5Screate_simple (2, count, NULL);
	assert(mem_dataspace != FAIL);
	
	status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, fstart, NULL, 
				     in_count, NULL);
	assert(status != FAIL);
	
	status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace,
			 xfer_plist, xrM[k-1].getDataPointer());
	assert(status != FAIL);
	
	status = H5Sclose(mem_dataspace);
	status = H5Sclose(file_dataspace);
	status = H5Dclose(dset_id);
	
	xsM[k-1].partition(ddM[k-1]);
	xsM[k-1].redim(rdimM[k-1], sdimM[k-1]);
	
	//----------------------------------------
	// Open and read xs dataset
	//----------------------------------------
	dset_id = H5Dopen(group_id, "xs", H5P_DEFAULT);
	
	file_dataspace = H5Dget_space (dset_id);                               
	assert(file_dataspace != FAIL);
	
	status=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, sstart, NULL,
				   in_count, NULL);
	assert(status != FAIL);
	
	mem_dataspace = H5Screate_simple (2, count, NULL);
	assert(mem_dataspace != FAIL);
	
	status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, fstart, NULL, 
				     in_count, NULL);
	assert(status != FAIL);
	
	status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace,
			 xfer_plist, xsM[k-1].getDataPointer());
	assert(status != FAIL);
	
	status = H5Sclose(mem_dataspace);
	status = H5Sclose(file_dataspace);
	status = H5Dclose(dset_id);
	
	yrM[k-1].partition(ddM[k-1]);
	yrM[k-1].redim(rdimM[k-1], sdimM[k-1]);
	
	//----------------------------------------
	// Open and read yr dataset
	//----------------------------------------
	dset_id = H5Dopen(group_id, "yr", H5P_DEFAULT);
	
	file_dataspace = H5Dget_space (dset_id);                               
	assert(file_dataspace != FAIL);
	
	status=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, sstart, NULL,
				   in_count, NULL);
	assert(status != FAIL);
	
	mem_dataspace = H5Screate_simple (2, count, NULL);
	assert(mem_dataspace != FAIL);
	
	status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, fstart, NULL, 
				     in_count, NULL);
	assert(status != FAIL);
	
	status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace,
			 xfer_plist, yrM[k-1].getDataPointer());
	assert(status != FAIL);
	
	status = H5Sclose(mem_dataspace);
	status = H5Sclose(file_dataspace);
	status = H5Dclose(dset_id);
	
	ysM[k-1].partition(ddM[k-1]);
	ysM[k-1].redim(rdimM[k-1], sdimM[k-1]);
	
	//----------------------------------------
	// Open and read ys dataset
	//----------------------------------------
	dset_id = H5Dopen(group_id, "ys", H5P_DEFAULT);
	
	file_dataspace = H5Dget_space (dset_id);                               
	assert(file_dataspace != FAIL);
	
	status=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, sstart, NULL,
				   in_count, NULL);
	assert(status != FAIL);
	
	mem_dataspace = H5Screate_simple (2, count, NULL);
	assert(mem_dataspace != FAIL);
	
	status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, fstart, NULL, 
				     in_count, NULL);
	assert(status != FAIL);
	
	status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace,
			 xfer_plist, ysM[k-1].getDataPointer());
	assert(status != FAIL);
	
	status = H5Sclose(mem_dataspace);
	status = H5Sclose(file_dataspace);
	status = H5Dclose(dset_id);
      }
      
      H5Pclose(xfer_plist);
      status = H5Gclose(group_id);
      
    }
    //----------------------------------------
    // Close the file
    //----------------------------------------
    status = H5Fclose(file_id);
  }
  
  int kk;
  
  for (int cuG=0;cuG<nmbrOfGridsM; cuG++) {
    maskM[cuG].partition(ddM[cuG]);
    maskM[cuG].redim(rdimM[cuG],sdimM[cuG]);
    intSerialArray *localFlags;
    localFlags = flagValuesM[cuG].getSerialArrayPointer();
    doubleSerialArray *localMask;
    localMask = maskM[cuG].getSerialArrayPointer();
    for (int jy=(*localFlags).getBase(1); jy<=(*localFlags).getBound(1); jy++) {
      for (int ix=(*localFlags).getBase(0); ix<=(*localFlags).getBound(0); ix++) {
	if ((*localFlags)(ix,jy) == 0)
	  (*localMask)(ix,jy) = 0;
	else
	  (*localMask)(ix,jy) = 1;
      }
    }
    if (gridTypeM[cuG] != 1) {
      /*
	Compute all second order metric derivatives
      */
      xM[cuG].updateGhostBoundaries();
      yM[cuG].updateGhostBoundaries();
      
      Range inner_x(1,rdimM[cuG]-2), inner_y(1,sdimM[cuG]-2);
      Index all;
      
      xrrM[cuG].partition(ddM[cuG]);
      xssM[cuG].partition(ddM[cuG]);
      yrrM[cuG].partition(ddM[cuG]);
      yssM[cuG].partition(ddM[cuG]);

      xrrM[cuG].redim(rdimM[cuG],sdimM[cuG]);
      xssM[cuG].redim(rdimM[cuG],sdimM[cuG]);
      yrrM[cuG].redim(rdimM[cuG],sdimM[cuG]);
      yssM[cuG].redim(rdimM[cuG],sdimM[cuG]);
      
      xrrM[cuG](inner_x,all) = (xM[cuG](inner_x-1,all) - 2.*xM[cuG](inner_x,all) + xM[cuG](inner_x+1,all)) / (pow(r_stepM[cuG],2.));
      xrrM[cuG](0,all) = (2.*xM[cuG](0,all) - 5.*xM[cuG](1,all) + 4.*xM[cuG](2,all) - xM[cuG](3,all)) / (pow(r_stepM[cuG],2.));
      xrrM[cuG](rdimM[cuG]-1,all) = (2.*xM[cuG](rdimM[cuG]-1,all) - 5.*xM[cuG](rdimM[cuG]-2,all) + 4.*xM[cuG](rdimM[cuG]-3,all) - xM[cuG](rdimM[cuG]-4,all)) / (pow(r_stepM[cuG],2.));
      
      xssM[cuG](all,inner_y) = (xM[cuG](all,inner_y-1) - 2.*xM[cuG](all,inner_y) + xM[cuG](all,inner_y+1)) / (pow(s_stepM[cuG],2.));
      xssM[cuG](all,0) = (2.*xM[cuG](all,0) - 5.*xM[cuG](all,1) + 4.*xM[cuG](all,2) - xM[cuG](all,3)) / (pow(s_stepM[cuG],2.));
      xssM[cuG](all,sdimM[cuG]-1) = (2.*xM[cuG](all,sdimM[cuG]-1) - 5.*xM[cuG](all,sdimM[cuG]-2) + 4.*xM[cuG](all,sdimM[cuG]-3) - xM[cuG](all,sdimM[cuG]-4)) / (pow(s_stepM[cuG],2.));
      
      yrrM[cuG](inner_x,all) = (yM[cuG](inner_x-1,all) - 2.*yM[cuG](inner_x,all) + yM[cuG](inner_x+1,all)) / (pow(r_stepM[cuG],2.));
      yrrM[cuG](0,all) = (2.*yM[cuG](0,all) - 5.*yM[cuG](1,all) + 4.*yM[cuG](2,all) - yM[cuG](3,all)) / (pow(r_stepM[cuG],2.));
      yrrM[cuG](rdimM[cuG]-1,all) = (2.*yM[cuG](rdimM[cuG]-1,all) - 5.*yM[cuG](rdimM[cuG]-2,all) + 4.*yM[cuG](rdimM[cuG]-3,all) - yM[cuG](rdimM[cuG]-4,all)) / (pow(r_stepM[cuG],2.));
      
      yssM[cuG](all,inner_y) = (yM[cuG](all,inner_y-1) - 2.*yM[cuG](all,inner_y) + yM[cuG](all,inner_y+1)) / (pow(s_stepM[cuG],2.));
      yssM[cuG](all,0) = (2.*yM[cuG](all,0) - 5.*yM[cuG](all,1) + 4.*yM[cuG](all,2) - yM[cuG](all,3)) / (pow(s_stepM[cuG],2.));
      yssM[cuG](all,sdimM[cuG]-1) = (2.*yM[cuG](all,sdimM[cuG]-1) - 5.*yM[cuG](all,sdimM[cuG]-2) + 4.*yM[cuG](all,sdimM[cuG]-3) - yM[cuG](all,sdimM[cuG]-4)) / (pow(s_stepM[cuG],2.));
      
      sqrtOfGM[cuG].partition(ddM[cuG]);
      
      sqrtOfGM[cuG].redim(rdimM[cuG],sdimM[cuG]);
      sqrtOfGM[cuG] = xrM[cuG]*ysM[cuG] - xsM[cuG]*yrM[cuG];
    }
  }
  
  /*
    correct indices read from xcog to suit normal c-style i.e. 0...
  */
  for (kk=0; kk<nmbrOfGridsM; kk++){
    i_point[kk] -= 1;
    j_point[kk] -= 1;
    i_loc[kk] -= 1;
    j_loc[kk] -= 1;
    gridLoc[kk] -= 1;
  }
  
  /*
    Setup interpolation information for all
    processors for all grids
  */
  /*
    First sort the interpolation points in ascending grid order, this
    simplifies things
  */
  gridSort(i_point, j_point, i_loc, j_loc, gridLoc, r_loc, s_loc);
  
  getLocalInterp(nr_int_point, i_point, j_point, i_loc, j_loc, gridLoc, r_loc, s_loc);
  
  setBoundaryType();

  delete[] nr_int_point;
  delete[] i_point;
  delete[] j_point;
  delete[] i_loc;
  delete[] j_loc;
  delete[] gridLoc;
  delete[] r_loc;
  delete[] s_loc;
  
  lowerbound_0M.resize(nmbrOfGridsM);
  lowerbound_1M.resize(nmbrOfGridsM);
  upperbound_0M.resize(nmbrOfGridsM);
  upperbound_1M.resize(nmbrOfGridsM);
  for (kk=0; kk<nmbrOfGridsM; kk++) {
    intSerialArray *Flag = flagValuesM[kk].getSerialArrayPointer();
    
    if (Flag->getBase(0) == flagValuesM[kk].getBase(0))
      lowerbound_0M[kk] = Flag->getBase(0);
    else
      lowerbound_0M[kk] = Flag->getBase(0) + 1;
    
    if (Flag->getBase(1) == flagValuesM[kk].getBase(1))
      lowerbound_1M[kk] = Flag->getBase(1);
    else
      lowerbound_1M[kk] = Flag->getBase(1) + 1;
    
    if (Flag->getBound(0) == flagValuesM[kk].getBound(0))
      upperbound_0M[kk] = Flag->getBound(0);
    else
      upperbound_0M[kk] = Flag->getBound(0) - 1;
    
    if (Flag->getBound(1) == flagValuesM[kk].getBound(1))
      upperbound_1M[kk] = Flag->getBound(1);
    else
      upperbound_1M[kk] = Flag->getBound(1) - 1;
  }
}

CompositeGrid::CompositeGrid(int dim1, int dim2, Partitioning_Type gridDistribution[]) {
  fromFileM = false;
  int *nr_int_point;
  
  boundaryWidthM = 2;
  
  nmbrOfGridsM = 1;

  ddM = new Partitioning_Type[nmbrOfGridsM];
  if (gridDistribution != 0) {
    for (int dist=0; dist<nmbrOfGridsM; dist++) {
      ddM[dist] = gridDistribution[dist];
    }
  }
  
  nr_int_point = new int[nmbrOfGridsM];
  myIntPointsM = new intSerialArray[nmbrOfGridsM];
  
  xM = new doubleArray[nmbrOfGridsM];
  yM = new doubleArray[nmbrOfGridsM];
  flagValuesM = new intArray[nmbrOfGridsM];
  maskM = new doubleArray[nmbrOfGridsM];
  xrM = new doubleArray[nmbrOfGridsM];
  xsM = new doubleArray[nmbrOfGridsM];
  yrM = new doubleArray[nmbrOfGridsM];
  ysM = new doubleArray[nmbrOfGridsM];
  xrrM = new doubleArray[nmbrOfGridsM];
  xssM = new doubleArray[nmbrOfGridsM];
  yrrM = new doubleArray[nmbrOfGridsM];
  yssM = new doubleArray[nmbrOfGridsM];
  sqrtOfGM = new doubleArray[nmbrOfGridsM];
  
  r_stepM = new double[nmbrOfGridsM];
  s_stepM = new double[nmbrOfGridsM];
  gridTypeM.resize(nmbrOfGridsM);
  rdimM.resize(nmbrOfGridsM);
  sdimM.resize(nmbrOfGridsM);
  bcsM = new BC[nmbrOfGridsM];
  
  
  for (int k = nmbrOfGridsM-1; k >= 0; k--) {
    gridTypeM[k] = 1;
    nr_int_point[k] = 0;
    rdimM[k] = dim1+2;
    sdimM[k] = dim2+2;
    r_stepM[k] = 1./(dim1-1);
    s_stepM[k] = 1./(dim2-1);
    
    bcsM[k].bP[0] = 1;
    bcsM[k].bP[1] = 1;
    bcsM[k].bP[2] = 1;
    bcsM[k].bP[3] = 1;
    
    Index all;
    flagValuesM[k].partition(ddM[k]);
    flagValuesM[k].redim(rdimM[k],sdimM[k]);
    
    flagValuesM[k] = 1;
    flagValuesM[k](0,all) = 0;
    flagValuesM[k](rdimM[k]-1,all) = 0;
    flagValuesM[k](all,0) = 0;
    flagValuesM[k](all,sdimM[k]-1) = 0;
    
    xM[k].partition(ddM[k]);
    xM[k].redim(rdimM[k],sdimM[k]);
    
    for (int i=0; i<rdimM[k]; i++) {
      xM[k](i,all) = (i-1)*r_stepM[k];
    }
    
    yM[k].partition(ddM[k]);
    yM[k].redim(rdimM[k],sdimM[k]);
    
    for (int j=0; j<sdimM[k]; j++) {
      yM[k](all,j) = (j-1)*s_stepM[k];
    }	
    
  } // for ... nmbrOfGridsM
  
  for (int cuG=0;cuG<nmbrOfGridsM; cuG++) {
    maskM[cuG].partition(ddM[cuG]);
    maskM[cuG].redim(rdimM[cuG],sdimM[cuG]);
    
    maskM[cuG] = flagValuesM[cuG].convertTo_doubleArray();
  }
  
  setBoundaryType();
  
  lowerbound_0M.resize(nmbrOfGridsM);
  lowerbound_1M.resize(nmbrOfGridsM);
  upperbound_0M.resize(nmbrOfGridsM);
  upperbound_1M.resize(nmbrOfGridsM);
  for (int kk = 0; kk < nmbrOfGridsM; kk++) {
    intSerialArray *Flag = flagValuesM[kk].getSerialArrayPointer();
    
    if (Flag->getBase(0) == flagValuesM[kk].getBase(0))
      lowerbound_0M[kk] = Flag->getBase(0);
    else
      lowerbound_0M[kk] = Flag->getBase(0) + 1;
    
    if (Flag->getBase(1) == flagValuesM[kk].getBase(1))
      lowerbound_1M[kk] = Flag->getBase(1);
    else
      lowerbound_1M[kk] = Flag->getBase(1) + 1;
    
    if (Flag->getBound(0) == flagValuesM[kk].getBound(0))
      upperbound_0M[kk] = Flag->getBound(0);
    else
      upperbound_0M[kk] = Flag->getBound(0) - 1;
    
    if (Flag->getBound(1) == flagValuesM[kk].getBound(1))
      upperbound_1M[kk] = Flag->getBound(1);
    else
      upperbound_1M[kk] = Flag->getBound(1) - 1;
  }
}

CompositeGrid::~CompositeGrid() {
  delete []xM;
  delete []yM;
  delete []flagValuesM;
  delete []maskM;
  delete []xrM;
  delete []xsM;
  delete []yrM;
  delete []ysM;
  delete []xrrM;
  delete []xssM;
  delete []yrrM;
  delete []yssM;
  delete []sqrtOfGM;
  delete []myIntPointsM;
  delete []r_stepM;
  delete []s_stepM;
  delete []bcsM;
  if (fromFileM) {
    delete []receiveTypeM;
    delete []sendTypeM;
    delete []My_i_LM;
    delete []My_j_LM;
    delete []intInterpolationLocationM;
    delete []doubleInterpolationLocationM;
    delete []interpolationCoordinatesM;
  }
  delete[] ddM;
}

void CompositeGrid::setBoundaryType() {
  for (int k=0; k<nrGrids(); k++) {
    if (bcsM[k].bP[0] == 0) {  // interpolation
      bcsM[k].lowR = INTERPOLATION;
    } else if (bcsM[k].bP[0] == 1) {  // no-slip or inflow
      bcsM[k].lowR = NOSLIP;
    } else if (bcsM[k].bP[0] == 3) {  // periodic
      bcsM[k].lowR = PERIODIC;
    } else if (bcsM[k].bP[0] == 5) {  // slip wall
      bcsM[k].lowR = SLIP;
    } else if (bcsM[k].bP[0] == 20) {  // inflow
      bcsM[k].lowR = INFLOW;
    } else if (bcsM[k].bP[0] == 21) {  // outflow
      bcsM[k].lowR = OUTFLOW;
    } else if (bcsM[k].bP[0] == 22) {  // outflow with neumann pressure
      bcsM[k].lowR = OUTFLOW_NEUM_P;
    } else if (bcsM[k].bP[0] == 50) {  // twilight-zone
      bcsM[k].lowR = DIRICHLET;
    }
    if (bcsM[k].bP[1] == 0) {  // interpolation
      bcsM[k].hiR = INTERPOLATION;
    } else if (bcsM[k].bP[1] == 1) {  // no-slip or inflow
      bcsM[k].hiR = NOSLIP;
    } else if (bcsM[k].bP[1] == 3) {  // periodic
      bcsM[k].hiR = PERIODIC;
    } else if (bcsM[k].bP[1] == 5) {  // slip wall
      bcsM[k].hiR = SLIP;
    } else if (bcsM[k].bP[1] == 20) {  // inflow
      bcsM[k].hiR = INFLOW;
    } else if (bcsM[k].bP[1] == 21) {  // outflow
      bcsM[k].hiR = OUTFLOW;
    } else if (bcsM[k].bP[1] == 22) {  // outflow with neumann pressure
      bcsM[k].hiR = OUTFLOW_NEUM_P;
    } else if (bcsM[k].bP[1] == 50) {  // twilight-zone
      bcsM[k].hiR = DIRICHLET;
    }
    if (bcsM[k].bP[2] == 0) {  // interpolation
      bcsM[k].lowS = INTERPOLATION;
    } else if (bcsM[k].bP[2] == 1) {  // no-slip or inflow
      bcsM[k].lowS = NOSLIP;
    } else if (bcsM[k].bP[2] == 3) {  // periodic
      bcsM[k].lowS = PERIODIC;
    } else if (bcsM[k].bP[2] == 5) {  // slip wall
      bcsM[k].lowS = SLIP;
    } else if (bcsM[k].bP[2] == 20) {  // outflow
      bcsM[k].lowS = INFLOW;
    } else if (bcsM[k].bP[2] == 21) {  // outflow
      bcsM[k].lowS = OUTFLOW;
    } else if (bcsM[k].bP[2] == 22) {  // outflow with neumann pressure
      bcsM[k].lowS = OUTFLOW_NEUM_P;
    } else if (bcsM[k].bP[2] == 50) {  // twilight-zone
      bcsM[k].lowS = DIRICHLET;
    }
    if (bcsM[k].bP[3] == 0) {  // interpolation
      bcsM[k].hiS = INTERPOLATION;
    } else if (bcsM[k].bP[3] == 1) {  // no-slip or inflow
      bcsM[k].hiS = NOSLIP;
    } else if (bcsM[k].bP[3] == 3) {  // periodic
      bcsM[k].hiS = PERIODIC;
    } else if (bcsM[k].bP[3] == 5) {  // slip wall
      bcsM[k].hiS = SLIP;
    } else if (bcsM[k].bP[3] == 20) {  // outflow
      bcsM[k].hiS = INFLOW;
    } else if (bcsM[k].bP[3] == 21) {  // outflow
      bcsM[k].hiS = OUTFLOW;
    } else if (bcsM[k].bP[3] == 22) {  // outflow with neumann pressure
      bcsM[k].hiS = OUTFLOW_NEUM_P;
    } else if (bcsM[k].bP[3] == 50) {  // twilight-zone
      bcsM[k].hiS = DIRICHLET;
    }
  }
}

boundaryType CompositeGrid::getBoundaryType(int k, Side side) const {
  if (side == lowR)
    return bcsM[k].lowR;
  else if (side == hiR)
    return bcsM[k].hiR;
  else if (side == lowS)
    return bcsM[k].lowS;
  else if (side == hiS)
    return bcsM[k].hiS;
  else {
    std::cerr << "INCORRECT SIDE TYPE IN CompositeGrid::getBoundaryType(int k, Side side) " 
	      << endl;
    return (boundaryType) 0;
  }
}

void CompositeGrid::setGhostCellWidth() {
  for (int k=0; k<nmbrOfGridsM; k++) {
    ddM[k].SpecifyInternalGhostBoundaryWidths(1,1);
  }
}

void CompositeGrid::saveCoordinatesToFile() const {
  ofstream xC, yC, maskFile;

  int myid = Communication_Manager::My_Process_Number;
  int grid;
  char filename[120];

  for (grid=0; grid < nmbrOfGridsM; grid++) {
    sprintf(filename, "%s%d%d", "./matlabData/xCoordinates", myid, grid);
    xC.open(filename, ios::out);
    sprintf(filename, "%s%d%d", "./matlabData/yCoordinates", myid, grid);
    yC.open(filename, ios::out);
    sprintf(filename, "%s%d%d", "./matlabData/gridMask", myid, grid);
    maskFile.open(filename, ios::out);
    doubleSerialArray *localX = xM[grid].getSerialArrayPointer();
    doubleSerialArray *localY = yM[grid].getSerialArrayPointer();
    doubleSerialArray *localMask = maskM[grid].getSerialArrayPointer();
    
    for (int j=lb_1(grid); j<=ub_1(grid); j++) {
      for (int i=lb_0(grid); i<=ub_0(grid); i++) {
	xC.write((const char*) &(*localX)(i, j), sizeof(double));
	yC.write((const char*) &(*localY)(i, j), sizeof(double));
	maskFile.write((const char*) &(*localMask)(i, j), sizeof(double));
      }
    }
    
    xC.close();
    yC.close();
    maskFile.close();
  }
  
  sprintf(filename, "%s%d%s", "./matlabData/readx", myid, ".m");
  xC.open(filename, ios::out);
  sprintf(filename, "%s%d%s", "./matlabData/ready", myid, ".m");
  yC.open(filename, ios::out);
  sprintf(filename, "%s%d%s", "./matlabData/readmask", myid, ".m");
  maskFile.open(filename, ios::out);
  for (grid=0; grid<nmbrOfGridsM; grid++) {
    sprintf(filename, "%s%d%d", "./xCoordinates", myid, grid);
    xC << "fid = fopen('" << filename << "');" << endl;
    
    xC << "x" << myid << grid << "= fread(fid,[" << ub_0(grid)-lb_0(grid)+1 << " " << ub_1(grid)-lb_1(grid)+1 << "],'double');" << endl;
    xC << "st = fclose(fid);" << endl;
    
    sprintf(filename,"%s%d%d","./yCoordinates",myid,grid);
    yC << "fid = fopen('" << filename << "');" << endl;
    
    yC << "y" << myid << grid << "= fread(fid,[" << ub_0(grid)-lb_0(grid)+1 << " " << ub_1(grid)-lb_1(grid)+1 << "],'double');" << endl;
    yC << "st = fclose(fid);" << endl;
    
    sprintf(filename,"%s%d%d","./gridMask",myid,grid);
    maskFile << "fid = fopen('" << filename << "');" << endl;
    
    maskFile << "mask" << myid << grid << "= fread(fid,[" << ub_0(grid)-lb_0(grid)+1 << " " << ub_1(grid)-lb_1(grid)+1 << "],'double');" << endl;
    maskFile << "st = fclose(fid);" << endl;
  }
  xC.close();
  yC.close();
  maskFile.close();
  
  if (myid == 0) {
    xC.open("./matlabData/readx.m", ios::out);
    yC.open("./matlabData/ready.m", ios::out);
    maskFile.open("./matlabData/readmask.m", ios::out);
    for (int proc = 0; proc<Communication_Manager::numberOfProcessors(); proc++) {
      sprintf(filename, "%s%d", "readx", proc);
      xC << filename << endl;
      sprintf(filename, "%s%d", "ready", proc);
      yC << filename << endl;
      sprintf(filename, "%s%d", "readmask", proc);
      maskFile << filename << endl;
    }
  }
}

void CompositeGrid::saveCoordinatesToHDF5File(const char* name) const {
  hid_t       file_id, file_props, group_id;
  herr_t      status;
  
  file_props = H5Pcreate (H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(file_props, MPI_COMM_WORLD, MPI_INFO_NULL);
  assert(status != FAIL);
  
  //----------------------------------------
  // Open an existing hdf5-file 
  //----------------------------------------
  file_id = H5Fopen(name, H5F_ACC_RDWR, file_props);
  assert(file_id != FAIL);
  
  status = H5Pclose(file_props);
  assert(status != FAIL);  
  
  //----------------------------------------
  // Write coordinates for every component
  // grid
  //----------------------------------------
  char word[200];
  for (int k = 0; k < nmbrOfGridsM; k++) {
    //----------------------------------------
    // Open the group for component grid [k]
    // and use as offset for coordinates
    //----------------------------------------
    sprintf(word, "/root/overlapping grid/component grid %d",k);
    
    group_id = H5Gopen(file_id, word, H5P_DEFAULT);
    
    //----------------------------------------
    // Prepare for collective data transfer
    //----------------------------------------
    hid_t xfer_plist;         
    xfer_plist = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist != FAIL);
    status = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    assert(status != FAIL);
    
    hid_t file_dataspace, file_dataspace2, mem_dataspace, mem_dataspace2, sid1, sid2, dataset1, dataset2;
    hid_t lcpl = H5P_DEFAULT;
    hsize_t count[2], in_count[2], tot_count[2];
    hsize_t start[2];
    
    //----------------------------------------
    // x is used as a template 
    // for distribution of data over processes
    //----------------------------------------
    
    int bottom0, bottom1, top0, top1;
    
    bottom0 = ( flagValuesM[k].getLocalBase(0) == flagValuesM[k].getBase(0) )?0:1;
    bottom1 = ( flagValuesM[k].getLocalBase(1) == flagValuesM[k].getBase(1) )?0:1;
    top0    = ( flagValuesM[k].getLocalBound(0) == flagValuesM[k].getBound(0) )?0:1;
    top1    = ( flagValuesM[k].getLocalBound(1) == flagValuesM[k].getBound(1) )?0:1;
    
    start[1] = xM[k].getLocalBase(0);
    start[0] = xM[k].getLocalBase(1);
    count[1] = xM[k].getLocalBound(0) - start[1] + 1;
    count[0] = xM[k].getLocalBound(1) - start[0] + 1;
    
    in_count[0] = xM[k].getLocalLength(1) - bottom1 - top1;
    in_count[1] = xM[k].getLocalLength(0) - bottom0 - top0;
    
    tot_count[0] = xM[k].getLength(1);
    tot_count[1] = xM[k].getLength(0);
    //----------------------------------------
    // Create data [set,space] for x and y
    //----------------------------------------
    sid1 = H5Screate_simple (2, tot_count, NULL);
    assert(sid1 != FAIL);
    sid2 = H5Screate_simple (2, tot_count, NULL);
    assert(sid2 != FAIL);
    
    dataset1 = H5Dcreate(group_id, "x", H5T_NATIVE_DOUBLE, sid1, lcpl, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset1 != FAIL);
    dataset2 = H5Dcreate(group_id, "y", H5T_NATIVE_DOUBLE, sid1, lcpl, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset2 != FAIL);
    
    file_dataspace = H5Dget_space (dataset1);   
    assert(file_dataspace != FAIL);
    file_dataspace2 = H5Dget_space (dataset2);   
    assert(file_dataspace2 != FAIL);
    
    start[0] += bottom1; start[1] += bottom0;
    
    status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, NULL,
				 in_count, NULL);
    assert(status != FAIL);
    status = H5Sselect_hyperslab(file_dataspace2, H5S_SELECT_SET, start, NULL,
				 in_count, NULL);
    assert(status != FAIL);
    
    if (bottom0 == 0) count[1]++; if (bottom1 == 0) count[0]++;
    if (top0 == 0)    count[1]++; if (top1 == 0)    count[0]++;
    
    mem_dataspace = H5Screate_simple (2, count, NULL);
    assert(mem_dataspace != FAIL);
    mem_dataspace2 = H5Screate_simple (2, count, NULL);
    assert(mem_dataspace2 != FAIL);
    
    start[0] = 1; start[1] = 1;
    
    status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, start, NULL, 
				 in_count, NULL);
    assert(status != FAIL);
    status = H5Sselect_hyperslab(mem_dataspace2, H5S_SELECT_SET, start, NULL, 
				 in_count, NULL);
    assert(status != FAIL);
    //----------------------------------------
    // Write x and y to file
    //----------------------------------------
    status = H5Dwrite(dataset1, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace,
		      xfer_plist, xM[k].getDataPointer());
    assert(status != FAIL);
    
    status = H5Dwrite(dataset2, H5T_NATIVE_DOUBLE, mem_dataspace2, file_dataspace2,
		      xfer_plist, yM[k].getDataPointer());
    assert(status != FAIL);
    
    //----------------------------------------
    // Release resources
    //----------------------------------------
    H5Sclose(sid1);
    H5Dclose(dataset1);
    H5Sclose(file_dataspace);
    H5Sclose(mem_dataspace);
    H5Pclose(xfer_plist);
    H5Sclose(sid2);
    H5Dclose(dataset2);
    H5Sclose(file_dataspace2);
    H5Sclose(mem_dataspace2);
    
    //----------------------------------------
    // Close component grid group
    //----------------------------------------
    status = H5Gclose(group_id);
    assert(status != FAIL);
  }
  status = H5Fclose(file_id);
  assert(status != FAIL);
}

int CompositeGrid::nrGrids() const {
  return nmbrOfGridsM;
}

Range CompositeGrid::getBounds(int direction, int k) const {
  int low, hi;
  
  low = 0;

  if (direction == 0)
    hi = rdimM[k] - 1;
  else
    hi = sdimM[k] - 1;

  Range retRange(low,hi);
  
  return retRange;
}
  
void CompositeGrid::updateGhostBoundaries() {
  for (int k=0; k<nmbrOfGridsM; k++) {
    xM[k].updateGhostBoundaries();
    yM[k].updateGhostBoundaries();
    if (gridTypeM[k] != 1) {
      xM[k].updateGhostBoundaries();
      yM[k].updateGhostBoundaries();
      xrM[k].updateGhostBoundaries();
      xsM[k].updateGhostBoundaries();
      yrM[k].updateGhostBoundaries();
      ysM[k].updateGhostBoundaries();
      xrrM[k].updateGhostBoundaries();
      xssM[k].updateGhostBoundaries();
      yrrM[k].updateGhostBoundaries();
      yssM[k].updateGhostBoundaries();
      sqrtOfGM[k].updateGhostBoundaries();
    }
    
    flagValuesM[k].updateGhostBoundaries();
    maskM[k].updateGhostBoundaries();
  }
}

void CompositeGrid::getLocalInterp(int *nr_int_point,
				   intSerialArray *i_point,
				   intSerialArray *j_point,
				   intSerialArray *i_loc,
				   intSerialArray *j_loc,
				   intSerialArray *gridLoc,
				   doubleSerialArray *r_loc,
				   doubleSerialArray *s_loc) {
  int nrProcs, currGrid;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nrProcs);
  
  int i,j,wG,intProc,proc;
  int point[2];
  int low_i,low_j,hi_i,hi_j,count;
  
  myNrOfIntPointsM = 0;
  nrOfPointsToComputeM = 0;
  
  for (int k=0; k<nmbrOfGridsM; k++) {
    low_i = xM[k].getLocalBase(0) + xM[k].getInternalGhostCellWidth(0);
    hi_i = xM[k].getLocalBound(0) - xM[k].getInternalGhostCellWidth(0);
    low_j = xM[k].getLocalBase(1) + xM[k].getInternalGhostCellWidth(1);
    hi_j = xM[k].getLocalBound(1) - xM[k].getInternalGhostCellWidth(1);
    
    intSerialArray *localIndex = new intSerialArray[nrProcs];
    intSerialArray **nonlocalIndex = new intSerialArray*[nrProcs];
    for (i=0; i<nrProcs; i++) {
      localIndex[i].redim(nr_int_point[k]);localIndex[i] = 0;
      nonlocalIndex[i] = new intSerialArray[nmbrOfGridsM];
      for (j=0; j<nmbrOfGridsM; j++) {
	nonlocalIndex[i][j].redim(nr_int_point[j]);
	nonlocalIndex[i][j] = 0;
      }
    }
    
    int sumMy = 0, sumOther = 0;
    for (i=0; i<nr_int_point[k]; i++)
      if (nrProcs == 1) {
	localIndex[0]( myNrOfIntPointsM(k,0) ) = i;
	myNrOfIntPointsM(k,0)++;
	sumMy++;
      } else if ((i_point[k](i) >= low_i && i_point[k](i) <= hi_i) && (j_point[k](i) >= low_j && j_point[k](i) <= hi_j)) {
	point[0] = i_loc[k](i) + 1;   // add 1 to get center of interpolation stencil
	point[1] = j_loc[k](i) + 1;
	wG = gridLoc[k](i);
	intProc = xM[wG].Array_Descriptor.findProcNum(point);
	localIndex[intProc]( myNrOfIntPointsM(k,intProc) ) = i;
	myNrOfIntPointsM(k,intProc)++;
	sumMy++;
      }
      for (currGrid=0; currGrid<nmbrOfGridsM; currGrid++)
	for (i=0; i<nr_int_point[currGrid]; i++)
	  if (nrProcs == 1) {
	    if ( gridLoc[currGrid](i) == k) {
	      nonlocalIndex[0][currGrid]( nrOfPointsToComputeM(k,currGrid,0) ) = i;
	      nrOfPointsToComputeM(k,currGrid,0)++;
	      sumOther++;
	    }
	  } else if ((i_loc[currGrid](i)+1 >= low_i && i_loc[currGrid](i)+1 <= hi_i) && (j_loc[currGrid](i)+1 >= low_j && j_loc[currGrid](i)+1 <= hi_j) && gridLoc[currGrid](i) == k) {  // add 1 to [ij]_loc to get center of interpolation stencil
	    point[0] = i_point[currGrid](i);
	    point[1] = j_point[currGrid](i);
	    intProc = xM[currGrid].Array_Descriptor.findProcNum(point);
	    nonlocalIndex[intProc][currGrid]( nrOfPointsToComputeM(k,currGrid,intProc) ) = i;
	    nrOfPointsToComputeM(k,currGrid,intProc)++;
	    sumOther++;
	  }

      myIntPointsM[k].redim(sumMy,2);
      My_i_LM[k].redim(sumMy);
      My_j_LM[k].redim(sumMy);
      
      intInterpolationLocationM[k].redim(sumOther,2);
      doubleInterpolationLocationM[k].redim(sumOther,2);
      interpolationCoordinatesM[k].redim(sumOther,2);
      
      count = 0;
      for (proc=0; proc<nrProcs; proc++)
	for (i=0; i<myNrOfIntPointsM(k,proc); i++) {
	  myIntPointsM[k](count,I) = i_point[k](localIndex[proc](i));
	  myIntPointsM[k](count,J) = j_point[k](localIndex[proc](i));
	  My_i_LM[k](count) = i_loc[k](localIndex[proc](i));
	  My_j_LM[k](count) = j_loc[k](localIndex[proc](i));
	  count++;
	}
      
      count = 0;
      for (proc=0; proc<nrProcs; proc++)
	for (currGrid=0; currGrid<nmbrOfGridsM; currGrid++)
	  for (i=0; i<nrOfPointsToComputeM(k,currGrid,proc); i++) {
	    intInterpolationLocationM[k](count,I) = i_loc[currGrid](nonlocalIndex[proc][currGrid](i));
	    intInterpolationLocationM[k](count,J) = j_loc[currGrid](nonlocalIndex[proc][currGrid](i));
	    doubleInterpolationLocationM[k](count,I) = (double) i_loc[currGrid](nonlocalIndex[proc][currGrid](i));
	    doubleInterpolationLocationM[k](count,J) = (double) j_loc[currGrid](nonlocalIndex[proc][currGrid](i));
	    interpolationCoordinatesM[k](count,I) = r_loc[currGrid](nonlocalIndex[proc][currGrid](i));
	    interpolationCoordinatesM[k](count,J) = s_loc[currGrid](nonlocalIndex[proc][currGrid](i));
	    count++;
	  }
      
      for (int del=0;del<nrProcs;del++)
	delete[] nonlocalIndex[del];
      
      delete[] localIndex;
      delete[] nonlocalIndex;
    }
  
  int theSendBufferSize = 0, theReceiveBufferSize = 0;
  
  for (currGrid = 0; currGrid < nmbrOfGridsM; currGrid++) {
    theSendBufferSize += (intInterpolationLocationM[currGrid]).getLength(0);
    theReceiveBufferSize += (myIntPointsM[currGrid]).getLength(0);
  }
  
  theSendBufferM.redim(theSendBufferSize);
  theReceiveBufferM.redim(theReceiveBufferSize);
}

void CompositeGrid::printIntInf(int k) const {
  int nrProcs;
  
  Partitioning_Type local_DD = xM[k].getPartition();
  intSerialArray procSet = local_DD.getProcessorSet();
  nrProcs = procSet.getLength(0);
  
  int myid, low_p, allProcs;
  low_p = procSet.getBase(0);

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &allProcs);
  
  if (myid == procSet(low_p))
    std::cout << "\n===== Printing interpolation info. for component grid "<< k << " =====\n\n";
  
  for (int i=0; i<nrProcs; i++) {
    if (myid == procSet(low_p + i)) {
      std::cout << "Proc. nr " << myid << " has :\n";
      
      for (int p=0; p<allProcs; p++) {
	std::cout << " " << myNrOfIntPointsM(k,p) << " interpolation points in other grids, to be computed by proc. " << p << endl;
	for (int currGrid=0; currGrid<nmbrOfGridsM; currGrid++)
	  std::cout << " It computes interpolation info. for " << nrOfPointsToComputeM(k,currGrid,p) << " points in grid " << currGrid << " for proc. " << p << endl;
	std::cout << endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void CompositeGrid::gridSort(intSerialArray *i_point,
			     intSerialArray *j_point,
			     intSerialArray *i_loc,
			     intSerialArray *j_loc,
			     intSerialArray *gridLoc,
			     doubleSerialArray *r_loc,
			     doubleSerialArray *s_loc) {
  int i, nrOfSwaps;
  
  for (int cG = 0; cG < nmbrOfGridsM; cG++) {
    do {
      nrOfSwaps = 0;
      i = gridLoc[cG].getBase(0);
      while (i < gridLoc[cG].getBound(0)) {
	if (gridLoc[cG](i) > gridLoc[cG](i+1)) {
	  swap(gridLoc[cG](i),gridLoc[cG](i+1));
	  swap(i_point[cG](i),i_point[cG](i+1));
	  swap(j_point[cG](i),j_point[cG](i+1));
	  swap(i_loc[cG](i),i_loc[cG](i+1));
	  swap(j_loc[cG](i),j_loc[cG](i+1));
	  swap(r_loc[cG](i),r_loc[cG](i+1));
	  swap(s_loc[cG](i),s_loc[cG](i+1));
	  nrOfSwaps++;
	}
	i++;
      }
    } while (nrOfSwaps > 0);
  }
}

int CompositeGrid::getInterpIndex(int i, int j, int grid) const {
  int index;

  for (index = 0; index < myIntPointsM[grid].getLength(0); index++)
    if ((myIntPointsM[grid](index,I) == i) && (myIntPointsM[grid](index,J) == j))
      return index;

  std::cerr << "Error in Grid:: getInterpIndex i=" << i << ", j=" << j << ", grid=" << grid << "\n";
  return -1;
}

doubleArray CompositeGrid::normalVector_R_Component(int grid, int side, Index ixb, Index low_r, Index hi_r, Index iyb, Index low_s, Index hi_s) const {
  if (ixb.getBase() == -1000) 
    ixb = Index(1, rdimM[grid]-2);
  if (low_r.getBase() == -1000) 
    low_r = Index(1);
  if (hi_r.getBase() == -1000) 
    hi_r = Index(rdimM[grid]-2);
  if (iyb.getBase() == -1000) 
    iyb = Index(1, sdimM[grid]-2);
  if (low_s.getBase() == -1000) 
    low_s = Index(1);
  if (hi_s.getBase() == -1000) 
    hi_s = Index(sdimM[grid]-2); 

  if (gridTypeM[grid] == 1) {
    //----------------------------------------
    // Regular cartesian grid, easy
    //----------------------------------------
    if (side==0) {
      doubleArray one(low_r,iyb);
      one = 1;
      return one;
    } else if (side==1) {
      doubleArray one(hi_r,iyb);
      one = 1;
      return one;
    } else if (side==2) {
      doubleArray zero(ixb,low_s);
      zero = 0.;
      return zero;
    } else if (side==3) {
      doubleArray zero(ixb,hi_s);
      zero = 0.;
      return zero;
    } else {
      std::cerr << "INCORRECT SIDE IN CompositeGrid::normalVector_R_Component(int grid, int side, Index ixb, Index low_r, Index hi_r, Index iyb, Index low_s, Index hi_s)" 
		<< endl;
      return xM[grid];
    }
  } else {
    //----------------------------------------
    // Curvilinear grid, not so easy
    //----------------------------------------
    if (side==0)
      return ysM[grid](low_r,iyb) / pow(xsM[grid](low_r,iyb)*xsM[grid](low_r,iyb)+ysM[grid](low_r,iyb)*ysM[grid](low_r,iyb),0.5);
    else if (side==1)
      return ysM[grid](hi_r,iyb) / pow(xsM[grid](hi_r,iyb)*xsM[grid](hi_r,iyb)+ysM[grid](hi_r,iyb)*ysM[grid](hi_r,iyb),0.5) ;
    else if (side==2)
      return -yrM[grid](ixb,low_s) / pow(xrM[grid](ixb,low_s)*xrM[grid](ixb,low_s)+yrM[grid](ixb,low_s)*yrM[grid](ixb,low_s),0.5);
    else if (side==3)
      return -yrM[grid](ixb,hi_s) / pow(xrM[grid](ixb,hi_s)*xrM[grid](ixb,hi_s)+yrM[grid](ixb,hi_s)*yrM[grid](ixb,hi_s),0.5);
    else {
      std::cerr << "INCORRECT SIDE IN CompositeGrid::normalVector_R_Component(int grid, int side, Index ixb, Index low_r, Index hi_r, Index iyb, Index low_s, Index hi_s)" 
		<< endl;
      return xrM[grid];
    }
  }
}

doubleArray CompositeGrid::normalVector_S_Component(int grid, int side, Index ixb, Index low_r, Index hi_r, Index iyb, Index low_s, Index hi_s) const {
  if (ixb.getBase() == -1000) 
    ixb = Index(1, rdimM[grid]-2);
  if (low_r.getBase() == -1000) 
    low_r = Index(1);
  if (hi_r.getBase() == -1000) 
    hi_r = Index(rdimM[grid]-2);
  if (iyb.getBase() == -1000) 
    iyb = Index(1, sdimM[grid]-2);
  if (low_s.getBase() == -1000) 
    low_s = Index(1);
  if (hi_s.getBase() == -1000) 
    hi_s = Index(sdimM[grid]-2); 
  
  if (gridTypeM[grid] == 1) {
    //----------------------------------------
    // Regular cartesian grid, easy
    //----------------------------------------
    if (side==0) {
      doubleArray zero(low_r,iyb);
      zero = 0;
      return zero;
    } else if (side == 1) {
      doubleArray zero(hi_r,iyb);
      zero = 0;
      return zero;
    } else if (side == 2) {
      doubleArray one(ixb,low_s);
      one = 1;
      return one;
    } else if (side == 3) {
      doubleArray one(ixb,hi_s);
      one = 1;
      return one;
    } else {
      std::cerr << "INCORRECT SIDE IN CompositeGrid::normalVector_R_Component(int grid, int side, Index ixb, Index low_r, Index hi_r, Index iyb, Index low_s, Index hi_s)" 
		<< endl;
      return xM[grid];
    }
  } else {
    //----------------------------------------
    // Curvilinear grid, not so easy
    //----------------------------------------
    if (side==0)
      return  -xsM[grid](low_r,iyb) / pow(xsM[grid](low_r,iyb)*xsM[grid](low_r,iyb)+ysM[grid](low_r,iyb)*ysM[grid](low_r,iyb),0.5);
    else if (side == 1)
      return  -xsM[grid](hi_r,iyb) / pow(xsM[grid](hi_r,iyb)*xsM[grid](hi_r,iyb)+ysM[grid](hi_r,iyb)*ysM[grid](hi_r,iyb),0.5) ;
    else if (side == 2)
      return  xrM[grid](ixb,low_s) / pow(xrM[grid](ixb,low_s)*xrM[grid](ixb,low_s)+yrM[grid](ixb,low_s)*yrM[grid](ixb,low_s),0.5);
    else if (side == 3)
      return  xrM[grid](ixb,hi_s) / pow(xrM[grid](ixb,hi_s)*xrM[grid](ixb,hi_s)+yrM[grid](ixb,hi_s)*yrM[grid](ixb,hi_s),0.5);
    else {
      std::cerr << "INCORRECT SIDE IN CompositeGrid::normalVector_S_Component(int grid, int side, Index ixb, Index low_r, Index hi_r, Index iyb, Index low_s, Index hi_s)" 
		<< endl;
      return xrM[grid];
    }
  }
}

doubleArray CompositeGrid::gridSize(int grid, int direction1, int direction2, Index ix, Index iy) const {
  if (ix.getBase() == -1000) 
    ix = Index(1, rdimM[grid]-2);
  if (iy.getBase() == -1000) 
    iy = Index(1, sdimM[grid]-2);
  
  doubleArray result(ix,iy);
  
  if (direction2 == 0) {
    if (direction1 == 0)
      result = abs(xM[grid](ix+1,iy) - xM[grid](ix-1,iy))/2.;
    else if (direction1 == 1)
      result = abs(yM[grid](ix,iy+1) - yM[grid](ix,iy-1))/2.;
  } else if (direction2 == 1) {
    if (direction1 == 0)
      result = abs(xM[grid](ix,iy+1) - xM[grid](ix,iy-1))/2.;
    else if (direction1 == 1)
      result = abs(yM[grid](ix+1,iy) - yM[grid](ix-1,iy))/2.;
  }
  return result;
}

void CompositeGrid::checkIfGridIsRegular(int grid) {
  //----------------------------------------
  // Get bounds for local part of array
  //----------------------------------------

  int i_low_corner = yM[grid].getLocalBase(0);
  int i_hi_corner =  yM[grid].getLocalBound(0);
  int j_low_corner = yM[grid].getLocalBase(1);
  int j_hi_corner =  yM[grid].getLocalBound(1);

  doubleSerialArray *localX, *localY;

  localX = xM[grid].getSerialArrayPointer();
  localY = yM[grid].getSerialArrayPointer();

  Range ix(i_low_corner + 1, i_hi_corner - 1), iy(j_low_corner + 1, j_hi_corner - 1);

  int xrr = 0, xs = 0, yr = 0, yss = 0;
  double small = 0.001;

  //----------------------------------------
  // Test if xs == 0 everywhere
  //----------------------------------------

  if (max(abs(((*localX)(ix,iy+1) - (*localX)(ix,iy-1))/(2.*s_stepM[grid]))) > small)
    xs = 1;

  //----------------------------------------
  // Test if yr == 0 everywhere
  //----------------------------------------

  if (max(abs(((*localY)(ix+1,iy) - (*localY)(ix-1,iy))/(2.*r_stepM[grid]))) > small)
    yr = 1;

  //----------------------------------------
  // Test if xrr == 0 everywhere
  //----------------------------------------

  if (max(abs(((*localX)(ix+1,iy) - 2.*(*localX)(ix,iy) + (*localX)(ix-1,iy))/(r_stepM[grid]*r_stepM[grid]))) > small)
    xrr = 1;

  //----------------------------------------
  // Test if yss == everywhere
  //----------------------------------------

  if (max(abs(((*localY)(ix,iy+1) - 2.*(*localY)(ix,iy) + (*localY)(ix,iy-1))/(s_stepM[grid]*s_stepM[grid]))) > small) {
    yss = 1;
  }
  if ( (xs + xrr + yr + yss) == 0 ) {
    gridTypeM[grid] = 1;
  } else if ( (xs + yr) == 0 ) {
    gridTypeM[grid] = 2;
  } else {
    gridTypeM[grid] = 3;
  }
}

double CompositeGrid::hSquare() const {
  double result = 100000.;
  for (int cG = 0; cG < nmbrOfGridsM; cG++) {
    result = min(result, pow(r_stepM[cG],2.) + pow(s_stepM[cG],2.));
  }
  return result;
}

int CompositeGrid::rDim(int grid) const { return rdimM[grid]; };
int CompositeGrid::sDim(int grid) const { return sdimM[grid]; };
int CompositeGrid::gType(int grid) const { return gridTypeM[grid]; };
double CompositeGrid::rStep(int grid) const { return r_stepM[grid]; };
double CompositeGrid::sStep(int grid) const { return s_stepM[grid]; };

void CompositeGrid::simplifyFlagValues() {
  Index all;
  for (int cG = 0; cG < nmbrOfGridsM; cG++) {
    if (getBoundaryType(cG, lowR) != PERIODIC) {
      flagValuesM[cG](0, all) = 0;
      flagValuesM[cG](rdimM[cG] - 1, all) = 0;
    }
    if (getBoundaryType(cG, lowS) != PERIODIC) {
      flagValuesM[cG](all, 0) = 0;
      flagValuesM[cG](all, sdimM[cG] - 1) = 0;
    }
  }
}
