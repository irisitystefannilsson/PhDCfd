#include <mpi.h>


#include "GridFunc.hh"
#include <cmath>

#undef min
#undef max

#include <algorithm>

#include "hdf5.h"
#include <assert.h>
#define FAIL -1

doubleArray gridFunction::boundaryDiff(int i, int gType, const double t, const doubleArray & xIn, const doubleArray & yIn, const doubleArray & xr = doubleArray(), const doubleArray & xs = doubleArray(), const doubleArray & yr = doubleArray(), const doubleArray & ys = doubleArray())
{
  
  if (i == 0) // low r 
    {
      if (gType == 1)
	{
	  doubleArray resArray1 = - bc_x(xIn, yIn, t, 0, 0);
	  return resArray1;
	}
      else
	{
	  doubleArray resArray2 = ( ys*bc_x(xIn,yIn,t,0,0) - xs*bc_y(xIn,yIn,t,0,0) ) / pow(xs*xs+ys*ys,0.5);
	  return resArray2;
	}
    }
  else if (i == 1) // high r 
    {
      if (gType == 1)
	{
	  doubleArray resArray3 = bc_x(xIn,yIn,t,0,0);
	  return resArray3;
	}
      else
	{
	  doubleArray resArray4 = ( ys*bc_x(xIn,yIn,t,0,0) - xs*bc_y(xIn,yIn,t,0,0) ) / pow(xs*xs+ys*ys,0.5);
	  return resArray4;
	}
    }
  else if (i == 2) // low s 
    {
      if (gType == 1)
	{
	  doubleArray resArray5 = - bc_y(xIn,yIn,t,0,0);
	  return resArray5;
	}
      else
	{
	  doubleArray resArray6 = ( - yr*bc_x(xIn,yIn,t,0,0) + xr*bc_y(xIn,yIn,t,0,0) ) / pow(xr*xr+yr*yr,0.5);
	  return resArray6;
	}
    }
  else // high s 
    {
      if (gType == 1)
	{
	  doubleArray resArray7 = bc_y(xIn,yIn,t,0,0);
	  return resArray7;
	}
      else
	{
	  doubleArray resArray8 = ( - yr*bc_x(xIn,yIn,t,0,0) + xr*bc_y(xIn,yIn,t,0,0) ) / pow(xr*xr+yr*yr,0.5);
	  return resArray8;
	}
    }
}

gridFunction::gridFunction(CompositeGrid & Kimera, double present)
{
  myGridM = &Kimera;
  timeM = present;
  
  conservativeDifferenceOperatorsM = false;

  iTypeM = quadratic;

  fieldValuesM = new doubleArray[myGridM->nrGrids()];
  flagValuesM = new intArray[myGridM->nrGrids()];

  bcsM = new BC[myGridM->nrGrids()];

  boundaryValue = &GridFunction::zero;
  boundaryValue_d_dt = &GridFunction::zero;
  for (int k = 0; k < myGridM->nrGrids(); k++) 
    {
      int rdim, sdim;
      rdim = myGridM->rDim(k);
      sdim = myGridM->sDim(k);
      
      fieldValuesM[k].partition(myGridM->ddM[k]);
      flagValuesM[k].partition(myGridM->ddM[k]);
      
      fieldValuesM[k].redim(rdim, sdim);
      flagValuesM[k].redim(rdim, sdim);
      flagValuesM[k] = myGridM->flagValuesM[k];
      //        myGridM->flagValuesM[k].display();
      
    }
}

gridFunction::gridFunction(const gridFunction & fun)
{
  myGridM = fun.myGridM;
  timeM = fun.timeM;

  iTypeM = fun.iTypeM;

  fieldValuesM = new doubleArray[myGridM->nrGrids()];
  flagValuesM = new intArray[myGridM->nrGrids()];

  bcsM = new BC[myGridM->nrGrids()];

  for (int k = 0; k<myGridM->nrGrids(); k++) 
    {
      //    fieldValues[k].partition(fun.fieldValues[k]);
      //    flagValuesM[k].partition(fun.flagValuesM[k]);
      
      Index ix = fun.fieldValuesM[k].getFullRange(0);
      Index iy = fun.fieldValuesM[k].getFullRange(1);
      fieldValuesM[k].redim(ix,iy);
      
      //    ix = fun.flagValuesM[k].getFullRange(0);
      //    iy = fun.flagValuesM[k].getFullRange(1);
      flagValuesM[k].redim(ix,iy);
      
      fieldValuesM[k](ix,iy) = fun.fieldValuesM[k](ix,iy);
      flagValuesM[k](ix,iy) = fun.flagValuesM[k](ix,iy);
      
    }
}

gridFunction::gridFunction(int i)
{
}

gridFunction::~gridFunction()
{
  delete []fieldValuesM;
  delete []flagValuesM;
  delete []bcsM;
}

void gridFunction::printValues(char *info)
{
  int i,j,pxadd,pyadd;
  char filename[64];
  int myid = Communication_Manager::My_Process_Number;
  sprintf(filename,"gnuplotData/P++%s.%d",info,myid);

  FILE *fp;
  if ((fp = fopen(filename,"w")) == NULL)
    {
      std::cerr << "Could not open" << filename << endl;
    }
  for (int k=0; k<(myGridM->nmbrOfGridsM); k++)
    {
      pxadd = 0;
      pyadd = 0;
      
      if (bcsM[k].hiR == 3)
	pxadd = 1;
      if (bcsM[k].hiS == 3)
	pyadd = 1;
      doubleSerialArray *localXCoords = myGridM->xM[k].getSerialArrayPointer();
      doubleSerialArray *localYCoords = myGridM->yM[k].getSerialArrayPointer();
      doubleSerialArray *localField = fieldValuesM[k].getSerialArrayPointer();
      intSerialArray *localFlags = myGridM->flagValuesM[k].getSerialArrayPointer();
      
      for (i=(*localXCoords).getBase(0)+1;i<=(*localXCoords).getBound(0)-1+pxadd;i++) {
	for (j=(*localXCoords).getBase(1)+1;j<=(*localXCoords).getBound(1)-1+pyadd;j++) {
	  if ((*localFlags)(i,j))
	    fprintf(fp,"%f %f %f\n",(*localXCoords)(i,j),(*localYCoords)(i,j),(*localField)(i,j));
	}
      }
    }
  
  fclose(fp);
  
  sprintf(filename,"matlabData/PPP%s%d.m",info,myid);
  
  
  if ((fp = fopen(filename,"w")) == NULL) 
    {
      std::cerr << "Could not open" << filename << endl;
    }
  for (int cG=0; cG<(myGridM->nrGrids()); cG++)
    {
      pxadd = 0;
      pyadd = 0;
      
      if (bcsM[cG].hiR == PERIODIC)
	pxadd = 1;
      if (bcsM[cG].hiS == PERIODIC)
	pyadd = 1;
      doubleSerialArray *localXCoords = myGridM->xM[cG].getSerialArrayPointer();
      doubleSerialArray *localYCoords = myGridM->yM[cG].getSerialArrayPointer();
      doubleSerialArray *localField = fieldValuesM[cG].getSerialArrayPointer();
      intSerialArray *localFlags = myGridM->flagValuesM[cG].getSerialArrayPointer();
      
      fprintf(fp,"X%d%d=[",cG,myid);
      for (i=(*localXCoords).getBase(0)+1;i<=(*localXCoords).getBound(0)-1+pxadd;i++) {
	//    for (i=(*localXCoords).getBase(0);i<=(*localXCoords).getBound(0)+pxadd;i++) {
	fprintf(fp,"[");
	for (j=(*localXCoords).getBase(1)+1;j<=(*localXCoords).getBound(1)-1+pyadd;j++) {
	  //      for (j=(*localXCoords).getBase(1);j<=(*localXCoords).getBound(1)+pyadd;j++) {
	  fprintf(fp,"%f ",(*localXCoords)(i,j));
	}
	fprintf(fp,"]\n");
      }
      fprintf(fp,"];\n");
      
      fprintf(fp,"Y%d%d=[",cG,myid);
      for (i=(*localXCoords).getBase(0)+1;i<=(*localXCoords).getBound(0)-1+pxadd;i++) {
	//    for (i=(*localXCoords).getBase(0);i<=(*localXCoords).getBound(0)+pxadd;i++) {
	fprintf(fp,"[");
	for (j=(*localXCoords).getBase(1)+1;j<=(*localXCoords).getBound(1)-1+pyadd;j++) {
	  //      for (j=(*localXCoords).getBase(1);j<=(*localXCoords).getBound(1)+pyadd;j++) {
	  fprintf(fp,"%f ",(*localYCoords)(i,j));
	}
	fprintf(fp,"]\n");
      }
      fprintf(fp,"];\n");
      
      fprintf(fp,"%s%d%d=[",info,cG,myid);
      for (i=(*localXCoords).getBase(0)+1;i<=(*localXCoords).getBound(0)-1+pxadd;i++) {
	//    for (i=(*localXCoords).getBase(0);i<=(*localXCoords).getBound(0)+pxadd;i++) {
	fprintf(fp,"[");
	for (j=(*localXCoords).getBase(1)+1;j<=(*localXCoords).getBound(1)-1+pyadd;j++) {
	  //      for (j=(*localXCoords).getBase(1);j<=(*localXCoords).getBound(1)+pyadd;j++) {
	  fprintf(fp,"%f ",(*localField)(i,j));
	}
	fprintf(fp,"]\n");
      }
      fprintf(fp,"];\n");
      
      fprintf(fp,"Mask%d%d=[",cG,myid);
      for (i=(*localXCoords).getBase(0)+1;i<=(*localXCoords).getBound(0)-1+pxadd;i++) {
	//    for (i=(*localXCoords).getBase(0);i<=(*localXCoords).getBound(0)+pxadd;i++) {
	fprintf(fp,"[");
	for (j=(*localXCoords).getBase(1)+1;j<=(*localXCoords).getBound(1)-1+pyadd;j++) {
	  //      for (j=(*localXCoords).getBase(1);j<=(*localXCoords).getBound(1)+pyadd;j++) {
	  if ((*localFlags)(i,j) == 0)
	    fprintf(fp,"NaN ");
	  else
	    fprintf(fp,"%d ",1);
	}
	fprintf(fp,"]\n");
      }
      fprintf(fp,"];\n");
      
      fprintf(fp,"%s%d%d=%s%d%d.*Mask%d%d ;\n",info,cG,myid,info,cG,myid,cG,myid);
    }
  
  
  fclose(fp);
}

void gridFunction::saveToFile(const char *name)
{
  ofstream varfile;
  char filename[120];

  int myid = Communication_Manager::My_Process_Number;
  int grid;

  for (grid=0; grid<myGridM->nrGrids(); grid++)
    {
      sprintf(filename,"%s%s%d%d","./matlabData/",name,myid,grid);
      varfile.open(filename,ios::out);
      doubleSerialArray *local = fieldValuesM[grid].getSerialArrayPointer();
      intSerialArray *localFlags = flagValuesM[grid].getSerialArrayPointer();

      for (int j=myGridM->lb_1(grid); j<=myGridM->ub_1(grid); j++)
	{
	  for (int i=myGridM->lb_0(grid); i<=myGridM->ub_0(grid); i++)
	    {
	      if ((*localFlags)(i,j) == 0)
		{
		  double notanumber = HUGE_VAL;
		  varfile.write((const char*) &notanumber,sizeof(double));
		}
	      else
		varfile.write((const char*) &(*local)(i,j),sizeof(double));
	    }
	}

      varfile.close();
    }
  sprintf(filename,"%s%s%d%s","./matlabData/read",name,myid,".m");
  varfile.open(filename,ios::out);
  for (grid=0; grid<myGridM->nrGrids(); grid++)
    {
      sprintf(filename,"%s%s%d%d","./",name,myid,grid);
      varfile << "fid = fopen('" << filename << "');" << endl;
      
      varfile << name << myid << grid << "= fread(fid,[" << myGridM->ub_0(grid)-myGridM->lb_0(grid)+1 << " " << myGridM->ub_1(grid)-myGridM->lb_1(grid)+1 << "],'double');" << endl;
      varfile << name << myid << grid << "=" << name << myid << grid << "*2/2;" << endl;
      varfile << name << myid << grid << "=" << name << myid << grid << ".*isfinite(" << name << myid << grid << ");" << endl;
      varfile << "st = fclose(fid);" << endl;
      
    }
  varfile.close();

  if (myid == 0)
    {
      sprintf(filename,"%s%s%s","./matlabData/read",name,".m");
      varfile.open(filename,ios::out);
      for (int proc=0; proc<Communication_Manager::numberOfProcessors(); proc++)
	{
	  sprintf(filename,"%s%s%d","read",name,proc);
	  varfile << filename << endl;
	}
    }
}

void gridFunction::readFromFile(const char *name)
{
  ifstream varfile;
  char filename[120];

  int myid = Communication_Manager::My_Process_Number;
  int grid;

  for (grid=0; grid<myGridM->nrGrids(); grid++)
    {
      sprintf(filename,"%s%s%d%d","./matlabData/",name,myid,grid);
      varfile.open(filename,ios::in);
      doubleSerialArray *local = fieldValuesM[grid].getSerialArrayPointer();
      intSerialArray *localFlags = flagValuesM[grid].getSerialArrayPointer();

      for (int j=(*local).getBase(1); j<=(*local).getBound(1); j++)
	{
	  for (int i=(*local).getBase(0); i<=(*local).getBound(0); i++)
	    {
	      if ((*localFlags)(i,j) == 0)
		{
		  varfile.read((char*) &(*local)(i,j),sizeof(double));
		  (*local)(i,j) = 0.0001;
		}
	      else
		varfile.read((char*) &(*local)(i,j),sizeof(double));
	    }
	}

      varfile.close();
    }
}

void gridFunction::createHDF5File(const char *name)
{
  hid_t       file_id, file_props, group_id, group_id2, dataset_id, dataspace;
  hid_t lcpl = H5P_DEFAULT;
  herr_t      status;
  size_t s_hint = 100;
  hsize_t h_size;

  file_props = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(file_props, MPI_COMM_WORLD, MPI_INFO_NULL);
  assert(status != FAIL);
  
  //----------------------------------------
  // Create a new hdf5-file.
  //----------------------------------------
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, file_props);
  assert(file_id != FAIL);

  status = H5Pclose(file_props);
  assert(status != FAIL);
  //----------------------------------------
  // Create a root group
  //----------------------------------------
  group_id = H5Gcreate(file_id, "/root", lcpl, H5P_DEFAULT, H5P_DEFAULT);
  assert(group_id != FAIL);
  //----------------------------------------
  // ...and a overlapping grid group under 
  // that
  //----------------------------------------
  group_id2 = H5Gcreate(group_id, "overlapping grid", lcpl, H5P_DEFAULT, H5P_DEFAULT);
  assert(group_id2 != FAIL);
  //----------------------------------------
  // Write the number of component grids
  //----------------------------------------
  h_size = 1;
  dataspace = H5Screate_simple(1, &h_size, NULL);
  dataset_id = H5Dcreate(group_id2, "n_components", H5T_NATIVE_INT, dataspace, lcpl, H5P_DEFAULT, H5P_DEFAULT);
  int nOfGrids = myGridM->nrGrids();
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nOfGrids);
  assert(status != FAIL);

  status = H5Sclose(dataspace);
  assert(status != FAIL);
  status = H5Dclose(dataset_id);
  assert(status != FAIL);
  //----------------------------------------
  // Write the number of fields/grid,
  // 0 when we start
  //----------------------------------------
  h_size = 1;
  dataspace = H5Screate_simple(1, &h_size, NULL);
  dataset_id = H5Dcreate(group_id2, "n_fields", H5T_NATIVE_INT, dataspace, lcpl, H5P_DEFAULT, H5P_DEFAULT);
  int nrFields = 0;
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nrFields);
  assert(status != FAIL);

  status = H5Sclose(dataspace);
  assert(status != FAIL);
  status = H5Dclose(dataset_id);
  assert(status != FAIL);

  //----------------------------------------
  // Create a group for every component grid
  // and write flag array to file
  //----------------------------------------
  hid_t cgrid_group;
  char cgridName[200];
  for (int grid=0; grid<nOfGrids; grid++)
    {
      sprintf(cgridName, "component grid %d",grid);
      cgrid_group = H5Gcreate(group_id2, cgridName, lcpl, H5P_DEFAULT, H5P_DEFAULT);

      //----------------------------------------
      // Prepare for collective data transfer
      //----------------------------------------
      hid_t xfer_plist;         
      xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      assert(xfer_plist != FAIL);
      status = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
      assert(status != FAIL);

      hid_t file_dataspace, mem_dataspace, dataset1;
      hsize_t count[2], in_count[2], tot_count[2];
      hsize_t start[2];
      
      start[1] = flagValuesM[grid].getLocalBase(0);
      start[0] = flagValuesM[grid].getLocalBase(1);
      count[1] = flagValuesM[grid].getLocalBound(0) - start[1] + 1;
      count[0] = flagValuesM[grid].getLocalBound(1) - start[0] + 1;
      
      int bottom0, bottom1, top0, top1;

      bottom0 = ( flagValuesM[grid].getLocalBase(0) == flagValuesM[grid].getBase(0) ) ? 0 : 1;
      bottom1 = ( flagValuesM[grid].getLocalBase(1) == flagValuesM[grid].getBase(1) ) ? 0 : 1;
      top0    = ( flagValuesM[grid].getLocalBound(0) == flagValuesM[grid].getBound(0) ) ? 0 : 1;
      top1    = ( flagValuesM[grid].getLocalBound(1) == flagValuesM[grid].getBound(1) ) ? 0 : 1;
      
      in_count[0] = flagValuesM[grid].getLocalLength(1) - bottom1 - top1;
      in_count[1] = flagValuesM[grid].getLocalLength(0) - bottom0 - top0;

      tot_count[0] = flagValuesM[grid].getLength(1);
      tot_count[1] = flagValuesM[grid].getLength(0);
 
      //----------------------------------------
      // Create data [set,space] for variable
      //----------------------------------------
      file_dataspace = H5Screate_simple(2, tot_count, NULL);
      assert(file_dataspace != FAIL);

      start[0] += bottom1; start[1] += bottom0;
      status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, NULL,
				   in_count, NULL);

      assert(status != FAIL);
      
      if (bottom0 == 0) count[1]++; if (bottom1 == 0) count[0]++;
      if (top0 == 0)    count[1]++; if (top1 == 0)    count[0]++;

      mem_dataspace = H5Screate_simple (2, count, NULL);
      assert(mem_dataspace != FAIL);

      start[0] = 1; start[1] = 1;

      status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, start, NULL, 
      				   in_count, NULL);
      assert(status != FAIL);
      //----------------------------------------
      // Write variable to file
      //----------------------------------------
      dataset1 = H5Dcreate(cgrid_group, "mask", H5T_NATIVE_INT, file_dataspace, lcpl, H5P_DEFAULT, H5P_DEFAULT);
      assert(dataset1 != FAIL);

      status = H5Dwrite(dataset1, H5T_NATIVE_INT, mem_dataspace, file_dataspace,
			xfer_plist, myGridM->flagValuesM[grid].getDataPointer());
      assert(status != FAIL);

      //----------------------------------------
      // Release resources
      //----------------------------------------
      H5Dclose(dataset1);
      H5Sclose(file_dataspace);
      H5Sclose(mem_dataspace);
      H5Pclose(xfer_plist);

      status = H5Gclose(cgrid_group);
      assert(status != FAIL);
    }
  //----------------------------------------
  // Close groups
  //----------------------------------------

  status = H5Gclose(group_id2);
  assert(status != FAIL);
  status = H5Gclose(group_id);
  assert(status != FAIL);

  //----------------------------------------
  // Close hdf5-file
  //----------------------------------------
  status = H5Fclose(file_id);
  assert(status != FAIL);
}

void gridFunction::saveToHDF5File(const char *fileName, const char *varName)
{
  hid_t       file_id, file_props, group_id;
  herr_t      status;
  
  //----------------------------------------
  // Setup for collective file access
  //----------------------------------------
  file_props = H5Pcreate (H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(file_props, MPI_COMM_WORLD, MPI_INFO_NULL);
  assert(status != FAIL);
  
  //----------------------------------------
  // Open an existing hdf5-file 
  //----------------------------------------
  file_id = H5Fopen(fileName, H5F_ACC_RDWR, file_props);
  assert(file_id != FAIL);
  
  status = H5Pclose(file_props);
  assert(status != FAIL);  

  //----------------------------------------
  // Write the variable for every component
  // grid
  //----------------------------------------
  char word[200];
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
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

      hid_t file_dataspace, mem_dataspace, sid1 = FAIL, dataset1;
      hsize_t count[2], in_count[2], tot_count[2];
      hsize_t start[2];
      
      int bottom0, bottom1, top0, top1;

      bottom0 = ( flagValuesM[k].getLocalBase(0) == flagValuesM[k].getBase(0) )?0:1;
      bottom1 = ( flagValuesM[k].getLocalBase(1) == flagValuesM[k].getBase(1) )?0:1;
      top0    = ( flagValuesM[k].getLocalBound(0) == flagValuesM[k].getBound(0) )?0:1;
      top1    = ( flagValuesM[k].getLocalBound(1) == flagValuesM[k].getBound(1) )?0:1;

      start[1] = fieldValuesM[k].getLocalBase(0);
      start[0] = fieldValuesM[k].getLocalBase(1);
      count[1] = fieldValuesM[k].getLocalBound(0) - start[1] + 1;
      count[0] = fieldValuesM[k].getLocalBound(1) - start[0] + 1;
      
      in_count[0] = fieldValuesM[k].getLocalLength(1) - bottom1 - top1;
      in_count[1] = fieldValuesM[k].getLocalLength(0) - bottom0 - top0;

      tot_count[0] = fieldValuesM[k].getLength(1);
      tot_count[1] = fieldValuesM[k].getLength(0);

      //----------------------------------------
      // Check if dataset already exists
      //----------------------------------------

      //----------------------------------------
      // Turn off error printing
      //----------------------------------------
      herr_t (*old_func)(void*);
      void *old_client_data;
      H5Eget_auto1(&old_func, &old_client_data);
      hid_t error_stack;
      hid_t lcpl = H5P_DEFAULT;
      H5Eset_auto(error_stack, NULL, NULL);

      dataset1 = H5Dopen(group_id, varName, H5P_DEFAULT);      

      //----------------------------------------
      // Restore previous error handler 
      //----------------------------------------
      H5Eset_auto1(old_func, old_client_data);
      
      bool closesid1 = false;
      if (dataset1 < 0)
	{
	  //----------------------------------------
	  // Create data [set,space] for variable
	  //----------------------------------------
	  sid1 = H5Screate_simple(2, tot_count, NULL);
	  assert(sid1 != FAIL);
	  
	  dataset1 = H5Dcreate(group_id, varName, H5T_NATIVE_DOUBLE, sid1, lcpl, H5P_DEFAULT, H5P_DEFAULT);
	  assert(dataset1 != FAIL);

	  closesid1 = true;
	}
      
      file_dataspace = H5Dget_space(dataset1);   
      assert(file_dataspace != FAIL);

      start[0] += bottom1; start[1] += bottom0;
      status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, NULL,
				   in_count, NULL);
      assert(status != FAIL);
      
      if (bottom0 == 0) count[1]++; if (bottom1 == 0) count[0]++;
      if (top0 == 0)    count[1]++; if (top1 == 0)    count[0]++;

      mem_dataspace = H5Screate_simple (2, count, NULL);
      assert(mem_dataspace != FAIL);

      start[0] = 1; start[1] = 1;
      status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, start, NULL, 
				   in_count, NULL);
      assert(status != FAIL);
      //----------------------------------------
      // Write variable to file
      //----------------------------------------
      status = H5Dwrite(dataset1, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace,
      			xfer_plist, fieldValuesM[k].getDataPointer());
      assert(status != FAIL);

      //----------------------------------------
      // Release resources
      //----------------------------------------
      if (closesid1)
	{
	  H5Sclose(sid1);
	}
      H5Dclose(dataset1);
      H5Sclose(file_dataspace);
      H5Sclose(mem_dataspace);
      H5Pclose(xfer_plist);

      //----------------------------------------
      // Close component grid group
      //----------------------------------------
      status = H5Gclose(group_id);
      assert(status != FAIL);
    }
  status = H5Fclose(file_id);
  assert(status != FAIL);
}

void gridFunction::readFromHDF5File(const char *fileName, const char *varName)
{
  hid_t       file_id, file_props, group_id;
  herr_t      status;
  
  //----------------------------------------
  // Setup for collective file access
  //----------------------------------------
  file_props = H5Pcreate (H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(file_props, MPI_COMM_WORLD, MPI_INFO_NULL);
  assert(status != FAIL);
  
  //----------------------------------------
  // Open an existing hdf5-file 
  //----------------------------------------
  file_id = H5Fopen(fileName, H5F_ACC_RDONLY, file_props);
  assert(file_id != FAIL);
  
  status = H5Pclose(file_props);
  assert(status != FAIL);  

  //----------------------------------------
  // Read the variable for every component
  // grid
  //----------------------------------------
  char word[200];
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
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

      hid_t file_dataspace, mem_dataspace, dataset1;
      hsize_t count[2], in_count[2], tot_count[2];
      hsize_t start[2];

      int bottom0, bottom1, top0, top1;

      bottom0 = ( flagValuesM[k].getLocalBase(0) == flagValuesM[k].getBase(0) )?0:1;
      bottom1 = ( flagValuesM[k].getLocalBase(1) == flagValuesM[k].getBase(1) )?0:1;
      top0    = ( flagValuesM[k].getLocalBound(0) == flagValuesM[k].getBound(0) )?0:1;
      top1    = ( flagValuesM[k].getLocalBound(1) == flagValuesM[k].getBound(1) )?0:1;

      start[1] = fieldValuesM[k].getLocalBase(0);
      start[0] = fieldValuesM[k].getLocalBase(1);
      count[1] = fieldValuesM[k].getLocalBound(0) - start[1] + 1;
      count[0] = fieldValuesM[k].getLocalBound(1) - start[0] + 1;
      
      in_count[0] = fieldValuesM[k].getLocalLength(1) - bottom1 - top1;
      in_count[1] = fieldValuesM[k].getLocalLength(0) - bottom0 - top0;

      tot_count[0] = fieldValuesM[k].getLength(1);
      tot_count[1] = fieldValuesM[k].getLength(0);
      //----------------------------------------
      // Open data set for variable
      //----------------------------------------
      dataset1 = H5Dopen(group_id, varName, H5P_DEFAULT);
      assert(dataset1 != FAIL);

      file_dataspace = H5Dget_space (dataset1);   
      assert(file_dataspace != FAIL);

      start[0] += bottom1; start[1] += bottom0;
      status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, NULL,
				   in_count, NULL);
      assert(status != FAIL);

      if (bottom0 == 0) count[1]++; if (bottom1 == 0) count[0]++;
      if (top0 == 0)    count[1]++; if (top1 == 0)    count[0]++;
      
      mem_dataspace = H5Screate_simple (2, count, NULL);
      assert(mem_dataspace != FAIL);

      start[0] = 1; start[1] = 1;
      status = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, start, NULL, 
				   in_count, NULL);
      assert(status != FAIL);
      //----------------------------------------
      // Write variable to file
      //----------------------------------------
      status = H5Dread(dataset1, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace,
		       xfer_plist, fieldValuesM[k].getDataPointer());
      assert(status != FAIL);

      //----------------------------------------
      // Release resources
      //----------------------------------------
      H5Dclose(dataset1);
      H5Sclose(file_dataspace);
      H5Sclose(mem_dataspace);
      H5Pclose(xfer_plist);

      //----------------------------------------
      // Close component grid group
      //----------------------------------------
      status = H5Gclose(group_id);
      assert(status != FAIL);
    }
  status = H5Fclose(file_id);
  assert(status != FAIL);
}

void gridFunction::initializeGridData()
{
  if (boundaryValue != 0)
    for (int k=0; k<(myGridM->nrGrids()); k++)
      fieldValuesM[k] = boundaryValue(myGridM->xM[k], myGridM->yM[k], timeM, -1, -1);
  else
    {
      for (int k=0; k<(myGridM->nrGrids()); k++)
        fieldValuesM[k] = sin(0*myGridM->xM[k]);
    std::cerr << "Cannot initializeGridData, no boundaryValue function defined!!" << std::endl;
    }
}

doubleArray gridFunction::x(int k, Index ix, Index iy) const
{
  if (ix.getBase() == -1000) 
    ix = getDiscrBounds(0,k);
  
  if (iy.getBase() == -1000)
    iy = getDiscrBounds(1,k); // Indices for discretization points in the grid

  if (myGridM->gType(k) == 1)
    {
      doubleArray resArray1 = ( ( fieldValuesM[k](ix+1,iy) - fieldValuesM[k](ix-1,iy) ) / ( myGridM->xM[k](ix+1,iy) - myGridM->xM[k](ix-1,iy) ) );
      return resArray1;
    }
  else
    {
      if (conservativeDifferenceOperatorsM)
	{
	  doubleArray resArray2 = ( ( ( myGridM->ysM[k]*fieldValuesM[k] )(ix+1,iy)  - ( myGridM->ysM[k]*fieldValuesM[k] )(ix-1,iy) ) / (2*myGridM->rStep(k)) - ( ( myGridM->yrM[k]*fieldValuesM[k] )(ix,iy+1) - ( myGridM->yrM[k]*fieldValuesM[k] )(ix,iy-1) ) / (2*myGridM->sStep(k)) ) / myGridM->sqrtOfGM[k](ix,iy);
	  return resArray2;
	}
      else
	{
	  doubleArray resArray3 = ( ( myGridM->ysM[k](ix,iy)*(fieldValuesM[k](ix+1,iy) - fieldValuesM[k](ix-1,iy)) / (2*myGridM->rStep(k)) - myGridM->yrM[k](ix,iy)*(fieldValuesM[k](ix,iy+1) - fieldValuesM[k](ix,iy-1)) / (2*myGridM->sStep(k))) / myGridM->sqrtOfGM[k](ix,iy) );
	  return resArray3;
	}
    }

}

gridFunction gridFunction::x() const
{
  gridFunction result = *this;
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ix = getDiscrBounds(0,k);
      Index iy = getDiscrBounds(1,k);
      result.fieldValuesM[k].redim(ix,iy);
      result.fieldValuesM[k] = x(k);
    }
  return result;
}

doubleArray gridFunction::y(int k, Index ix, Index iy) const
{
  if (ix.getBase() == -1000) 
    ix = getDiscrBounds(0,k);
  
  if (iy.getBase() == -1000)
    iy = getDiscrBounds(1,k); // Indices for discretization points in the grid

  if (myGridM->gType(k) == 1)
    {
      doubleArray resArray1 = ( fieldValuesM[k](ix,iy+1) - fieldValuesM[k](ix,iy-1) ) / ( myGridM->yM[k](ix,iy+1) - myGridM->yM[k](ix,iy-1) );
      return resArray1;
    }
  else
    {
      if (conservativeDifferenceOperatorsM)
	{
	  doubleArray resArray2 = ( ( ( -myGridM->xsM[k]*fieldValuesM[k] )(ix+1,iy) - ( -myGridM->xsM[k]*fieldValuesM[k] )(ix-1,iy) ) / (2*myGridM->rStep(k)) + ( ( myGridM->xrM[k]*fieldValuesM[k] )(ix,iy+1) - ( myGridM->xrM[k]*fieldValuesM[k] )(ix,iy-1) ) / (2*myGridM->sStep(k)) ) / myGridM->sqrtOfGM[k](ix,iy);
	  return resArray2;
	}
      else
	{
	  doubleArray resArray3 = ( -myGridM->xsM[k](ix,iy)*(fieldValuesM[k](ix+1,iy) - fieldValuesM[k](ix-1,iy)) / (2*myGridM->rStep(k)) + myGridM->xrM[k](ix,iy)*(fieldValuesM[k](ix,iy+1) - fieldValuesM[k](ix,iy-1)) / (2*myGridM->sStep(k))) / myGridM->sqrtOfGM[k](ix,iy);
	  return resArray3;
	}
    }
  
}  

gridFunction gridFunction::y() const
{
  gridFunction result = *this;
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ix = getDiscrBounds(0,k);
      Index iy = getDiscrBounds(1,k);
      result.fieldValuesM[k].redim(ix,iy);
      result.fieldValuesM[k] = y(k);
    }
  return result;
}

doubleArray gridFunction::xx(int k, Index ix, Index iy) const
{
#define ys(x,y) (myGridM->ysM[k](x,y))
#define yr(x,y) (myGridM->yrM[k](x,y))
#define xs(x,y) (myGridM->xsM[k](x,y))
#define xr(x,y) (myGridM->xrM[k](x,y))
#define xrr(x,y) (myGridM->xrrM[k](x,y))
#define xss(x,y) (myGridM->xssM[k](x,y))
#define yrr(x,y) (myGridM->yrrM[k](x,y))
#define yss(x,y) (myGridM->yssM[k](x,y))
#define r_step (myGridM->rStep(k))
#define s_step (myGridM->sStep(k))
#define sqrtOfG(x,y) (myGridM->sqrtOfGM[k](x,y))

#define xr_A(x,y) ((myGridM->xrM[k]*fieldValuesM[k])(x,y))
#define xs_A(x,y) ((myGridM->xsM[k]*fieldValuesM[k])(x,y))
#define yr_A(x,y) ((myGridM->yrM[k]*fieldValuesM[k])(x,y))
#define ys_A(x,y) ((myGridM->ysM[k]*fieldValuesM[k])(x,y))

  if (ix.getBase() == -1000) 
    ix = getDiscrBounds(0,k);
  
  if (iy.getBase() == -1000)
    iy = getDiscrBounds(1,k); // Indices for discretization points in the grid

  if (myGridM->gType(k) == 1)
    {
      doubleArray resArray1 = ( fieldValuesM[k](ix+1,iy) - 2*fieldValuesM[k](ix,iy) + fieldValuesM[k](ix-1,iy) ) / (0.25*pow( myGridM->xM[k](ix+1,iy) - myGridM->xM[k](ix-1,iy),2. ));
      return resArray1;
    }
  else
    {
      /*    if (conservativeDifferenceOperatorsM)
	    {
	    doubleArray resArray2 = 
	    + ( + (
	    ( (ys_A(ix+1,iy) - ys_A(ix-1,iy))/(2.*r_step) - (yr_A(ix,iy+1) - yr_A(ix,iy-1))/(2.*s_step) )*ys(ix,iy) / sqrtOfG(ix,iy)
	    )(ix+1,iy)
	    - ( 
	    ( (ys_A(ix+1,iy) - ys_A(ix-1,iy))/(2.*r_step) - (yr_A(ix,iy+1) - yr_A(ix,iy-1))/(2.*s_step) )*ys(ix,iy) / sqrtOfG(ix,iy)
	    )(ix-1,iy)
	    ) / (2.*r_step)
	    + ( - (
	    ( (ys_A(ix+1,iy) - ys_A(ix-1,iy))/(2.*r_step) - (yr_A(ix,iy+1) - yr_A(ix,iy-1))/(2.*s_step) )*yr(ix,iy) / sqrtOfG(ix,iy)
	    )(ix,iy+1)
	      + ( 
	      ( (ys_A(ix+1,iy) - ys_A(ix-1,iy))/(2.*r_step) - (yr_A(ix,iy+1) - yr_A(ix,iy-1))/(2.*s_step) )*yr(ix,iy) / sqrtOfG(ix,iy)
	      )(ix,iy-1)
	      ) / (2.*s_step)
	      ;
	      return resArray2;
	      }
	      else */
      {
	doubleArray resArray3 = ( pow(ys(ix,iy),2.)*(fieldValuesM[k](ix+1,iy) - 2*fieldValuesM[k](ix,iy) + fieldValuesM[k](ix-1,iy))/pow(r_step,2.) - 2*yr(ix,iy)*ys(ix,iy)*((fieldValuesM[k](ix+1,iy+1) - fieldValuesM[k](ix+1,iy-1)) - (fieldValuesM[k](ix-1,iy+1) - fieldValuesM[k](ix-1,iy-1)))/(4.*r_step*s_step) + pow(yr(ix,iy),2.)*(fieldValuesM[k](ix,iy+1) - 2*fieldValuesM[k](ix,iy) + fieldValuesM[k](ix,iy-1))/pow(s_step,2.)) / pow(sqrtOfG(ix,iy),2.)
	  + ( (pow(ys(ix,iy),2.)*yrr(ix,iy) - 2.*yr(ix,iy)*ys(ix,iy)*(yr(ix,iy+1)-yr(ix,iy-1))/(2.*s_step) + pow(yr(ix,iy),2.)*yss(ix,iy))*(xs(ix,iy)*(fieldValuesM[k](ix+1,iy) - fieldValuesM[k](ix-1,iy))/(2.*r_step) - xr(ix,iy)*(fieldValuesM[k](ix,iy+1) - fieldValuesM[k](ix,iy-1))/(2.*s_step)) + (pow(ys(ix,iy),2.)*xrr(ix,iy) - 2.*yr(ix,iy)*ys(ix,iy)*(xr(ix,iy+1)-xr(ix,iy-1))/(2.*s_step) + pow(yr(ix,iy),2.)*xss(ix,iy))*(yr(ix,iy)*(fieldValuesM[k](ix,iy+1) - fieldValuesM[k](ix,iy-1))/(2.*s_step) - ys(ix,iy)*(fieldValuesM[k](ix+1,iy) - fieldValuesM[k](ix-1,iy))/(2.*r_step)) ) / pow(sqrtOfG(ix,iy),3.);
	return resArray3;
      }
    }
  
}
#undef ys
#undef yr
#undef xs
#undef xr
#undef xrr
#undef xss
#undef yrr
#undef yss
#undef r_step
#undef s_step
#undef sqrtOfG

#undef xr_A
#undef xs_A
#undef yr_A
#undef ys_A

gridFunction gridFunction::xx() const
{
  gridFunction result = *this;
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ix = getDiscrBounds(0,k);
      Index iy = getDiscrBounds(1,k);
      result.fieldValuesM[k].redim(ix,iy);
      result.fieldValuesM[k] = xx(k);
    }
  return result;
}

doubleArray gridFunction::yy(int k, Index ix, Index iy) const
{
#define ys(x,y) (myGridM->ysM[k](x,y))
#define yr(x,y) (myGridM->yrM[k](x,y))
#define xs(x,y) (myGridM->xsM[k](x,y))
#define xr(x,y) (myGridM->xrM[k](x,y))
#define xrr(x,y) (myGridM->xrrM[k](x,y))
#define xss(x,y) (myGridM->xssM[k](x,y))
#define yrr(x,y) (myGridM->yrrM[k](x,y))
#define yss(x,y) (myGridM->yssM[k](x,y))
#define r_step (myGridM->rStep(k))
#define s_step (myGridM->sStep(k))
#define sqrtOfG(x,y) (myGridM->sqrtOfGM[k](x,y))

  if (ix.getBase() == -1000) 
    ix = getDiscrBounds(0,k);
  
  if (iy.getBase() == -1000)
    iy = getDiscrBounds(1,k); // Indices for discretization points in the grid
  
  if (myGridM->gridTypeM[k] == 1)
    {
      doubleArray resArray1 = ( fieldValuesM[k](ix,iy+1) - 2*fieldValuesM[k](ix,iy) + fieldValuesM[k](ix,iy-1) ) / (0.25*pow( myGridM->yM[k](ix,iy+1) - myGridM->yM[k](ix,iy-1),2. ));
      return resArray1;
    }
  else
    {
      doubleArray resArray2 = ( pow(xs(ix,iy),2.)*(fieldValuesM[k](ix+1,iy) - 2*fieldValuesM[k](ix,iy) + fieldValuesM[k](ix-1,iy))/pow(r_step,2.) - 2*xr(ix,iy)*xs(ix,iy)*((fieldValuesM[k](ix+1,iy+1) - fieldValuesM[k](ix+1,iy-1)) - (fieldValuesM[k](ix-1,iy+1) - fieldValuesM[k](ix-1,iy-1)))/(4.*r_step*s_step) + pow(xr(ix,iy),2.)*(fieldValuesM[k](ix,iy+1) - 2*fieldValuesM[k](ix,iy) + fieldValuesM[k](ix,iy-1))/pow(s_step,2.)) / pow(sqrtOfG(ix,iy),2.)
	+ ( (pow(xs(ix,iy),2.)*yrr(ix,iy) - 2.*xr(ix,iy)*xs(ix,iy)*(yr(ix,iy+1)-yr(ix,iy-1))/(2.*s_step) + pow(xr(ix,iy),2.)*yss(ix,iy))*(xs(ix,iy)*(fieldValuesM[k](ix+1,iy) - fieldValuesM[k](ix-1,iy))/(2.*r_step) - xr(ix,iy)*(fieldValuesM[k](ix,iy+1) - fieldValuesM[k](ix,iy-1))/(2.*s_step)) + (pow(xs(ix,iy),2.)*xrr(ix,iy) - 2.*xr(ix,iy)*xs(ix,iy)*(xr(ix,iy+1)-xr(ix,iy-1))/(2.*s_step) + pow(xr(ix,iy),2.)*xss(ix,iy))*(yr(ix,iy)*(fieldValuesM[k](ix,iy+1) - fieldValuesM[k](ix,iy-1))/(2.*s_step) - ys(ix,iy)*(fieldValuesM[k](ix+1,iy) - fieldValuesM[k](ix-1,iy))/(2.*r_step)) ) / pow(sqrtOfG(ix,iy),3.);
      return resArray2;
    }
  
}  
#undef ys
#undef yr
#undef xs
#undef xr
#undef xrr
#undef xss
#undef yrr
#undef yss
#undef r_step
#undef s_step
#undef sqrtOfG

gridFunction gridFunction::yy() const
{
  gridFunction result = *this;
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ix = getDiscrBounds(0,k);
      Index iy = getDiscrBounds(1,k);
      result.fieldValuesM[k].redim(ix,iy);
      result.fieldValuesM[k] = yy(k);
    }
  return result;
}

doubleArray gridFunction::xy(int k, Index ix, Index iy) const
{
#define ys(x,y) (myGridM->ysM[k](x,y))
#define yr(x,y) (myGridM->yrM[k](x,y))
#define xs(x,y) (myGridM->xsM[k](x,y))
#define xr(x,y) (myGridM->xrM[k](x,y))
#define xrr(x,y) (myGridM->xrrM[k](x,y))
#define xss(x,y) (myGridM->xssM[k](x,y))
#define yrr(x,y) (myGridM->yrrM[k](x,y))
#define yss(x,y) (myGridM->yssM[k](x,y))
#define r_step (myGridM->r_stepM[k])
#define s_step (myGridM->s_stepM[k])
#define sqrtOfG(x,y) (myGridM->sqrtOfGM[k](x,y))

  if (ix.getBase() == -1000) 
    ix = getDiscrBounds(0,k);
  
  if (iy.getBase() == -1000)
    iy = getDiscrBounds(1,k); // Indices for discretization points in the grid
  
  if (myGridM->gridTypeM[k] == 1)
    {
      doubleArray resArray1 = ( (fieldValuesM[k](ix+1,iy+1) - fieldValuesM[k](ix+1,iy-1)) - (fieldValuesM[k](ix-1,iy+1) - fieldValuesM[k](ix-1,iy-1)) ) / ( (myGridM->yM[k](ix,iy+1) - myGridM->yM[k](ix,iy-1))*(myGridM->xM[k](ix+1,iy) - myGridM->xM[k](ix-1,iy)) );
      return resArray1;
    }
  else
    {
      doubleArray resArray2 = ( (xr(ix,iy)*ys(ix,iy) + xs(ix,iy)*yr(ix,iy))*((fieldValuesM[k](ix+1,iy+1) - fieldValuesM[k](ix+1,iy-1)) - (fieldValuesM[k](ix-1,iy+1) - fieldValuesM[k](ix-1,iy-1)))/(4.*r_step*s_step) - xr(ix,iy)*yr(ix,iy)*(fieldValuesM[k](ix,iy+1) - 2*fieldValuesM[k](ix,iy) + fieldValuesM[k](ix,iy-1))/pow(s_step,2.) - xs(ix,iy)*ys(ix,iy)*(fieldValuesM[k](ix+1,iy) - 2*fieldValuesM[k](ix,iy) + fieldValuesM[k](ix-1,iy))/pow(r_step,2.) ) / pow(sqrtOfG(ix,iy),2.) 
	+ ( (xr(ix,iy)*yss(ix,iy) - xs(ix,iy)*(yr(ix,iy+1)-yr(ix,iy-1))/(2.*s_step)) / pow(sqrtOfG(ix,iy),2.) + ( xs(ix,iy)*ys(ix,iy)*( (xrr(ix,iy)*ys(ix,iy) + xr(ix,iy)*(yr(ix,iy+1)-yr(ix,iy-1))/(2.*s_step)) - ((xr(ix,iy+1)-xr(ix,iy-1))/(2.*s_step)*yr(ix,iy) + xs(ix,iy)*yrr(ix,iy)) ) - xr(ix,iy)*ys(ix,iy)*( ((xr(ix,iy+1)-xr(ix,iy-1))/(2.*s_step)*ys(ix,iy) + xr(ix,iy)*yss(ix,iy)) - (xss(ix,iy)*yr(ix,iy) + xs(ix,iy)*(yr(ix,iy+1)-yr(ix,iy-1))/(2.*s_step)) ) ) / pow(sqrtOfG(ix,iy),3.) )*(fieldValuesM[k](ix+1,iy) - fieldValuesM[k](ix-1,iy))/(2.*r_step) 
	+ ( (xs(ix,iy)*yrr(ix,iy) - xr(ix,iy)*(yr(ix,iy+1)-yr(ix,iy-1))/(2.*s_step)) / pow(sqrtOfG(ix,iy),2.) + ( xr(ix,iy)*yr(ix,iy)*( ((xr(ix,iy+1)-xr(ix,iy-1))/(2.*s_step)*ys(ix,iy) + xr(ix,iy)*yss(ix,iy)) - (xss(ix,iy)*yr(ix,iy) + xs(ix,iy)*(yr(ix,iy+1)-yr(ix,iy-1))/(2.*s_step)) ) - xs(ix,iy)*yr(ix,iy)*( (xrr(ix,iy)*ys(ix,iy) + xr(ix,iy)*(yr(ix,iy+1)-yr(ix,iy-1))/(2.*s_step)) - ((xr(ix,iy+1)-xr(ix,iy-1))/(2.*s_step)*yr(ix,iy) + xs(ix,iy)*yrr(ix,iy)) ) ) / pow(sqrtOfG(ix,iy),3.) )*(fieldValuesM[k](ix,iy+1) - fieldValuesM[k](ix,iy-1))/(2.*s_step); 
      return resArray2;
    }
  
}  
#undef ys
#undef yr
#undef xs
#undef xr
#undef xrr
#undef xss
#undef yrr
#undef yss
#undef r_step
#undef s_step
#undef sqrtOfG

gridFunction gridFunction::xy() const
{
  gridFunction result = *this;
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ix = getDiscrBounds(0,k);
      Index iy = getDiscrBounds(1,k);
      result.fieldValuesM[k].redim(ix,iy);
      result.fieldValuesM[k] = xy(k);
    }
  return result;
}

doubleArray gridFunction::laplacian(int k, Index ix, Index iy) const
{
  if (ix.getBase() == -1000) 
    ix = getDiscrBounds(0,k);
  
  if (iy.getBase() == -1000)
    iy = getDiscrBounds(1,k); // Indices for discretization points in the grid
  
  doubleArray resArray = xx(k,ix,iy) + yy(k,ix,iy);
  return resArray;
}  

gridFunction gridFunction::laplacian() const
{
  gridFunction result = *this;
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ix = getDiscrBounds(0,k);
      Index iy = getDiscrBounds(1,k);
      result.fieldValuesM[k].redim(ix,iy);
      result.fieldValuesM[k] = laplacian(k);
    }
  return result;
}

void gridFunction::updateBCs(boundaryType rv)
{
  int k, side;
  
  for (k = 0; k < myGridM->nrGrids(); k++)
    {
      for (side = Side(lowR); side <= Side(hiS); side++)
	{
	  //cout << "IN updateBCs(); values are: " << endl;
          //std::cout << "BC = " << (int) getBoundaryType(k,(Side) side) << endl;
	  if (getBoundaryType(k,(Side) side) == INTERPOLATION) // interpolation
	    {
	    }
	  else if (getBoundaryType(k,(Side) side) == DIRICHLET && (rv == DIRICHLET || rv == ALL)) // Dirichlet
	    {
	      gridFunction::dirichlet(side,k);
	    }
	  else if (getBoundaryType(k,(Side) side) == NEUMANN) // Neumann
	    {
	    }
	  else if (getBoundaryType(k,(Side) side) == PERIODIC) // Periodic
	    {
	    }
	}
    }
  for (k = 0; k < myGridM->nrGrids(); k++)
    {
      for (side = Side(lowR); side <= Side(hiS); side++)
        {
          if (getBoundaryType(k,(Side) side) == INTERPOLATION) // interpolation
            {
            }
          else if (getBoundaryType(k,(Side) side) == DIRICHLET) // Dirichlet
            {
            }
          else if (getBoundaryType(k,(Side) side) == NEUMANN) // Neumann
            {
            }
          else if (getBoundaryType(k,(Side) side) == PERIODIC && (rv == PERIODIC || rv == ALL)) // Periodic
            {
              gridFunction::periodic(side,k);
            }
        }
    }
}

boundaryType gridFunction::getBoundaryType(int k, Side side)
{
  if (side == lowR)
    return bcsM[k].lowR;
  else if (side == hiR)
    return bcsM[k].hiR;
  else if (side == lowS)
    return bcsM[k].lowS;
  else if (side == hiS)
    return bcsM[k].hiS;
  else
    {
      std::cerr << "INCORRECT SIDE TYPE IN gridFunction::getBoundaryType(int k, Side side) " 
	   << endl;
      return (boundaryType) 0;
    }
}

void gridFunction::finishBCs()
{
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ixb(1,(myGridM->rDim(k))-2);
      Index low_r(0), hi_r((myGridM->rDim(k))-1);
      Index iyb(1,(myGridM->sDim(k))-2);
      Index low_s(0), hi_s((myGridM->sDim(k))-1);
      
      for (int side=Side(lowR); side<=Side(hiS); side++)
	{
	  if (getBoundaryType(k,(Side) side) == INTERPOLATION) // interpolation
	    {
	      if (side==0)
		{
		  fieldValuesM[k](low_r,iyb) = 3.*fieldValuesM[k](low_r+1,iyb) - 3.*fieldValuesM[k](low_r+2,iyb) + 1.*fieldValuesM[k](low_r+3,iyb);
		}
	      else if (side==1)
		{
		  fieldValuesM[k](hi_r,iyb) = 3.*fieldValuesM[k](hi_r-1,iyb) - 3.*fieldValuesM[k](hi_r-2,iyb) + 1.*fieldValuesM[k](hi_r-3,iyb);
		} 
	      else if (side==2)
		{
		  fieldValuesM[k](ixb,low_s) = 3.*fieldValuesM[k](ixb,low_s+1) - 3.*fieldValuesM[k](ixb,low_s+2) + 1.*fieldValuesM[k](ixb,low_s+3);
		}
	      else if (side==3)
		{
		  fieldValuesM[k](ixb,hi_s) = 3.*fieldValuesM[k](ixb,hi_s-1) - 3.*fieldValuesM[k](ixb,hi_s-2) + 1.*fieldValuesM[k](ixb,hi_s-3);
		}
	    }
	  else if (getBoundaryType(k,(Side) side) == DIRICHLET || getBoundaryType(k,(Side) side) == EXTRAPOLATION) // Dirichlet or extrapolation boundary
	    {
	      if (side==0)
		{
		  fieldValuesM[k](low_r,iyb) = 3.*fieldValuesM[k](low_r+1,iyb) - 3.*fieldValuesM[k](low_r+2,iyb) + 1.*fieldValuesM[k](low_r+3,iyb);
		}
	      else if (side==1)
		{
		  fieldValuesM[k](hi_r,iyb) = 3.*fieldValuesM[k](hi_r-1,iyb) - 3.*fieldValuesM[k](hi_r-2,iyb) + 1.*fieldValuesM[k](hi_r-3,iyb);
		} 
	      else if (side==2)
		{
		  fieldValuesM[k](ixb,low_s) = 3.*fieldValuesM[k](ixb,low_s+1) - 3.*fieldValuesM[k](ixb,low_s+2) + 1.*fieldValuesM[k](ixb,low_s+3);
		}
	      else if (side==3)
		{
		  fieldValuesM[k](ixb,hi_s) = 3.*fieldValuesM[k](ixb,hi_s-1) - 3.*fieldValuesM[k](ixb,hi_s-2) + 1.*fieldValuesM[k](ixb,hi_s-3);
		}
	      
	    }
	  else if (getBoundaryType(k,(Side) side) == NEUMANN) // Neumann
	    {
	      gridFunction::neumann(side,k);
	    }
	  else if (getBoundaryType(k,(Side) side) == PERIODIC) // Periodic
	    {
	    }
	  else
	    {
	    }
	}
      fieldValuesM[k](0,0) = 3.*fieldValuesM[k](1,1) - 3.*fieldValuesM[k](2,2) + 1.*fieldValuesM[k](3,3);
      
      fieldValuesM[k](myGridM->rdimM[k]-1,0) = 3.*fieldValuesM[k](myGridM->rdimM[k]-2,1) - 3.*fieldValuesM[k](myGridM->rdimM[k]-3,2) + 1.*fieldValuesM[k](myGridM->rdimM[k]-4,3);
      
      fieldValuesM[k](0,myGridM->sdimM[k]-1) = 3.*fieldValuesM[k](1,myGridM->sdimM[k]-2) - 3.*fieldValuesM[k](2,myGridM->sdimM[k]-3) + 1.*fieldValuesM[k](3,myGridM->sdimM[k]-4);
      
      fieldValuesM[k](myGridM->rdimM[k]-1,myGridM->sdimM[k]-1) = 3.*fieldValuesM[k](myGridM->rdimM[k]-2,myGridM->sdimM[k]-2) -3.*fieldValuesM[k](myGridM->rdimM[k]-3,myGridM->sdimM[k]-3) + 1.*fieldValuesM[k](myGridM->rdimM[k]-4,myGridM->sdimM[k]-4);
    }
}

void gridFunction::dirichlet(int i, int k)
{
  Index ixb(0, myGridM->rDim(k));
  int low_r(1), hi_r((myGridM->rDim(k)) - 2);
  Index iyb(0, myGridM->sDim(k));
  int low_s(1), hi_s((myGridM->sDim(k)) - 2);

  //cout << "In dirichlet; values are: i = " << i << ", k = " << k << endl;
  //cout << "low_r = " << low_r << endl;
  //cout << "low_s = " << low_s << endl;
  //cout << "hi_r = " << hi_r << endl;
  //cout << "hi_s = " << hi_s << endl;
  if (i == 0) // low r 
    {
      //where (myGridM->flagValuesM[k](low_r, iyb) > 0)
        fieldValuesM[k](low_r, iyb) = boundaryValue((myGridM->xM[k])(low_r, iyb), (myGridM->yM[k])(low_r, iyb), timeM, k, i);  
    }
  else if (i == 1) // high r
    {
      //where (myGridM->flagValuesM[k](hi_r, iyb) > 0)
        fieldValuesM[k](hi_r, iyb) = boundaryValue((myGridM->xM[k])(hi_r, iyb), (myGridM->yM[k])(hi_r, iyb), timeM, k, i);
    }
  else if (i == 2) // low s
    {
      //where (myGridM->flagValuesM[k](ixb, low_s) > 0)
        fieldValuesM[k](ixb, low_s) = boundaryValue((myGridM->xM[k])(ixb, low_s), (myGridM->yM[k])(ixb, low_s), timeM, k, i);
    }
  else if (i == 3) // high s
    {
      //where (myGridM->flagValuesM[k](ixb, hi_s) > 0)
        fieldValuesM[k](ixb, hi_s) = boundaryValue((myGridM->xM[k])(ixb, hi_s), (myGridM->yM[k])(ixb, hi_s), timeM, k, i);
    }
}

void gridFunction::dirichlet(int i, int k, gridFunction & F)
{
  Index ixb(0, myGridM->rDim(k));
  int low_r(1), hi_r((myGridM->rDim(k))-2);
  Index iyb(0, myGridM->sDim(k));
  int low_s(1), hi_s((myGridM->sDim(k))-2);

  if (i==0) // low r 
    {
      //where (myGridM->flagValuesM[k](low_r,iyb) > 0) 
        F.fieldValuesM[k](low_r,iyb) = boundaryValue((myGridM->xM[k])(low_r,iyb), (myGridM->yM[k])(low_r,iyb), timeM, k, i);  
    }
  else if (i==1) // high r
    {
      //where (myGridM->flagValuesM[k](hi_r,iyb) > 0)
        F.fieldValuesM[k](hi_r,iyb) = boundaryValue((myGridM->xM[k])(hi_r,iyb), (myGridM->yM[k])(hi_r,iyb), timeM, k, i) * myGridM->maskM[k](hi_r,iyb);
    }
  else if (i==2) // low s
    {
      //where (myGridM->flagValuesM[k](ixb,low_s) > 0)
        F.fieldValuesM[k](ixb,low_s) = boundaryValue((myGridM->xM[k])(ixb,low_s), (myGridM->yM[k])(ixb,low_s), timeM, k, i);
    }
  else if (i==3) // high s
    {
      //where (myGridM->flagValuesM[k](ixb,hi_s) > 0)
        F.fieldValuesM[k](ixb,hi_s) = boundaryValue((myGridM->xM[k])(ixb,hi_s), (myGridM->yM[k])(ixb,hi_s), timeM, k, i);
    }
}

void gridFunction::neumann(int i, int k)
{
#define ys(x,y) (myGridM->ysM[k](x,y))
#define yr(x,y) (myGridM->yrM[k](x,y))
#define xs(x,y) (myGridM->xsM[k](x,y))
#define xr(x,y) (myGridM->xrM[k](x,y))
#define r_step (myGridM->rStep(k))
#define s_step (myGridM->sStep(k))
#define XC(a,b) (myGridM->xM[k](a,b))
#define YC(a,b) (myGridM->yM[k](a,b))
#define sqrtOfG(x,y) (myGridM->sqrtOfGM[k](x,y))

  double dx = myGridM->xM[k](3,1) - myGridM->xM[k](1,1);
  double dy = myGridM->yM[k](1,3) - myGridM->yM[k](1,1);

  Index ixb(1,(myGridM->rDim(k))-2);
  int low_r(0), hi_r((myGridM->rDim(k))-1);
  Index iyb(1,(myGridM->sDim(k))-2);
  int low_s(0), hi_s((myGridM->sDim(k))-1);


  if (i==0) // low r 
    {
      if (myGridM->gType(k) == 1)
	{
          where (myGridM->flagValuesM[k](low_r,iyb) > 0) 
            fieldValuesM[k](low_r,iyb) = fieldValuesM[k](low_r+2,iyb) + 
	    dx*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(low_r+1,iyb),YC(low_r+1,iyb));
	}
      else
	{
          where (myGridM->flagValuesM[k](low_r,iyb) > 0) 
            fieldValuesM[k](low_r,iyb) = ( fieldValuesM[k](low_r+1,iyb-1)*r_step*(xr(low_r+1,iyb)*xs(low_r+1,iyb) + yr(low_r+1,iyb)*ys(low_r+1,iyb)) - fieldValuesM[k](low_r+1,iyb+1)*r_step*(xr(low_r+1,iyb)*xs(low_r+1,iyb) + yr(low_r+1,iyb)*ys(low_r+1,iyb)) - 2.*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(low_r+1,iyb),YC(low_r+1,iyb),xr(low_r+1,iyb),xs(low_r+1,iyb),yr(low_r+1,iyb),ys(low_r+1,iyb))*sqrtOfG(low_r+1,iyb)*r_step*s_step*sqrt(xs(low_r+1,iyb)*xs(low_r+1,iyb) + ys(low_r+1,iyb)*ys(low_r+1,iyb)) + fieldValuesM[k](low_r+2,iyb)*s_step*(xs(low_r+1,iyb)*xs(low_r+1,iyb) + ys(low_r+1,iyb)*ys(low_r+1,iyb)) )
	    / (s_step*(xs(low_r+1,iyb)*xs(low_r+1,iyb) + ys(low_r+1,iyb)*ys(low_r+1,iyb)));
	}
    }
  
  else if (i==1) // high r 
    {
      if (myGridM->gType(k) == 1)
	{
          where (myGridM->flagValuesM[k](hi_r,iyb) > 0)
            fieldValuesM[k](hi_r,iyb) = fieldValuesM[k](hi_r-2,iyb) +
	    dx*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(hi_r-1,iyb),YC(hi_r-1,iyb));
	}
      else
	{
          where (myGridM->flagValuesM[k](hi_r,iyb) > 0)
            fieldValuesM[k](hi_r,iyb) = ( -fieldValuesM[k](hi_r-1,iyb-1)*r_step*(xr(hi_r-1,iyb)*xs(hi_r-1,iyb) + yr(hi_r-1,iyb)*ys(hi_r-1,iyb)) + fieldValuesM[k](hi_r-1,iyb+1)*r_step*(xr(hi_r-1,iyb)*xs(hi_r-1,iyb) + yr(hi_r-1,iyb)*ys(hi_r-1,iyb)) + 2.*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(hi_r-1,iyb),YC(hi_r-1,iyb),xr(hi_r-1,iyb),xs(hi_r-1,iyb),yr(hi_r-1,iyb),ys(hi_r-1,iyb))*sqrtOfG(hi_r-1,iyb)*r_step*s_step*sqrt(xs(hi_r-1,iyb)*xs(hi_r-1,iyb) + ys(hi_r-1,iyb)*ys(hi_r-1,iyb)) + fieldValuesM[k](hi_r-2,iyb)*s_step*(xs(hi_r-1,iyb)*xs(hi_r-1,iyb) + ys(hi_r-1,iyb)*ys(hi_r-1,iyb)) )
	    / ( s_step*(xs(hi_r-1,iyb)*xs(hi_r-1,iyb) + ys(hi_r-1,iyb)*ys(hi_r-1,iyb)) );
	}
    }
  
  else if (i==2) // low s 
    {
      if (myGridM->gType(k) == 1)
	{
          where (myGridM->flagValuesM[k](ixb,low_s) > 0)
            fieldValuesM[k](ixb,low_s) = fieldValuesM[k](ixb,low_s+2) + 
	    dy*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(ixb,low_s+1),YC(ixb,low_s+1));
	}
      else
	{      
          where (myGridM->flagValuesM[k](ixb,low_s) > 0)
            fieldValuesM[k](ixb,low_s) = ( -2.*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(ixb,low_s+1),YC(ixb,low_s+1),xr(ixb,low_s+1),xs(ixb,low_s+1),yr(ixb,low_s+1),ys(ixb,low_s+1))*sqrtOfG(ixb,low_s+1)*r_step*s_step*sqrt(xr(ixb,low_s+1)*xr(ixb,low_s+1) + yr(ixb,low_s+1)*yr(ixb,low_s+1)) + fieldValuesM[k](ixb,low_s+2)*r_step*(xr(ixb,low_s+1)*xr(ixb,low_s+1) + yr(ixb,low_s+1)*yr(ixb,low_s+1)) + (fieldValuesM[k](ixb+1,low_s+1)-fieldValuesM[k](ixb-1,low_s+1))*s_step*(xr(ixb,low_s+1)*xs(ixb,low_s+1) + yr(ixb,low_s+1)*ys(ixb,low_s+1)) )
	    / ( r_step*(xr(ixb,low_s+1)*xr(ixb,low_s+1) + yr(ixb,low_s+1)*yr(ixb,low_s+1)) );
	}
    }
  
  else if (i==3) // high s 
    {
      if (myGridM->gType(k) == 1)
	{
          where (myGridM->flagValuesM[k](ixb,hi_s) > 0)
            fieldValuesM[k](ixb,hi_s) = fieldValuesM[k](ixb,hi_s-2) +
	    dy*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(ixb,hi_s-1),YC(ixb,hi_s-1));
	}
      else
	{
          where (myGridM->flagValuesM[k](ixb,hi_s) > 0)
            fieldValuesM[k](ixb,hi_s) = ( 2.*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(ixb,hi_s-1),YC(ixb,hi_s-1),xr(ixb,hi_s-1),xs(ixb,hi_s-1),yr(ixb,hi_s-1),ys(ixb,hi_s-1))*sqrtOfG(ixb,hi_s-1)*r_step*s_step*sqrt(xr(ixb,hi_s-1)*xr(ixb,hi_s-1) + yr(ixb,hi_s-1)*yr(ixb,hi_s-1)) + fieldValuesM[k](ixb,hi_s-2)*r_step*(xr(ixb,hi_s-1)*xr(ixb,hi_s-1) + yr(ixb,hi_s-1)*yr(ixb,hi_s-1)) - (fieldValuesM[k](ixb+1,hi_s-1)-fieldValuesM[k](ixb-1,hi_s-1))*s_step*(xr(ixb,hi_s-1)*xs(ixb,hi_s-1) + yr(ixb,hi_s-1)*ys(ixb,hi_s-1)) )
	    / ( r_step*(xr(ixb,hi_s-1)*xr(ixb,hi_s-1) + yr(ixb,hi_s-1)*yr(ixb,hi_s-1)) );
	}
    }
#undef ys
#undef yr
#undef xs
#undef xr
#undef r_step
#undef s_step
#undef XC
#undef YC
#undef sqrtOfG
}

void gridFunction::robin(int i, int k)
{
#define ys(x,y) (myGridM->ysM[k](x,y))
#define yr(x,y) (myGridM->yrM[k](x,y))
#define xs(x,y) (myGridM->xsM[k](x,y))
#define xr(x,y) (myGridM->xrM[k](x,y))
#define r_step (myGridM->rStep(k))
#define s_step (myGridM->sStep(k))
#define XC(a,b) (myGridM->xM[k](a,b))
#define YC(a,b) (myGridM->yM[k](a,b))
#define sqrtOfG(x,y) (myGridM->sqrtOfGM[k](x,y))

  double dx = myGridM->xM[k](3,1) - myGridM->xM[k](1,1);
  double dy = myGridM->yM[k](1,3) - myGridM->yM[k](1,1);

  Index ixb(1,(myGridM->rDim(k))-2);
  int low_r(0), hi_r((myGridM->rDim(k))-1);
  Index iyb(1,(myGridM->sDim(k))-2);
  int low_s(0), hi_s((myGridM->sDim(k))-1);


  if (i==0) // low r 
    {
      if (myGridM->gType(k) == 1)
	{
          where (myGridM->flagValuesM[k](low_r,iyb) > 0) 
            fieldValuesM[k](low_r,iyb) = fieldValuesM[k](low_r+2,iyb) + 
	    dx*(
		+ boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(low_r+1,iyb),YC(low_r+1,iyb))
		- fieldValuesM[k](low_r+1,iyb)
		);
	}
      else
	{
          where (myGridM->flagValuesM[k](low_r,iyb) > 0) 
            fieldValuesM[k](low_r,iyb) = ( fieldValuesM[k](low_r+1,iyb-1)*r_step*(xr(low_r+1,iyb)*xs(low_r+1,iyb) + yr(low_r+1,iyb)*ys(low_r+1,iyb)) - fieldValuesM[k](low_r+1,iyb+1)*r_step*(xr(low_r+1,iyb)*xs(low_r+1,iyb) + yr(low_r+1,iyb)*ys(low_r+1,iyb)) + (- 2.*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(low_r+1,iyb),YC(low_r+1,iyb),xr(low_r+1,iyb),xs(low_r+1,iyb),yr(low_r+1,iyb),ys(low_r+1,iyb)) + fieldValuesM[k](low_r+1,iyb) )*sqrtOfG(low_r+1,iyb)*r_step*s_step*sqrt(xs(low_r+1,iyb)*xs(low_r+1,iyb) + ys(low_r+1,iyb)*ys(low_r+1,iyb)) + fieldValuesM[k](low_r+2,iyb)*s_step*(xs(low_r+1,iyb)*xs(low_r+1,iyb) + ys(low_r+1,iyb)*ys(low_r+1,iyb)) )
	    / (s_step*(xs(low_r+1,iyb)*xs(low_r+1,iyb) + ys(low_r+1,iyb)*ys(low_r+1,iyb)));
	}
    }
  
  else if (i==1) // high r 
    {
      if (myGridM->gType(k) == 1)
	{
          where (myGridM->flagValuesM[k](hi_r,iyb) > 0)
            fieldValuesM[k](hi_r,iyb) = fieldValuesM[k](hi_r-2,iyb) +
	    dx*(
		+ boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(hi_r-1,iyb),YC(hi_r-1,iyb))
		- fieldValuesM[k](hi_r-1,iyb)
		);
	}
      else
	{
          where (myGridM->flagValuesM[k](hi_r,iyb) > 0)
            fieldValuesM[k](hi_r,iyb) = ( -fieldValuesM[k](hi_r-1,iyb-1)*r_step*(xr(hi_r-1,iyb)*xs(hi_r-1,iyb) + yr(hi_r-1,iyb)*ys(hi_r-1,iyb)) + fieldValuesM[k](hi_r-1,iyb+1)*r_step*(xr(hi_r-1,iyb)*xs(hi_r-1,iyb) + yr(hi_r-1,iyb)*ys(hi_r-1,iyb)) + (+ 2.*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(hi_r-1,iyb),YC(hi_r-1,iyb),xr(hi_r-1,iyb),xs(hi_r-1,iyb),yr(hi_r-1,iyb),ys(hi_r-1,iyb)) - fieldValuesM[k](hi_r-1,iyb) )*sqrtOfG(hi_r-1,iyb)*r_step*s_step*sqrt(xs(hi_r-1,iyb)*xs(hi_r-1,iyb) + ys(hi_r-1,iyb)*ys(hi_r-1,iyb)) + fieldValuesM[k](hi_r-2,iyb)*s_step*(xs(hi_r-1,iyb)*xs(hi_r-1,iyb) + ys(hi_r-1,iyb)*ys(hi_r-1,iyb)) )
	    / ( s_step*(xs(hi_r-1,iyb)*xs(hi_r-1,iyb) + ys(hi_r-1,iyb)*ys(hi_r-1,iyb)) );
	}
    }
  
  else if (i==2) // low s 
    {
      if (myGridM->gType(k) == 1)
	{
          where (myGridM->flagValuesM[k](ixb,low_s) > 0)
            fieldValuesM[k](ixb,low_s) = fieldValuesM[k](ixb,low_s+2) + 
	    dy*(
		+ boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(ixb,low_s+1),YC(ixb,low_s+1))
		- fieldValuesM[k](ixb,low_s+1)
		);
	}
      else
	{      
          where (myGridM->flagValuesM[k](ixb,low_s) > 0)
            fieldValuesM[k](ixb,low_s) = ( ( -2.*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(ixb,low_s+1),YC(ixb,low_s+1),xr(ixb,low_s+1),xs(ixb,low_s+1),yr(ixb,low_s+1),ys(ixb,low_s+1)) + fieldValuesM[k](ixb,low_s+1) )*sqrtOfG(ixb,low_s+1)*r_step*s_step*sqrt(xr(ixb,low_s+1)*xr(ixb,low_s+1) + yr(ixb,low_s+1)*yr(ixb,low_s+1)) + fieldValuesM[k](ixb,low_s+2)*r_step*(xr(ixb,low_s+1)*xr(ixb,low_s+1) + yr(ixb,low_s+1)*yr(ixb,low_s+1)) + (fieldValuesM[k](ixb+1,low_s+1)-fieldValuesM[k](ixb-1,low_s+1))*s_step*(xr(ixb,low_s+1)*xs(ixb,low_s+1) + yr(ixb,low_s+1)*ys(ixb,low_s+1)) )
	    / ( r_step*(xr(ixb,low_s+1)*xr(ixb,low_s+1) + yr(ixb,low_s+1)*yr(ixb,low_s+1)) );
	}
    }

  else if (i==3) // high s 
    {
      if (myGridM->gType(k) == 1)
	{
          where (myGridM->flagValuesM[k](ixb,hi_s) > 0)
            fieldValuesM[k](ixb,hi_s) = fieldValuesM[k](ixb,hi_s-2) +
	    dy*(
		+ boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(ixb,hi_s-1),YC(ixb,hi_s-1))
		- fieldValuesM[k](ixb,hi_s-1)
		);
	}
      else
	{
          where (myGridM->flagValuesM[k](ixb,hi_s) > 0)
            fieldValuesM[k](ixb,hi_s) = ( ( + 2.*boundaryDiff(i,myGridM->gridTypeM[k], timeM, XC(ixb,hi_s-1),YC(ixb,hi_s-1),xr(ixb,hi_s-1),xs(ixb,hi_s-1),yr(ixb,hi_s-1),ys(ixb,hi_s-1)) - fieldValuesM[k](ixb,hi_s-1) )*sqrtOfG(ixb,hi_s-1)*r_step*s_step*sqrt(xr(ixb,hi_s-1)*xr(ixb,hi_s-1) + yr(ixb,hi_s-1)*yr(ixb,hi_s-1)) + fieldValuesM[k](ixb,hi_s-2)*r_step*(xr(ixb,hi_s-1)*xr(ixb,hi_s-1) + yr(ixb,hi_s-1)*yr(ixb,hi_s-1)) - (fieldValuesM[k](ixb+1,hi_s-1)-fieldValuesM[k](ixb-1,hi_s-1))*s_step*(xr(ixb,hi_s-1)*xs(ixb,hi_s-1) + yr(ixb,hi_s-1)*ys(ixb,hi_s-1)) )
	    / ( r_step*(xr(ixb,hi_s-1)*xr(ixb,hi_s-1) + yr(ixb,hi_s-1)*yr(ixb,hi_s-1)) );
	}
    }
#undef ys
#undef yr
#undef xs
#undef xr
#undef r_step
#undef s_step
#undef XC
#undef YC
#undef sqrtOfG
}

void gridFunction::periodic(int i, int k)
{
  Index ixb(0,(myGridM->rDim(k))-1);
  Index low_r(0), hi_r(myGridM->rDim(k)-1);
  Index iyb(0,(myGridM->sDim(k))-1);
  Index low_s(0), hi_s(myGridM->sDim(k)-1);

  if (i == 0) // low r 
    {
      fieldValuesM[k](low_r,iyb) = fieldValuesM[k](hi_r-1,iyb);
      //    fieldValuesM[k](low_r,iyb).display();
    }
  else if (i == 1) // high r
    {
      fieldValuesM[k](hi_r,iyb) = fieldValuesM[k](low_r+1,iyb);
      //    fieldValuesM[k](hi_r,iyb).display();
    }
  else if (i == 2) // low s
    {
      fieldValuesM[k](low_s,ixb) = fieldValuesM[k](hi_s-1,ixb);
      //    fieldValuesM[k](low_s,ixb).display();
    }
  else if (i == 3) // high s
    {
      fieldValuesM[k](hi_s,ixb) = fieldValuesM[k](low_s+1,ixb);
      //    fieldValuesM[k](hi_s,ixb).display();
    }
}

/*
  Compute the time-step restriction for a convective term
  on a overlapping grid, the cfl number is normalized with
  the length of the imaginary axis included in the stability
  region of the used time integration method
*/
double gridFunction::get_CFL_Dt(const double & cfl, const gridFunction & a, const gridFunction & b, CompositeGrid & Kimera)
{
  double invDt = 0.;
  double invGridDt;

  for (int k=0; k<Kimera.nrGrids(); k++) 
    {
      Index ix = getDiscrBounds(0,k), iy = getDiscrBounds(1,k);

      if (Kimera.gType(k) == 1) // cartesian (uniform) grid
	{
	  invGridDt = 
	    max(
		+ abs(a.fieldValuesM[k](ix,iy)) 
		/ (Kimera.xM[k](ix+1,iy) - Kimera.xM[k](ix-1,iy))
		+ abs(b.fieldValuesM[k](ix,iy)) 
		/ (Kimera.yM[k](ix,iy+1) - Kimera.yM[k](ix,iy-1))
		);
	}
      else // curvi-linear grid
	{
	  invGridDt = 
	    max(
		+ abs(
		      (
		       + a.fieldValuesM[k](ix,iy)*Kimera.ysM[k](ix,iy)
		       - b.fieldValuesM[k](ix,iy)*Kimera.xsM[k](ix,iy)
		       )
		      / (2.*Kimera.r_stepM[k]) / Kimera.sqrtOfGM[k](ix,iy)
		      )
		+ abs(
		      (
		       - a.fieldValuesM[k](ix,iy)*Kimera.yrM[k](ix,iy)
		       + b.fieldValuesM[k](ix,iy)*Kimera.xrM[k](ix,iy)
		       )
		      / (2.*Kimera.s_stepM[k]) / Kimera.sqrtOfGM[k](ix,iy)
		      )
		);
	}
      
      if (invGridDt > invDt)
	invDt = invGridDt;
    }
  
  return cfl/invDt;
}

double gridFunction::getDt(const double & cfl, const gridFunction & a, const gridFunction & b, const double visc, CompositeGrid & Kimera, const double alpha0, const double beta0)
{
  double invDt = 0.;
  double invGridDt;
  
  
  for (int k=0; k<Kimera.nrGrids(); k++) 
    {
      Index ix = getDiscrBounds(0,k), iy = getDiscrBounds(1,k);
      
      if (Kimera.gType(k) == 1) // cartesian (uniform) grid
	{
	  /*      invGridDt = cfl*min(
		  pow(
		  + pow(
		  + abs( a.fieldValuesM[k](ix,iy) )*(1./(beta0*(Kimera.xM[k](ix+1,iy)-Kimera.xM[k](ix-1,iy))))
		  + abs( b.fieldValuesM[k](ix,iy) )*(1./(beta0*(Kimera.yM[k](ix,iy+1)-Kimera.yM[k](ix,iy-1)))) 
		  ,2.)
		  + pow( 
		  + visc*(4./(alpha0*pow(0.5*(Kimera.xM[k](ix+1,iy)-Kimera.xM[k](ix-1,iy)),2.)))
		  + visc*(4./(alpha0*pow(0.5*(Kimera.yM[k](ix,iy+1)-Kimera.yM[k](ix,iy-1)),2.))) 
		  ,2.)
		  , -0.5)
		  ); */
	  invGridDt = max(
			  sqrt(
			       + pow(
				     + abs( a.fieldValuesM[k](ix,iy) )*(1./(beta0*(Kimera.xM[k](ix+1,iy)-Kimera.xM[k](ix-1,iy))))
				     + abs( b.fieldValuesM[k](ix,iy) )*(1./(beta0*(Kimera.yM[k](ix,iy+1)-Kimera.yM[k](ix,iy-1)))) 
				  ,2.)
			       + pow( 
				     + visc*(4./(alpha0*pow(0.5*(Kimera.xM[k](ix+1,iy)-Kimera.xM[k](ix-1,iy)),2.)))
				     + visc*(4./(alpha0*pow(0.5*(Kimera.yM[k](ix,iy+1)-Kimera.yM[k](ix,iy-1)),2.))) 
				     ,2.)
			       )
			  );
	  if (invGridDt < invDt)
	    invDt = invGridDt;
	}
      else // curvilinear grid 
	{
	  
	  doubleArray a1, b1, nu11, nu12, nu22;
	  a1.redim(ix,iy); b1.redim(ix,iy); nu11.redim(ix,iy); nu12.redim(ix,iy); nu22.redim(ix,iy);
	  
	  a1 = -a.fieldValuesM[k](ix,iy) / Kimera.sqrtOfGM[k](ix,iy) * Kimera.ysM[k](ix,iy) - b.fieldValuesM[k](ix,iy) / Kimera.sqrtOfGM[k](ix,iy) * (-Kimera.xsM[k](ix,iy)) 
	    + visc/(pow(Kimera.sqrtOfGM[k](ix,iy),3.))*(
						       (pow(Kimera.ysM[k](ix,iy),2.)*(Kimera.yrM[k](ix+1,iy) - Kimera.yrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy)*(Kimera.yrM[k](ix,iy+1) - Kimera.yrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.yrM[k](ix,iy),2.)*(Kimera.ysM[k](ix,iy+1) - Kimera.ysM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
						       * Kimera.xsM[k](ix,iy) 
						       + (pow(Kimera.ysM[k](ix,iy),2.)*(Kimera.xrM[k](ix+1,iy) - Kimera.xrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy)*(Kimera.xrM[k](ix,iy+1) - Kimera.xrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.yrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
						       * (-Kimera.ysM[k](ix,iy)) 
						       + (pow(Kimera.xsM[k](ix,iy),2.)*(Kimera.yrM[k](ix+1,iy) - Kimera.yrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy)*(Kimera.yrM[k](ix,iy+1) - Kimera.yrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.xrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
						       * Kimera.xsM[k](ix,iy)) 
	    + (pow(Kimera.xsM[k](ix,iy),2.)*(Kimera.xrM[k](ix+1,iy) - Kimera.xrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy)*(Kimera.xrM[k](ix,iy+1) - Kimera.xrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.xrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
	    *(-Kimera.ysM[k](ix,iy)
	      );
	  
	  b1 = -a.fieldValuesM[k](ix,iy) / Kimera.sqrtOfGM[k](ix,iy) * (-Kimera.yrM[k](ix,iy)) - b.fieldValuesM[k](ix,iy) / Kimera.sqrtOfGM[k](ix,iy) * Kimera.xrM[k](ix,iy) 
	    + visc/(pow(Kimera.sqrtOfGM[k](ix,iy),3.))*(
						       (pow(Kimera.ysM[k](ix,iy),2.)*(Kimera.yrM[k](ix+1,iy) - Kimera.yrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy)*(Kimera.yrM[k](ix,iy+1) - Kimera.yrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.yrM[k](ix,iy),2.)*(Kimera.ysM[k](ix,iy+1) - Kimera.ysM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
						       * (-Kimera.xrM[k](ix,iy)) 
						       + (pow(Kimera.ysM[k](ix,iy),2.)*(Kimera.xrM[k](ix+1,iy) - Kimera.xrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy)*(Kimera.xrM[k](ix,iy+1) - Kimera.xrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.yrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
						       * Kimera.yrM[k](ix,iy) 
						       + (pow(Kimera.xsM[k](ix,iy),2.)*(Kimera.yrM[k](ix+1,iy) - Kimera.yrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy)*(Kimera.yrM[k](ix,iy+1) - Kimera.yrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.xrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
						       * (-Kimera.xrM[k](ix,iy)) 
						       + (pow(Kimera.xsM[k](ix,iy),2.)*(Kimera.xrM[k](ix+1,iy) - Kimera.xrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy)*(Kimera.xrM[k](ix,iy+1) - Kimera.xrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.xrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
						       * Kimera.yrM[k](ix,iy)
						       );
	  
	  nu11 = visc/(pow(Kimera.sqrtOfGM[k](ix,iy),2.))*(pow(Kimera.ysM[k](ix,iy),2.) + pow(Kimera.xsM[k](ix,iy),2.));
	  
	  nu12 = visc/(pow(Kimera.sqrtOfGM[k](ix,iy),2.))*(-2.*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy) -2.*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy));
	  
	  nu22 = visc/(pow(Kimera.sqrtOfGM[k](ix,iy),2.))*(pow(Kimera.yrM[k](ix,iy),2.) + pow(Kimera.xrM[k](ix,iy),2.));
      
	  invGridDt = max(
			  sqrt(
			       pow( abs(a1)*(1./(beta0*2.*Kimera.r_stepM[k]))+abs(b1)*(1./(beta0*2.*Kimera.s_stepM[k])) ,2.)
			       + pow( nu11*(4./(alpha0*Kimera.r_stepM[k]*Kimera.r_stepM[k]))
				      + abs(nu12)*(1./(alpha0*Kimera.r_stepM[k]*Kimera.s_stepM[k]))
				      + nu22*(4./(alpha0*Kimera.s_stepM[k]*Kimera.s_stepM[k])) ,2.)
			       )
			  );
	}
      if (invGridDt > invDt)
	invDt = invGridDt;
    }
  
  return cfl/invDt;
}

double gridFunction::getDt(const double & cfl, const gridFunction & a, const gridFunction & b, const gridFunction & visc, CompositeGrid & Kimera, const double alpha0, const double beta0)
{
  double Dt = 10000.;
  double gridDt;
  
  
  for (int k=0; k<Kimera.nrGrids(); k++) 
    {
      Index ix = getDiscrBounds(0,k), iy = getDiscrBounds(1,k);
      
      if (Kimera.gType(k) == 1) // cartesian (uniform) grid
	{
	  gridDt = cfl*min(
			   pow(
			       pow( abs(a.fieldValuesM[k](ix,iy))*(1./(beta0*(Kimera.xM[k](ix+1,iy)-Kimera.xM[k](ix-1,iy)))) + abs(b.fieldValuesM[k](ix,iy))*(1./(beta0*(Kimera.yM[k](ix,iy+1)-Kimera.yM[k](ix,iy-1)))) ,2.)
	  + pow( visc.fieldValuesM[k](ix,iy)*(4./(alpha0*pow(0.5*(Kimera.xM[k](ix+1,iy)-Kimera.xM[k](ix-1,iy)),2.)))
		 + visc.fieldValuesM[k](ix,iy)*(4./(alpha0*pow(0.5*(Kimera.yM[k](ix,iy+1)-Kimera.yM[k](ix,iy-1)),2.))) ,2.)
			       , -0.5)
			   );
	  if (gridDt < Dt)
	    Dt = gridDt;
	}
      else // curvilinear grid 
	{
	  
	  doubleArray a1, b1, nu11, nu12, nu22;
	  a1.redim(ix,iy); b1.redim(ix,iy); nu11.redim(ix,iy); nu12.redim(ix,iy); nu22.redim(ix,iy);
	  
	  a1 = -a.fieldValuesM[k](ix,iy) / Kimera.sqrtOfGM[k](ix,iy) * Kimera.ysM[k](ix,iy) - b.fieldValuesM[k](ix,iy) / Kimera.sqrtOfGM[k](ix,iy) * (-Kimera.xsM[k](ix,iy)) 
	    + visc.fieldValuesM[k](ix,iy)/(pow(Kimera.sqrtOfGM[k](ix,iy),3.))*(
									     (pow(Kimera.ysM[k](ix,iy),2.)*(Kimera.yrM[k](ix+1,iy) - Kimera.yrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy)*(Kimera.yrM[k](ix,iy+1) - Kimera.yrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.yrM[k](ix,iy),2.)*(Kimera.ysM[k](ix,iy+1) - Kimera.ysM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
									     * Kimera.xsM[k](ix,iy) 
									     + (pow(Kimera.ysM[k](ix,iy),2.)*(Kimera.xrM[k](ix+1,iy) - Kimera.xrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy)*(Kimera.xrM[k](ix,iy+1) - Kimera.xrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.yrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
									     * (-Kimera.ysM[k](ix,iy)) 
									     + (pow(Kimera.xsM[k](ix,iy),2.)*(Kimera.yrM[k](ix+1,iy) - Kimera.yrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy)*(Kimera.yrM[k](ix,iy+1) - Kimera.yrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.xrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
									     * Kimera.xsM[k](ix,iy)) 
	    + (pow(Kimera.xsM[k](ix,iy),2.)*(Kimera.xrM[k](ix+1,iy) - Kimera.xrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy)*(Kimera.xrM[k](ix,iy+1) - Kimera.xrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.xrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
	    *(-Kimera.ysM[k](ix,iy)
	      );
	  
	  b1 = -a.fieldValuesM[k](ix,iy) / Kimera.sqrtOfGM[k](ix,iy) * (-Kimera.yrM[k](ix,iy)) - b.fieldValuesM[k](ix,iy) / Kimera.sqrtOfGM[k](ix,iy) * Kimera.xrM[k](ix,iy) 
	    + visc.fieldValuesM[k](ix,iy)/(pow(Kimera.sqrtOfGM[k](ix,iy),3.))*(
									     (pow(Kimera.ysM[k](ix,iy),2.)*(Kimera.yrM[k](ix+1,iy) - Kimera.yrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy)*(Kimera.yrM[k](ix,iy+1) - Kimera.yrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.yrM[k](ix,iy),2.)*(Kimera.ysM[k](ix,iy+1) - Kimera.ysM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
									     * (-Kimera.xrM[k](ix,iy)) 
									     + (pow(Kimera.ysM[k](ix,iy),2.)*(Kimera.xrM[k](ix+1,iy) - Kimera.xrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy)*(Kimera.xrM[k](ix,iy+1) - Kimera.xrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.yrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
									     * Kimera.yrM[k](ix,iy) 
									     + (pow(Kimera.xsM[k](ix,iy),2.)*(Kimera.yrM[k](ix+1,iy) - Kimera.yrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy)*(Kimera.yrM[k](ix,iy+1) - Kimera.yrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.xrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
									     * (-Kimera.xrM[k](ix,iy)) 
									     + (pow(Kimera.xsM[k](ix,iy),2.)*(Kimera.xrM[k](ix+1,iy) - Kimera.xrM[k](ix-1,iy))/(2.*Kimera.r_stepM[k]) - 2*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy)*(Kimera.xrM[k](ix,iy+1) - Kimera.xrM[k](ix,iy-1))/(2.*Kimera.s_stepM[k]) + pow(Kimera.xrM[k](ix,iy),2.)*(Kimera.xsM[k](ix,iy+1) - Kimera.xsM[k](ix,iy-1))/(2.*Kimera.s_stepM[k])) 
									     * Kimera.yrM[k](ix,iy)
									     );
	  
	  nu11 = visc.fieldValuesM[k](ix,iy)/(pow(Kimera.sqrtOfGM[k](ix,iy),2.))*(pow(Kimera.ysM[k](ix,iy),2.) + pow(Kimera.xsM[k](ix,iy),2.));
	  
	  nu12 = visc.fieldValuesM[k](ix,iy)/(pow(Kimera.sqrtOfGM[k](ix,iy),2.))*(-2.*Kimera.yrM[k](ix,iy)*Kimera.ysM[k](ix,iy) -2.*Kimera.xrM[k](ix,iy)*Kimera.xsM[k](ix,iy));
	  
	  nu22 = visc.fieldValuesM[k](ix,iy)/(pow(Kimera.sqrtOfGM[k](ix,iy),2.))*(pow(Kimera.yrM[k](ix,iy),2.) + pow(Kimera.xrM[k](ix,iy),2.));
      
	  gridDt = cfl*min(
			   pow(
			       pow( abs(a1)*(1./(beta0*Kimera.r_stepM[k]))+abs(b1)*(1./(beta0*Kimera.s_stepM[k])) ,2.)
			       + pow( nu11*(4./(alpha0*Kimera.r_stepM[k]*Kimera.r_stepM[k]))
				      + abs(nu12)*(1./(alpha0*Kimera.r_stepM[k]*Kimera.s_stepM[k]))
				      + nu22*(4./(alpha0*Kimera.s_stepM[k]*Kimera.s_stepM[k])) ,2.)
	  , -0.5)
			   );
	}
      if (gridDt < Dt)
	Dt = gridDt;
    }
  
  return Dt;
}


void gridFunction::updateGhostBoundaries()
{
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      fieldValuesM[k].updateGhostBoundaries();
    }
}


doubleSerialArray gridFunction::computeInterpolationPoints(int k)
#define i_L (myGridM->intInterpolationLocationM[k](all,0))
#define j_L (myGridM->intInterpolationLocationM[k](all,1))
{
  
  Index all;
  
  int sumOther = i_L.getLength(0);
  
  if (iTypeM == cubic)
    {
      if (sumOther > 0)
	{
	  doubleSerialArray U;
	  U = fieldValuesM[k].getLocalArrayWithGhostBoundaries();

	  doubleSerialArray epsilonM1 = epsilon(-1,k);
	  doubleSerialArray epsilon0 = epsilon(0,k);
	  doubleSerialArray epsilonP1 = epsilon(1,k);
	  doubleSerialArray epsilonP2 = epsilon(2,k);
	  doubleSerialArray zetaM1 = zeta(-1,k);
	  doubleSerialArray zeta0 = zeta(0,k);
	  doubleSerialArray zetaP1 = zeta(1,k);
	  doubleSerialArray zetaP2 = zeta(2,k);

	  doubleSerialArray result0 = epsilonM1*zetaM1*U(i_L,j_L) + epsilon0*zetaM1*U(i_L+1,j_L)
	    + epsilonP1*zetaM1*U(i_L+2,j_L) + epsilonP2*zetaM1*U(i_L+3,j_L)
	    + epsilonM1*zeta0*U(i_L,j_L+1) + epsilon0*zeta0*U(i_L+1,j_L+1) 
	    + epsilonP1*zeta0*U(i_L+2,j_L+1) + epsilonP2*zeta0*U(i_L+3,j_L+1) 
	    + epsilonM1*zetaP1*U(i_L,j_L+2) + epsilon0*zetaP1*U(i_L+1,j_L+2) 
	    + epsilonP1*zetaP1*U(i_L+2,j_L+2) + epsilonP2*zetaP1*U(i_L+3,j_L+2)
	    + epsilonM1*zetaP2*U(i_L,j_L+3) + epsilon0*zetaP2*U(i_L+1,j_L+3) 
	    + epsilonP1*zetaP2*U(i_L+2,j_L+3) + epsilonP2*zetaP2*U(i_L+3,j_L+3);
	  
	  return result0;
	}
    }
  else if (iTypeM == quadratic) 
    {
      if (sumOther > 0)
	{
	  doubleSerialArray U;
	  U = fieldValuesM[k].getLocalArrayWithGhostBoundaries();

	  doubleSerialArray alphaM1 = alpha(-1,k);
	  doubleSerialArray alpha0 = alpha(0,k);
	  doubleSerialArray alphaP1 = alpha(1,k);
	  doubleSerialArray betaM1 = beta(-1,k);
	  doubleSerialArray beta0 = beta(0,k);
	  doubleSerialArray betaP1 = beta(1,k);
	  

 	  doubleSerialArray result1 = alphaM1*betaM1*U(i_L,j_L) + alpha0*betaM1*U(i_L+1,j_L)
 	    + alphaP1*betaM1*U(i_L+2,j_L) + alphaM1*beta0*U(i_L,j_L+1) 
 	    + alpha0*beta0*U(i_L+1,j_L+1) + alphaP1*beta0*U(i_L+2,j_L+1) 
 	    + alphaM1*betaP1*U(i_L,j_L+2) + alpha0*betaP1*U(i_L+1,j_L+2) 
 	    + alphaP1*betaP1*U(i_L+2,j_L+2);
	
	  return result1;
	}
    }
  else if (iTypeM == linear)
    {
      doubleSerialArray U;
      U = fieldValuesM[k].getLocalArrayWithGhostBoundaries();

      doubleSerialArray result3 = gamma(-1,k)*delta(-1,k)*U(i_L,j_L) + gamma(0,k)*delta(-1,k)*U(i_L+1,j_L) + gamma(-1,k)*delta(0,k)*U(i_L,j_L+1) + gamma(0,k)*delta(0,k)*U(i_L+1,j_L+1);

      return result3;
    }
  doubleSerialArray temp;
  doubleSerialArray & result2 = temp;
  return result2;
}
#undef i_L
#undef j_L

doubleSerialArray gridFunction::alpha(int i,int k)
#define rstep (myGridM->rStep(k))
#define r_L (myGridM->interpolationCoordinatesM[k](all,0))
#define r (myGridM->doubleInterpolationLocationM[k](all,0))
{
  Index all;

  if (i == -1)
    {
      doubleSerialArray result1 = (r_L - (r+1-1)*rstep)*(r_L - (r+2-1)*rstep) / (2*rstep*rstep);
      return result1;
    }
  else if (i == 0)
    {
      doubleSerialArray result2 = - (r_L - (r-1)*rstep)*(r_L - (r+2-1)*rstep) / (rstep*rstep);
      return result2;
    }
  else 
    {
      doubleSerialArray result3 = (r_L - (r-1)*rstep)*(r_L - (r+1-1)*rstep) / (2*rstep*rstep);
      return result3;
    }
  
}
#undef rstep
#undef r_L
#undef r

doubleSerialArray gridFunction::beta(int i,int k)
#define sstep (myGridM->sStep(k))  
#define s_L (myGridM->interpolationCoordinatesM[k](all,1))
#define s (myGridM->doubleInterpolationLocationM[k](all,1))  
{
  Index all;

  if (i == -1) 
    {
      doubleSerialArray result1 = (s_L - (s+1-1)*sstep)*(s_L - (s+2-1)*sstep) / (2*sstep*sstep);
      return result1;
    }
  else if (i == 0)
    {
      doubleSerialArray result2 = - (s_L - (s-1)*sstep)*(s_L - (s+2-1)*sstep) / (sstep*sstep);
      return result2;
    }
  else 
    {
      doubleSerialArray result3 = (s_L - (s-1)*sstep)*(s_L - (s+1-1)*sstep) / (2*sstep*sstep);
      return result3;
    }

}
#undef sstep
#undef s_L
#undef s

doubleSerialArray gridFunction::gamma(int i,int k)
#define rstep (myGridM->rStep(k))
#define r_L (myGridM->interpolationCoordinatesM[k](all,0))
#define r (myGridM->doubleInterpolationLocationM[k](all,0))
{
  Index all;

  if (i == -1)
    {
      doubleSerialArray result1 = 1 + (r_L - (r-1)*rstep) / (-rstep);
      return result1;
    }
  else 
    {
      doubleSerialArray result2 = (r_L - (r-1)*rstep) / (rstep);
      return result2;
    }
}
#undef rstep
#undef r_L
#undef r

doubleSerialArray gridFunction::delta(int i,int k)
#define sstep (myGridM->sStep(k))  
#define s_L (myGridM->interpolationCoordinatesM[k](all,1))
#define s (myGridM->doubleInterpolationLocationM[k](all,1))  
{
  Index all;

  if (i == -1) 
    {
      doubleSerialArray result1 = 1 + (s_L - (s-1)*sstep) / (-sstep);
      return result1;
    }
  else 
    {
      doubleSerialArray result2 = (s_L - (s-1)*sstep) / (sstep);
      return result2;
    }
}
#undef sstep
#undef s_L
#undef s

doubleSerialArray gridFunction::epsilon(int i,int k)
#define rstep (myGridM->rStep(k))
#define r_L (myGridM->interpolationCoordinatesM[k](all,0))
#define r (myGridM->doubleInterpolationLocationM[k](all,0))
{
  Index all;

  if (i == -1)
    {
      doubleSerialArray result1 = - (r_L - (r+1-1)*rstep)*(r_L - (r+2-1)*rstep)*(r_L - (r+3-1)*rstep) / (6*rstep*rstep*rstep);
      return result1;
    }
  else if (i == 0)
    {
      doubleSerialArray result2 = (r_L - (r-1)*rstep)*(r_L - (r+2-1)*rstep)*(r_L - (r+3-1)*rstep) / (2*rstep*rstep*rstep);
      return result2;
    }
  else if (i == 1)
    {
      doubleSerialArray result3 = - (r_L - (r-1)*rstep)*(r_L - (r+1-1)*rstep)*(r_L - (r+3-1)*rstep) / (2*rstep*rstep*rstep);
      return result3;
    }
  else
    {
      doubleSerialArray result4 = (r_L - (r-1)*rstep)*(r_L - (r+1-1)*rstep)*(r_L - (r+2-1)*rstep) / (6*rstep*rstep*rstep);
      return result4;
    } 
}
#undef rstep
#undef r_L
#undef r

doubleSerialArray gridFunction::zeta(int i,int k)
#define sstep (myGridM->sStep(k))  
#define s_L (myGridM->interpolationCoordinatesM[k](all,1))
#define s (myGridM->doubleInterpolationLocationM[k](all,1))  
{
  Index all;

  if (i == -1) 
    {
      doubleSerialArray result1 = (s_L - (s+1-1)*sstep)*(s_L - (s+2-1)*sstep)*(s_L - (s+3-1)*sstep) / (6*sstep*sstep*sstep);
      return result1;
    }
  else if (i == 0)
    {
      doubleSerialArray result2 = (s_L - (s-1)*sstep)*(s_L - (s+2-1)*sstep)*(s_L - (s+3-1)*sstep) / (2*sstep*sstep*sstep);
      return result2;
    }
  else if (i == 1)
    {
      doubleSerialArray result3 = - (s_L - (s-1)*sstep)*(s_L - (s+1-1)*sstep)*(s_L - (s+3-1)*sstep) / (2*sstep*sstep*sstep);
      return result3;
    }
  else
    { 
      doubleSerialArray result4 = (s_L - (s-1)*sstep)*(s_L - (s+1-1)*sstep)*(s_L - (s+2-1)*sstep) / (6*sstep*sstep*sstep);
      return result4;
    }
}
#undef sstep
#undef s_L
#undef s

void gridFunction::interpolate(int gridToUpdate, gridFunction *target)
{
  int k, p, grid, proc, nrProcs, myRank, currGrid;
  Index all;

  /*
    First check if we really need to interpolate between component grids,
    i.e. is there more than one component grid
  */
  if (myGridM->nrGrids() == 1)
    return;

  MPI_Comm myComm;
  MPI_Group myGroup;
  
  MPI_Comm_group(MPI_COMM_WORLD, &myGroup);
  MPI_Comm_create(MPI_COMM_WORLD, myGroup, &myComm);
  MPI_Comm_size(myComm, &nrProcs);
  MPI_Comm_rank(myComm, &myRank);


  double *theSendArray, *theReceiveArray;
  theSendArray = myGridM->theSendBufferM.Array_Descriptor.Array_Data;
  theReceiveArray = myGridM->theReceiveBufferM.Array_Descriptor.Array_Data;
  
  MPI_Request *Request1 = new MPI_Request[nrProcs];
  MPI_Request *Request2 = new MPI_Request[nrProcs];

  /*
    First post all receives needed
  */

  int recHandle;
  int currIndex=0;
  int *block_lengths = new int[(myGridM->nmbrOfGridsM)*(myGridM->nmbrOfGridsM)];
  int *displacements = new int[(myGridM->nmbrOfGridsM)*(myGridM->nmbrOfGridsM)];

  for (p=0; p<nrProcs; p++)
    {
      for (k=0; k<myGridM->nrGrids(); k++)
	{
	  displacements[k] = 0;
	  for (grid=0; grid<k; grid++)
	    displacements[k] += myGridM->myIntPointsM[grid].getLength(0);
	  
	  for (proc=0; proc<p; proc++)
	    displacements[k] += myGridM->myNrOfIntPointsM(k,proc);
	  
	  block_lengths[k] = myGridM->myNrOfIntPointsM(k,p);
	}
      
      MPI_Type_indexed(myGridM->nrGrids(), block_lengths, displacements, MPI_DOUBLE, &(myGridM->receiveTypeM[p]));
      MPI_Type_commit(&(myGridM->receiveTypeM[p]));
      
      recHandle = 4;
      
      double *rpoint = 0x0;
      if (theReceiveArray == 0) // message of len=0
	MPI_Irecv(rpoint, 1, myGridM->receiveTypeM[p], p, recHandle, myComm, &(Request2[p]));
      else
	MPI_Irecv(theReceiveArray, 1, myGridM->receiveTypeM[p], p, recHandle, myComm, &(Request2[p]));
    }
  
  currIndex = 0;
  myGridM->offsetsM = 0;
  for (k=0; k<myGridM->nrGrids(); k++)
    {
      for (p=0; p<nrProcs; p++)
	{
	  for (currGrid=0; currGrid<myGridM->nrGrids(); currGrid++)
	    {
	      for (grid=0; grid<k; grid++)
		myGridM->offsetsM(p,k,currGrid) += myGridM->intInterpolationLocationM[grid].getLength(0);
	      
	      for (proc=0; proc<p; proc++)
		for (grid=0; grid<myGridM->nrGrids(); grid++)
		  myGridM->offsetsM(p,k,currGrid) += myGridM->nrOfPointsToComputeM(k,grid,proc);
	      
	      for (grid=0; grid<currGrid; grid++)
		myGridM->offsetsM(p,k,currGrid) += myGridM->nrOfPointsToComputeM(k,grid,p);
	      
	    }
      
	}
      Index sendIndex(currIndex, (myGridM->intInterpolationLocationM[k]).getLength(0));
      myGridM->theSendBufferM(sendIndex) = computeInterpolationPoints(k);
      currIndex += (myGridM->intInterpolationLocationM[k]).getLength(0);
    }
  
  /*
    Send all data to correct process
  */
  int sendHandle;
  for (p=0; p<nrProcs; p++)
    {
      int index = 0;
      for (int i=0; i<myGridM->nrGrids(); i++)
	for (int j=0; j<myGridM->nrGrids(); j++)
	  {
	    block_lengths[index] = myGridM->nrOfPointsToComputeM(j,i,p);
	    displacements[index] = myGridM->offsetsM(p,j,i);
	    index++;
	  }
      
      MPI_Type_indexed((myGridM->nrGrids())*(myGridM->nrGrids()),block_lengths,displacements,MPI_DOUBLE, &(myGridM->sendTypeM[p]));
      MPI_Type_commit(&(myGridM->sendTypeM[p]));
      
      sendHandle = 4;
      double *spoint = 0x0;
      if (theSendArray == 0) // message of len=0
	MPI_Isend(spoint, 1, myGridM->sendTypeM[p], p, sendHandle, myComm, &(Request1[p]));
      else
	MPI_Isend(theSendArray, 1, myGridM->sendTypeM[p], p, sendHandle, myComm, &(Request1[p]));
    }
  
  delete[] block_lengths;
  delete[] displacements;
  
  MPI_Status *status1 = new MPI_Status[nrProcs];
  MPI_Status *status2 = new MPI_Status[nrProcs];
  
  MPI_Waitall(nrProcs, Request2, status2);
  MPI_Waitall(nrProcs, Request1, status1);

  /*
    Finally insert interpolated elements into
    correct places in the gridFunction
  */
  if (target == 0)
    {
      target = this;
    }

  fillIntElements(myGridM->theReceiveBufferM, target, gridToUpdate);
  
  for (p=0; p<nrProcs; p++)
    {
      MPI_Type_free(&(myGridM->sendTypeM[p]));
      MPI_Type_free(&(myGridM->receiveTypeM[p]));
    }
  delete[] Request1;
  delete[] Request2;
  delete[] status1;
  delete[] status2;

  MPI_Comm_free(&myComm);
}

void gridFunction::fillIntElements(doubleSerialArray & elementArray, gridFunction *target, int gridToUpdate)
#define i_point (myGridM->myIntPointsM[k](all,0))
#define j_point (myGridM->myIntPointsM[k](all,1))
{
  Index all;
  int startElement, endElement, nrProcs;

  MPI_Comm_size(MPI_COMM_WORLD, &nrProcs);

  startElement=0;
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      if (k==gridToUpdate || gridToUpdate==-1)
	{
	  endElement = startElement;
	  doubleSerialArray *U;
	  U = target->fieldValuesM[k].getSerialArrayPointer();
	  for (int p=0; p<nrProcs; p++)
	    endElement += myGridM->myNrOfIntPointsM(k,p);
	  
	  if (endElement > startElement)
	    {
	      Range elementRange(startElement, endElement-1);
	      (*U)(i_point,j_point) = elementArray(elementRange);
	    }
	  startElement = endElement;
	}
    }
}
#undef i_point
#undef j_point

void gridFunction::applyMask(int k)
{
  fieldValuesM[k] *= myGridM->maskM[k];
}

void gridFunction::addTime(double add)
{
  timeM += add;
}

double gridFunction::getTime() const
{
  return timeM;
}

void gridFunction::checkErrors(doubleArray (*infunc) (const doubleArray & x, const doubleArray & y, const double t, int side, int grid))
{
  /*  double maxError = 0.0;

  for (int k=0; k<myGridM->nrGrids(); k++) 
  {
    Index ix = getDiscrBounds(0,k);
    Index iy = getDiscrBounds(1,k); 

    doubleArray gridError = (fabs((fieldValuesM[k](ix,iy) - infunc(myGridM->xM[k](ix,iy),myGridM->yM[k](ix,iy),time,-1,-1)) *myGridM->maskM[k](ix,iy)));

    double TempmaxError = max(gridError);
    maxError = (TempmaxError > maxError) ? TempmaxError : maxError;

  }
  if (Communication_Manager::My_Process_Number == 0) 
    std::cout << " MaxError is " << maxError << endl;

  */
  int k;

  for (k=0; k<myGridM->nrGrids(); k++) 
    {
      fieldValuesM[k] -= infunc(myGridM->xM[k], myGridM->yM[k], timeM, -1, -1);
    }
  double linf = l_inf();
  double lh = l_h();
  if (Communication_Manager::My_Process_Number == 0) 
    {
      std::cout << "  l_inf norm of Error is " << linf << endl;
      std::cout << "  l_h norm of Error is   " << lh << endl;
    }
  for (k=0; k<myGridM->nrGrids(); k++) 
    {
      fieldValuesM[k] += infunc(myGridM->xM[k], myGridM->yM[k], timeM, -1, -1);
    }
}

void gridFunction::PUMcheckErrors(doubleArray (*infunc) (const doubleArray & x, const doubleArray & y, const double t, int side, int grid),bool check)
{
  double R = 0.5;
  double overlap = 0.25;
  
  if (check)
    fieldValuesM[0] -= infunc(myGridM->xM[0], myGridM->yM[0], timeM, -1, -1);
  doubleArray radius0 = (sqrt(pow(myGridM->xM[0],2.)+pow(myGridM->yM[0],2.))-R)/overlap;
  where (radius0 <= 0)
    fieldValuesM[0] = 0;
  if (check)
    fieldValuesM[1] -= infunc(myGridM->xM[1],myGridM->yM[1], timeM, -1, -1);
  doubleArray radius1 = (sqrt(pow(myGridM->xM[1],2.)+pow(myGridM->yM[1],2.))-R)/overlap;
  where (radius1 >= 1)
    fieldValuesM[1] = 0;
  double linf = l_inf();
  double lh = l_h();
  if (Communication_Manager::My_Process_Number == 0 && check) 
    {
      std::cout << "  l_inf norm of Error is " << linf << endl;
      std::cout << "  l_h norm of Error is   " << lh << endl;
    }
  radius0 = (sqrt(pow(myGridM->xM[0],2.)+pow(myGridM->yM[0],2.))-R)/overlap;
  where (radius0 <= 0)
    fieldValuesM[0] = 10000000;

  radius1 = (sqrt(pow(myGridM->xM[1],2.)+pow(myGridM->yM[1],2.))-R)/overlap;
  where (radius1 >= 1)
    fieldValuesM[1] = 10000000;

}

doubleArray gridFunction::forcing(doubleArray (*infunc) (const doubleArray & x, const doubleArray & y, const double t, int grid, int side), int grid, double tadd, Index ix, Index iy)
{
  if (ix == -1000)
    ix = getDiscrBounds(0,grid);
  if (iy == -1000) 
    iy = getDiscrBounds(1,grid);  

  doubleArray force = infunc(myGridM->xM[grid](ix,iy), myGridM->yM[grid](ix,iy), timeM + tadd, -1, -1)*myGridM->maskM[grid](ix,iy);
  
  return force;
  
}

gridFunction gridFunction::forcing(doubleArray (*infunc) (const doubleArray & x, const doubleArray & y, const double t, int grid, int side), double tadd)
{
  gridFunction result = *this;
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ix = getDiscrBounds(0,k);
      Index iy = getDiscrBounds(1,k);
      result.fieldValuesM[k].redim(ix,iy);
      result.fieldValuesM[k] = infunc(myGridM->xM[k](ix,iy), myGridM->yM[k](ix,iy), timeM + tadd, -1, -1)*myGridM->maskM[k](ix,iy);
    }
  return result;
}

void gridFunction::setBoundaryConditionType(int k, Side side, boundaryType condition)
{
  if (side == lowR)
    bcsM[k].lowR = condition;
  else if (side == hiR)
    bcsM[k].hiR = condition;
  else if (side == lowS)
    bcsM[k].lowS = condition;
  else if (side == hiS)
    bcsM[k].hiS = condition;
}

void gridFunction::setBC(functionType fType)
{
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      for (int side=Side(lowR); side<=Side(hiS); side++)
	{
	  if (fType == DirichletAllOver && myGridM->getBoundaryType(k,(Side) side) != PERIODIC)
	    {
	      setBoundaryConditionType(k,(Side) side, DIRICHLET);
	    }
	  else if (fType == ExtrapolationAllOver && myGridM->getBoundaryType(k,(Side) side) != PERIODIC)
	    {
	      setBoundaryConditionType(k,(Side) side, EXTRAPOLATION);
	    }
	  else if (myGridM->getBoundaryType(k,(Side) side) == INTERPOLATION) // interpolation
	    {
	      setBoundaryConditionType(k,(Side) side, INTERPOLATION);
	    }
	  else if (myGridM->getBoundaryType(k,(Side) side) == NOSLIP) // no-slip or inflow
	    {
	      if (fType == scalar)
		setBoundaryConditionType(k,(Side) side, NEUMANN);
	      else if (fType == xComponentOfVector || fType == yComponentOfVector)
		setBoundaryConditionType(k,(Side) side, DIRICHLET);
	    }
	  else if (myGridM->getBoundaryType(k,(Side) side) == PERIODIC) // periodic
	    {
	      setBoundaryConditionType(k,(Side) side, PERIODIC);
	    }
	  else if (myGridM->getBoundaryType(k,(Side) side) == SLIP) // slip wall
	    {
	      if (fType == scalar)
		setBoundaryConditionType(k,(Side) side, NEUMANN);
	      else if (fType == xComponentOfVector && ( (Side) side == lowS || (Side) side == hiS )) 
                {
                  cout << "Setting Neumann condition on side " << side << endl;
                  setBoundaryConditionType(k,(Side) side, NEUMANN);
                }
	      else if (fType == xComponentOfVector && ( (Side) side == lowR || (Side) side == hiR ))
		setBoundaryConditionType(k,(Side) side, DIRICHLET);
	      else if (fType == yComponentOfVector && ( (Side) side == lowS || (Side) side == hiS ))
		setBoundaryConditionType(k,(Side) side, DIRICHLET);
	      else if (fType == yComponentOfVector && ( (Side) side == lowR || (Side) side == hiR ))
		setBoundaryConditionType(k,(Side) side, NEUMANN);
	    }
	  else if (myGridM->getBoundaryType(k,(Side) side) == INFLOW) // inflow
	    {
	      if (fType == scalar)
		setBoundaryConditionType(k,(Side) side, DIRICHLET);
	      else if (fType == xComponentOfVector || fType == yComponentOfVector)
		setBoundaryConditionType(k,(Side) side, DIRICHLET);
	    }
	  else if (myGridM->getBoundaryType(k,(Side) side) == OUTFLOW) // outflow
	    {
	      if (fType == scalar)
		setBoundaryConditionType(k,(Side) side, ROBIN);
	      else if (fType == xComponentOfVector || fType == yComponentOfVector)
		setBoundaryConditionType(k,(Side) side, EXTRAPOLATION);
	    }
	  else if (myGridM->getBoundaryType(k,(Side) side) == OUTFLOW_NEUM_P) // outflow with Neumann pressure
	    {
	      if (fType == scalar)
		setBoundaryConditionType(k,(Side) side, NEUMANN);
	      else if (fType == xComponentOfVector || fType == yComponentOfVector)
		setBoundaryConditionType(k,(Side) side, EXTRAPOLATION);
	    }
	  else if (myGridM->getBoundaryType(k,(Side) side) == DIRICHLET) // twilight-zone
	    {
	      setBoundaryConditionType(k,(Side) side, DIRICHLET);
	    }
	}
    }
}

void gridFunction::setFlagValues()
{
  int i, j, grid;

  for (grid=0; grid<myGridM->nrGrids(); grid++)
    {
      //      flagValuesM[grid].display();

      int low_i = flagValuesM[grid].getBase(0);
      int hi_i = flagValuesM[grid].getBound(0);
      int low_j = flagValuesM[grid].getBase(1);
      int hi_j = flagValuesM[grid].getBound(1);
      
      for (j=low_j+1; j<=hi_j-1; j++)
	{
	  if (flagValuesM[grid](low_i+1,j) != 0) // !=0 means it is not a hole point
	    {
	      if (getBoundaryType(grid,(Side) 0) == NEUMANN) // Neumann condition
		{
		  flagValuesM[grid](low_i,j) = grid + 2 + 1000;
		}
	      else if (getBoundaryType(grid,(Side) 0) == DIRICHLET) // Dirichlet condition
		{
		  flagValuesM[grid](low_i+1,j) = grid + 1 + 1000;
		  flagValuesM[grid](low_i,j) = grid + 1000;; // extrapolation from inner points
		}
	      else if (getBoundaryType(grid,(Side) 0) == PERIODIC) // periodic boundary
		{
		  flagValuesM[grid](low_i,j) = grid + 3 + 1000;
		}
	      else if (getBoundaryType(grid,(Side) 0) == EXTRAPOLATION) // extrapolation boundary
		{
		  flagValuesM[grid](low_i,j) = grid + 1000;
		}
	      else if (getBoundaryType(grid,(Side) 0) == ROBIN) // robin boundary
		{
		  flagValuesM[grid](low_i,j) = grid + 4 +1000;
		}
	    }
	  if (flagValuesM[grid](hi_i-1,j) != 0) // !=0 means it is not a hole point
	    {
	      if (getBoundaryType(grid,(Side) 1) == NEUMANN) // Neumann condition
		{
		  flagValuesM[grid](hi_i,j) = grid + 2 + 2000;
		}
	      else if (getBoundaryType(grid,(Side) 1) == DIRICHLET) // Dirichlet condition
		{
		  flagValuesM[grid](hi_i-1,j) = grid + 1 + 2000;
		  flagValuesM[grid](hi_i,j) = grid + 2000;; // extrapolation from inner points
		}
	      else if (getBoundaryType(grid,(Side) 1) == PERIODIC) // periodic boundary
		{
		  flagValuesM[grid](hi_i,j) = grid + 3 + 2000;
		}
	      else if (getBoundaryType(grid,(Side) 1) == EXTRAPOLATION) // extrapolation boundary
		{
		  flagValuesM[grid](hi_i,j) = grid + 2000;
		}
	      else if (getBoundaryType(grid,(Side) 1) == ROBIN) // robin boundary
		{
		  flagValuesM[grid](hi_i,j) = grid + 4 + 2000;
		}
	    }
	}
      
      for (i=low_i+1; i<=hi_i-1; i++)
	{
	  if (flagValuesM[grid](i,low_j+1) != 0) // !=0 means it is not a hole point
	    {
	      if (getBoundaryType(grid,(Side) 2) == NEUMANN) // Neumann condition
		{
		  flagValuesM[grid](i,low_j) = grid + 2 + 3000;
		}
	      else if (getBoundaryType(grid,(Side) 2) == DIRICHLET) // Dirichlet condition
		{
		  flagValuesM[grid](i,low_j+1) = grid + 1 + 3000;
		  flagValuesM[grid](i,low_j) = grid + 3000;; // extrapolation from inner points
		}
	      else if (getBoundaryType(grid,(Side) 2) == PERIODIC) // periodic boundary
		{
		  flagValuesM[grid](i,low_j) = grid + 3 + 3000;
		}
	      else if (getBoundaryType(grid,(Side) 2) == EXTRAPOLATION) // extrapolation boundary
		{
		  flagValuesM[grid](i,low_j) = grid + 3000;
		}
	      else if (getBoundaryType(grid,(Side) 2) == ROBIN) // robin boundary
		{
		  flagValuesM[grid](i,low_j) = grid + 4 + 3000;
		}
	    }
	  if (flagValuesM[grid](i,hi_j-1) != 0) // >0 means it is not a hole point
	    {
	      if (getBoundaryType(grid,(Side) 3) == NEUMANN) // Neumann condition
		{
		  flagValuesM[grid](i,hi_j) = grid + 2 + 4000;
		}
	      else if (getBoundaryType(grid,(Side) 3) == DIRICHLET) // Dirichlet condition
		{
		  flagValuesM[grid](i,hi_j-1) = grid + 1 + 4000;
		  flagValuesM[grid](i,hi_j) = grid + 4000;; // extrapolation from inner points
		}
	      else if (getBoundaryType(grid,(Side) 3) == PERIODIC) // periodic boundary
		{
		  flagValuesM[grid](i,hi_j) = grid + 3 + 4000;
		}
	      else if (getBoundaryType(grid,(Side) 3) == EXTRAPOLATION) // extrapolation boundary
		{
		  flagValuesM[grid](i,hi_j) = grid + 4000;
		}
	      else if (getBoundaryType(grid,(Side) 3) == ROBIN) // robin boundary
		{
		  flagValuesM[grid](i,hi_j) = grid + 4 + 4000;
		}
	    }
	}
      
      // (almost) lastly, mark the corners
      if (getBoundaryType(grid,(Side) 0) == PERIODIC)
	flagValuesM[grid](low_i,low_j) = grid + 3 + 1000;
      else if (getBoundaryType(grid,(Side) 2) == PERIODIC)
	flagValuesM[grid](low_i,low_j) = grid + 3 + 3000;
      else
	flagValuesM[grid](low_i,low_j) = 10000;
      
      if (getBoundaryType(grid,(Side) 0) == PERIODIC)
	flagValuesM[grid](low_i,hi_j) = grid + 3 + 1000;
      else if (getBoundaryType(grid,(Side) 3) == PERIODIC)
	flagValuesM[grid](low_i,hi_j) = grid + 3 + 4000;
      else
	flagValuesM[grid](low_i,hi_j) = 20000;
      
      if (getBoundaryType(grid,(Side) 2) == PERIODIC)
	flagValuesM[grid](hi_i,low_j) = grid + 3 + 3000;
      else if (getBoundaryType(grid,(Side) 1) == PERIODIC)
	flagValuesM[grid](hi_i,low_j) = grid + 3 + 2000;
      else
	flagValuesM[grid](hi_i,low_j) = 30000;
      
      if (getBoundaryType(grid,(Side) 1) == PERIODIC)
	flagValuesM[grid](hi_i,hi_j) = grid + 3 + 2000;
      else if (getBoundaryType(grid,(Side) 3) == PERIODIC)
	flagValuesM[grid](hi_i,hi_j) = grid + 3 + 4000;
      else
	flagValuesM[grid](hi_i,hi_j) = 40000;
      
      //----------------------------------------
      // Fix the no-slip (Dirichlet) extremal
      // points
      //----------------------------------------
      if (getBoundaryType(grid,(Side) 0) == DIRICHLET)
	{
	  flagValuesM[grid](low_i+1,low_j) = grid + 1 + 1000;
	  flagValuesM[grid](low_i+1,hi_j) = grid + 1 + 1000;
	}
      
      if (getBoundaryType(grid,(Side) 1) == DIRICHLET)
	{
	  flagValuesM[grid](hi_i-1,low_j) = grid + 1 + 2000;
	  flagValuesM[grid](hi_i-1,hi_j) = grid + 1 + 2000;
	}
      
      if (getBoundaryType(grid,(Side) 2) == DIRICHLET)
	{
	  flagValuesM[grid](low_i,low_j+1) = grid + 1 + 3000;
	  flagValuesM[grid](hi_i,low_j+1) = grid + 1 + 3000;
	}
      
      if (getBoundaryType(grid,(Side) 3) == DIRICHLET)
	{
	  flagValuesM[grid](low_i,hi_j-1) = grid + 1 + 4000;
	  flagValuesM[grid](hi_i,hi_j-1) = grid + 1 + 4000;
	}
      //     flagValuesM[grid].display();exit(1);
    }

  /*  for (int cuG=0;cuG<nmbrOfGrids; cuG++)
      {
      FDmask[cuG].redim(rdimM[cuG],sdimM[cuG]);
      intSerialArray *localFlags;
      localFlags = flagValuesM[cuG].getSerialArrayPointer();
      doubleSerialArray *localFDMask;
      localFDMask = FDmask[cuG].getSerialArrayPointer();
      for (int jy=(*localFlags).getBase(1); jy<=(*localFlags).getBound(1); jy++)
      {
      for (int ix=(*localFlags).getBase(0); ix<=(*localFlags).getBound(0); ix++)
      {
	if (((*localFlags)(ix,jy)) == (cuG + 1))
	  (*localFDMask)(ix,jy) = 1;
	else
	  (*localFDMask)(ix,jy) = 0;
	  }
	  }
	  }
  */

}

Range gridFunction::getDiscrBounds(int direction, int k) const
{
  int low = 0, hi = 0;
  
  if (direction == 0)
    {
      if (bcsM[k].lowR == NEUMANN || bcsM[k].lowR == ROBIN || bcsM[k].lowR == PERIODIC || bcsM[k].lowR == EXTRAPOLATION) 
	{
	  low = 1;
	}
      else
	{
	  low = 2;
	}
      if (bcsM[k].hiR == NEUMANN || bcsM[k].hiR == ROBIN || bcsM[k].hiR == PERIODIC || bcsM[k].hiR == EXTRAPOLATION) 
	{
	  hi = myGridM->rdimM[k] - 2;
	}
      else
	{
	  hi = myGridM->rdimM[k] - 3;
	}
    }
  else  if (direction==1)
    {
      if (bcsM[k].lowS == NEUMANN || bcsM[k].lowS == ROBIN || bcsM[k].lowS == PERIODIC || bcsM[k].lowS == EXTRAPOLATION) 
	{
	  low = 1;
	}
      else
	{
	  low = 2;
	}
      if (bcsM[k].hiS == NEUMANN || bcsM[k].hiS == ROBIN || bcsM[k].hiS == PERIODIC || bcsM[k].hiS == EXTRAPOLATION) 
	{
	  hi = myGridM->sdimM[k] - 2;
	}
      else
	{
	  hi = myGridM->sdimM[k] - 3;
	}
    }
  Range retRange(low,hi);
  
  return retRange;
}

void gridFunction::extrapolateToBoundaries()
{
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ixb(2,(myGridM->rDim(k))-4);
      Index low_r(1), hi_r((myGridM->rDim(k))-2);
      Index iyb(2,(myGridM->sDim(k))-4);
      Index low_s(1), hi_s((myGridM->sDim(k))-2);
      
      if (getBoundaryType(k,(Side) 0) != PERIODIC)
	{
	  fieldValuesM[k](low_r,iyb) = 3.*fieldValuesM[k](low_r+1,iyb) - 3.*fieldValuesM[k](low_r+2,iyb) + 1.*fieldValuesM[k](low_r+3,iyb);
	}
      
      if (getBoundaryType(k,(Side) 1) != PERIODIC)
	{
	  fieldValuesM[k](hi_r,iyb) = 3.*fieldValuesM[k](hi_r-1,iyb) - 3.*fieldValuesM[k](hi_r-2,iyb) + 1.*fieldValuesM[k](hi_r-3,iyb);
	}
      
      if (getBoundaryType(k,(Side) 2) != PERIODIC)
	{
	  fieldValuesM[k](ixb,low_s) = 3.*fieldValuesM[k](ixb,low_s+1) - 3.*fieldValuesM[k](ixb,low_s+2) + 1.*fieldValuesM[k](ixb,low_s+3);
	}
      
      if (getBoundaryType(k,(Side) 3) != PERIODIC)
	{
	  fieldValuesM[k](ixb,hi_s) = 3.*fieldValuesM[k](ixb,hi_s-1) - 3.*fieldValuesM[k](ixb,hi_s-2) + 1.*fieldValuesM[k](ixb,hi_s-3);
	}
      
      if (getBoundaryType(k,(Side) 0) == PERIODIC)
	fieldValuesM[k](1,1) = 3.*fieldValuesM[k](1,2) - 3*fieldValuesM[k](1,3) + 1*fieldValuesM[k](1,4);
      else if (getBoundaryType(k,(Side) 2) == PERIODIC)
	fieldValuesM[k](1,1) = 3.*fieldValuesM[k](2,1) - 3*fieldValuesM[k](3,1) + 1*fieldValuesM[k](4,1);
      else
	fieldValuesM[k](1,1) = 3.*fieldValuesM[k](2,2) - 3*fieldValuesM[k](3,3) + 1*fieldValuesM[k](4,4);
      
      if (getBoundaryType(k,(Side) 1) == PERIODIC)
	fieldValuesM[k](myGridM->rdimM[k]-2,1) = 3.*fieldValuesM[k](myGridM->rdimM[k]-2,2) - 3*fieldValuesM[k](myGridM->rdimM[k]-2,3) + 1*fieldValuesM[k](myGridM->rdimM[k]-2,4);
      else if (getBoundaryType(k,(Side) 2) == PERIODIC)
	fieldValuesM[k](myGridM->rdimM[k]-2,1) = 3.*fieldValuesM[k](myGridM->rdimM[k]-3,1) - 3*fieldValuesM[k](myGridM->rdimM[k]-4,1) + 1*fieldValuesM[k](myGridM->rdimM[k]-5,1);
      else
	fieldValuesM[k](myGridM->rdimM[k]-2,1) = 3.*fieldValuesM[k](myGridM->rdimM[k]-3,2) - 3*fieldValuesM[k](myGridM->rdimM[k]-4,3) + 1*fieldValuesM[k](myGridM->rdimM[k]-5,4);
      
      if (getBoundaryType(k,(Side) 0) == PERIODIC)
	fieldValuesM[k](1,myGridM->sdimM[k]-2) = 3.*fieldValuesM[k](1,myGridM->sdimM[k]-3) - 3*fieldValuesM[k](1,myGridM->sdimM[k]-4) + 1*fieldValuesM[k](1,myGridM->sdimM[k]-5);
      else if (getBoundaryType(k,(Side) 3) == PERIODIC)
	fieldValuesM[k](1,myGridM->sdimM[k]-2) = 3.*fieldValuesM[k](2,myGridM->sdimM[k]-2) - 3*fieldValuesM[k](3,myGridM->sdimM[k]-2) + 1*fieldValuesM[k](4,myGridM->sdimM[k]-2);
      else
	fieldValuesM[k](1,myGridM->sdimM[k]-2) = 3.*fieldValuesM[k](2,myGridM->sdimM[k]-3) - 3*fieldValuesM[k](3,myGridM->sdimM[k]-4) + 1*fieldValuesM[k](4,myGridM->sdimM[k]-5);
      
      if (getBoundaryType(k,(Side) 1) == PERIODIC)
	fieldValuesM[k](myGridM->rdimM[k]-2,myGridM->sdimM[k]-2) = 3.*fieldValuesM[k](myGridM->rdimM[k]-2,myGridM->sdimM[k]-3) - 3*fieldValuesM[k](myGridM->rdimM[k]-2,myGridM->sdimM[k]-4) + 1*fieldValuesM[k](myGridM->rdimM[k]-2,myGridM->sdimM[k]-5);
      else if (getBoundaryType(k,(Side) 3) == PERIODIC)
	fieldValuesM[k](myGridM->rdimM[k]-2,myGridM->sdimM[k]-2) = 3.*fieldValuesM[k](myGridM->rdimM[k]-3,myGridM->sdimM[k]-2) - 3*fieldValuesM[k](myGridM->rdimM[k]-4,myGridM->sdimM[k]-2) + 1*fieldValuesM[k](myGridM->rdimM[k]-5,myGridM->sdimM[k]-2);
      else
	fieldValuesM[k](myGridM->rdimM[k]-2,myGridM->sdimM[k]-2) = 3.*fieldValuesM[k](myGridM->rdimM[k]-3,myGridM->sdimM[k]-3) - 3*fieldValuesM[k](myGridM->rdimM[k]-4,myGridM->sdimM[k]-4) + 1*fieldValuesM[k](myGridM->rdimM[k]-5,myGridM->sdimM[k]-5);
    }
}

void gridFunction::setDivergenceToZeroAtBoundaries(gridFunction & otherVelocityComponent, functionType fType)
#define ys(x,y) (myGridM->ysM[k](x,y))
#define yr(x,y) (myGridM->yrM[k](x,y))
#define xs(x,y) (myGridM->xsM[k](x,y))
#define xr(x,y) (myGridM->xrM[k](x,y))
#define r_step (myGridM->rStep(k))
#define s_step (myGridM->sStep(k))
{
  otherVelocityComponent.updateGhostBoundaries();
  for (int k=0; k<myGridM->nrGrids(); k++)
    {
      Index ixb(1,(myGridM->rDim(k))-2);
      Index low_r(1), hi_r((myGridM->rDim(k))-2);
      Index iyb(1,(myGridM->sDim(k))-2);
      Index low_s(1), hi_s((myGridM->sDim(k))-2);
      //          std::cout << "GRIDTYPE IS " << myGridM->gType(k) << endl << endl;
      //----------------------------------------
      // Check gridtype to see how advanced the
      // divergence operator is
      //----------------------------------------
      if (myGridM->gType(k) == 1)
	//----------------------------------------
	// This is a regular (non-stretched)
	// Cartesian grid
	//----------------------------------------
	{
	  if (myGridM->getBoundaryType(k,(Side) 0) != SLIP && myGridM->getBoundaryType(k,(Side) 0) != PERIODIC && myGridM->getBoundaryType(k,(Side) 0) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 0) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 0) != OUTFLOW_NEUM_P && fType == xComponentOfVector)
	    {
	      fieldValuesM[k](low_r-1,iyb) = fieldValuesM[k](low_r+1,iyb) +
		(myGridM->xM[k](low_r+1,iyb)-myGridM->xM[k](low_r-1,iyb))*
		otherVelocityComponent.y(k,low_r,iyb);
	    }
	  if (myGridM->getBoundaryType(k,(Side) 1) != SLIP && myGridM->getBoundaryType(k,(Side) 1) != PERIODIC && myGridM->getBoundaryType(k,(Side) 1) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 1) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 1) != OUTFLOW_NEUM_P && fType == xComponentOfVector)
	    {
	      fieldValuesM[k](hi_r+1,iyb) = fieldValuesM[k](hi_r-1,iyb) +
		(myGridM->xM[k](hi_r-1,iyb)-myGridM->xM[k](hi_r+1,iyb))*
		otherVelocityComponent.y(k,hi_r,iyb);
	    }
	  if (myGridM->getBoundaryType(k,(Side) 2) != SLIP && myGridM->getBoundaryType(k,(Side) 2) != PERIODIC && myGridM->getBoundaryType(k,(Side) 2) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 2) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 2) != OUTFLOW_NEUM_P && fType == yComponentOfVector)
	    {
	      fieldValuesM[k](ixb,low_s-1) = fieldValuesM[k](ixb,low_s+1) +
		(myGridM->yM[k](ixb,low_s+1)-myGridM->yM[k](ixb,low_s-1))*
		otherVelocityComponent.x(k,ixb,low_s);
	    }
	  if (myGridM->getBoundaryType(k,(Side) 3) != SLIP && myGridM->getBoundaryType(k,(Side) 3) != PERIODIC && myGridM->getBoundaryType(k,(Side) 3) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 3) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 3) != OUTFLOW_NEUM_P && fType == yComponentOfVector)
	    {
	      fieldValuesM[k](ixb,hi_s+1) = fieldValuesM[k](ixb,hi_s-1) +
		(myGridM->yM[k](ixb,hi_s-1)-myGridM->yM[k](ixb,hi_s+1))*
		otherVelocityComponent.x(k,ixb,hi_s);
	    }
	}
      else if (myGridM->gType(k) == 2)
	//----------------------------------------
	// Cartesian grid
	//----------------------------------------
	{
	  if (myGridM->getBoundaryType(k,(Side) 0) != SLIP && myGridM->getBoundaryType(k,(Side) 0) != PERIODIC && myGridM->getBoundaryType(k,(Side) 0) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 0) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 0) != OUTFLOW_NEUM_P && fType == xComponentOfVector)
	    {
	      fieldValuesM[k](low_r-1,iyb) = fieldValuesM[k](low_r+1,iyb) +
		2.*r_step*xr(low_r,iyb)/ys(low_r,iyb)*
		( otherVelocityComponent(k)(low_r,iyb+1) - otherVelocityComponent(k)(low_r,iyb-1) ) / (2.*s_step);
	    }
	  if (myGridM->getBoundaryType(k,(Side) 1) != SLIP && myGridM->getBoundaryType(k,(Side) 1) != PERIODIC && myGridM->getBoundaryType(k,(Side) 1) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 1) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 1) != OUTFLOW_NEUM_P && fType == xComponentOfVector)
	    {
	      fieldValuesM[k](hi_r+1,iyb) = fieldValuesM[k](hi_r-1,iyb) -
		2.*r_step*xr(hi_r,iyb)/ys(hi_r,iyb)*
		( otherVelocityComponent(k)(hi_r,iyb+1) - otherVelocityComponent(k)(hi_r,iyb-1) ) / (2.*s_step);
	    }
	  if (myGridM->getBoundaryType(k,(Side) 2) != SLIP && myGridM->getBoundaryType(k,(Side) 2) != PERIODIC && myGridM->getBoundaryType(k,(Side) 2) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 2) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 2) != OUTFLOW_NEUM_P && fType == yComponentOfVector)
	    {
	      fieldValuesM[k](ixb,low_s-1) = fieldValuesM[k](ixb,low_s+1) +
		2.*s_step*ys(ixb,low_s)/xr(ixb,low_s)*
		( otherVelocityComponent(k)(ixb+1,low_s) - otherVelocityComponent(k)(ixb-1,low_s) ) / (2.*r_step);
	    }
	  if (myGridM->getBoundaryType(k,(Side) 3) != SLIP && myGridM->getBoundaryType(k,(Side) 3) != PERIODIC && myGridM->getBoundaryType(k,(Side) 3) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 3) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 3) != OUTFLOW_NEUM_P && fType == yComponentOfVector)
	    {
	      fieldValuesM[k](ixb,hi_s+1) = fieldValuesM[k](ixb,hi_s-1) -
		2.*s_step*ys(ixb,hi_s)/xr(ixb,hi_s)*
		( otherVelocityComponent(k)(ixb+1,hi_s) - otherVelocityComponent(k)(ixb-1,hi_s) ) / (2.*r_step);
	    }
	}
      else if (myGridM->gType(k) == 3)
	//----------------------------------------
	// Curvilinear grid
	//----------------------------------------
	{
	  if (myGridM->getBoundaryType(k,(Side) 0) != SLIP && myGridM->getBoundaryType(k,(Side) 0) != PERIODIC && myGridM->getBoundaryType(k,(Side) 0) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 0) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 0) != OUTFLOW_NEUM_P)
	    {
	      if (fType == xComponentOfVector)
		{
		  fieldValuesM[k](low_r-1,iyb) = 
		    (r_step*((-otherVelocityComponent(k)(low_r,iyb-1) + otherVelocityComponent(k)(low_r,iyb+1))*xr(low_r,iyb) + 
			     (fieldValuesM[k](low_r,iyb-1) - fieldValuesM[k](low_r,iyb+1))*yr(low_r,iyb))*ys(low_r,iyb) + 
		     s_step*(xs(low_r,iyb)*((3.*fieldValuesM[k](low_r,iyb) - 3.*fieldValuesM[k](low_r+1,iyb) + fieldValuesM[k](low_r+2,iyb))*xs(low_r,iyb) + 
					    otherVelocityComponent(k)(low_r+2,iyb)*ys(low_r,iyb)) + (3.*otherVelocityComponent(k)(low_r,iyb) - 4.*otherVelocityComponent(k)(low_r+1,iyb))*xs(low_r,iyb)*ys(low_r,iyb) + 
			     fieldValuesM[k](low_r+1,iyb)*pow(ys(low_r,iyb),2.)))/
		    (s_step*(pow(xs(low_r,iyb),2.) + pow(ys(low_r,iyb),2.)));
		}
	      else if (fType == yComponentOfVector)
		{
		  fieldValuesM[k](low_r-1,iyb) = 
		    (r_step*xs(low_r,iyb)*((fieldValuesM[k](low_r,iyb-1) - fieldValuesM[k](low_r,iyb+1))*xr(low_r,iyb) + 
					   (-otherVelocityComponent(k)(low_r,iyb-1) + otherVelocityComponent(k)(low_r,iyb+1))*yr(low_r,iyb)) + 
		     s_step*(ys(low_r,iyb)*((3.*otherVelocityComponent(k)(low_r,iyb) - 4.*otherVelocityComponent(k)(low_r+1,iyb) + otherVelocityComponent(k)(low_r+2,iyb))*xs(low_r,iyb) + 
					    fieldValuesM[k](low_r+2,iyb)*ys(low_r,iyb) + 3*fieldValuesM[k](low_r,iyb)*ys(low_r,iyb)) + 
			     fieldValuesM[k](low_r+1,iyb)*(pow(xs(low_r,iyb),2.) - 3*pow(ys(low_r,iyb),2.))))/
		    (s_step*(pow(xs(low_r,iyb),2.) + pow(ys(low_r,iyb),2.)));
		}
	    }
	  if (myGridM->getBoundaryType(k,(Side) 1) != SLIP && myGridM->getBoundaryType(k,(Side) 1) != PERIODIC && myGridM->getBoundaryType(k,(Side) 1) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 1) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 1) != OUTFLOW_NEUM_P)
	    {
	      if (fType == xComponentOfVector)
		{
		  fieldValuesM[k](hi_r+1,iyb) =
		    (r_step*(otherVelocityComponent(k)(hi_r,iyb-1)*xr(hi_r,iyb) - otherVelocityComponent(k)(hi_r,iyb+1)*xr(hi_r,iyb) + 
			     (-fieldValuesM[k](hi_r,iyb-1) + fieldValuesM[k](hi_r,iyb+1))*yr(hi_r,iyb))*ys(hi_r,iyb) + 
		     s_step*(fieldValuesM[k](hi_r-2,iyb)*pow(xs(hi_r,iyb),2.) + 
			     xs(hi_r,iyb)*(3.*fieldValuesM[k](hi_r,iyb)*xs(hi_r,iyb) + 
					   (otherVelocityComponent(k)(hi_r-2,iyb) - 4.*otherVelocityComponent(k)(hi_r-1,iyb) + 3.*otherVelocityComponent(k)(hi_r,iyb))*ys(hi_r,iyb)) + 
			     fieldValuesM[k](hi_r-1,iyb)*(-3.*pow(xs(hi_r,iyb),2.) + pow(ys(hi_r,iyb),2.))))/
		    (s_step*(pow(xs(hi_r,iyb),2.) + pow(ys(hi_r,iyb),2.)));
		}
	      else if (fType == yComponentOfVector)
		{
		  fieldValuesM[k](hi_r+1,iyb) =
		    (r_step*xs(hi_r,iyb)*(-(fieldValuesM[k](hi_r,iyb-1)*xr(hi_r,iyb)) + fieldValuesM[k](hi_r,iyb+1)*xr(hi_r,iyb) + 
					  (otherVelocityComponent(k)(hi_r,iyb-1) - otherVelocityComponent(k)(hi_r,iyb+1))*yr(hi_r,iyb)) + 
		     s_step*(ys(hi_r,iyb)*(otherVelocityComponent(k)(hi_r-2,iyb)*xs(hi_r,iyb) - 4.*otherVelocityComponent(k)(hi_r-1,iyb)*xs(hi_r,iyb) + 
					   3.*otherVelocityComponent(k)(hi_r,iyb)*xs(hi_r,iyb) + fieldValuesM[k](hi_r-2,iyb)*ys(hi_r,iyb) + 3.*fieldValuesM[k](hi_r,iyb)*ys(hi_r,iyb)) + 
			     fieldValuesM[k](hi_r-1,iyb)*(pow(xs(hi_r,iyb),2.) - 3.*pow(ys(hi_r,iyb),2.))))/
		    (s_step*(pow(xs(hi_r,iyb),2.) + pow(ys(hi_r,iyb),2.)));
		}
	    }
	  if (myGridM->getBoundaryType(k,(Side) 2) != SLIP && myGridM->getBoundaryType(k,(Side) 2) != PERIODIC && myGridM->getBoundaryType(k,(Side) 2) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 2) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 2) != OUTFLOW_NEUM_P)
	    {
	      if (fType == xComponentOfVector)
		{
		  fieldValuesM[k](ixb,low_s-1) = 
		    (r_step*((3.*fieldValuesM[k](ixb,low_s) - 3.*fieldValuesM[k](ixb,low_s+1) + fieldValuesM[k](ixb,low_s+2))*
			     pow(xr(ixb,low_s),2.) + 
			     (3.*otherVelocityComponent(k)(ixb,low_s) - 4.*otherVelocityComponent(k)(ixb,low_s+1) + otherVelocityComponent(k)(ixb,low_s+2))*xr(ixb,low_s)*yr(ixb,low_s) + 
			     fieldValuesM[k](ixb,low_s+1)*pow(yr(ixb,low_s),2.)) + 
		     s_step*yr(ixb,low_s)*((-otherVelocityComponent(k)(ixb-1,low_s) + otherVelocityComponent(k)(ixb+1,low_s))*xs(ixb,low_s) + 
					   (fieldValuesM[k](ixb-1,low_s) - fieldValuesM[k](ixb+1,low_s))*ys(ixb,low_s)))/
		    (r_step*(pow(xr(ixb,low_s),2.) + pow(yr(ixb,low_s),2.)));
		}
	      else if (fType == yComponentOfVector)
		{
		  fieldValuesM[k](ixb,low_s-1) =
		    (r_step*(fieldValuesM[k](ixb,low_s+1)*pow(xr(ixb,low_s),2.) + 
			     (3.*otherVelocityComponent(k)(ixb,low_s) - 4.*otherVelocityComponent(k)(ixb,low_s+1) + otherVelocityComponent(k)(ixb,low_s+2))*xr(ixb,low_s)*yr(ixb,low_s) + 
			     (3.*fieldValuesM[k](ixb,low_s) - 3.*fieldValuesM[k](ixb,low_s+1) + fieldValuesM[k](ixb,low_s+2))*pow(yr(ixb,low_s),2.)) + 
		     s_step*xr(ixb,low_s)*((fieldValuesM[k](ixb-1,low_s) - fieldValuesM[k](ixb+1,low_s))*xs(ixb,low_s) + 
					   (-otherVelocityComponent(k)(ixb-1,low_s) + otherVelocityComponent(k)(ixb+1,low_s))*ys(ixb,low_s)))/
		    (r_step*(pow(xr(ixb,low_s),2.) + pow(yr(ixb,low_s),2.)));
		}
	    }
	  if (myGridM->getBoundaryType(k,(Side) 3) != SLIP && myGridM->getBoundaryType(k,(Side) 3) != PERIODIC && myGridM->getBoundaryType(k,(Side) 3) != INTERPOLATION && myGridM->getBoundaryType(k,(Side) 3) != OUTFLOW && myGridM->getBoundaryType(k,(Side) 3) != OUTFLOW_NEUM_P)
	    {
	      if (fType == xComponentOfVector)
		{
		  fieldValuesM[k](ixb,hi_s+1) =
		    (r_step*((fieldValuesM[k](ixb,hi_s-2) - 3.*fieldValuesM[k](ixb,hi_s-1) + 3.*fieldValuesM[k](ixb,hi_s))*
			     pow(xr(ixb,hi_s),2.) + 
			     (otherVelocityComponent(k)(ixb,hi_s-2) - 4.*otherVelocityComponent(k)(ixb,hi_s-1) + 3.*otherVelocityComponent(k)(ixb,hi_s))*xr(ixb,hi_s)*yr(ixb,hi_s) + 
			     fieldValuesM[k](ixb,hi_s-1)*pow(yr(ixb,hi_s),2.)) + 
		     s_step*yr(ixb,hi_s)*((otherVelocityComponent(k)(ixb-1,hi_s) - otherVelocityComponent(k)(ixb+1,hi_s))*xs(ixb,hi_s) + 
					  (-fieldValuesM[k](ixb-1,hi_s) + fieldValuesM[k](ixb+1,hi_s))*ys(ixb,hi_s)))/
		    (r_step*(pow(xr(ixb,hi_s),2.) + pow(yr(ixb,hi_s),2.)));
		}
	      else if (fType == yComponentOfVector)
		{
		  fieldValuesM[k](ixb,hi_s+1) =
		    (r_step*(fieldValuesM[k](ixb,hi_s-1)*pow(xr(ixb,hi_s),2.) + 
			     (otherVelocityComponent(k)(ixb,hi_s-2) - 4.*otherVelocityComponent(k)(ixb,hi_s-1) + 3.*otherVelocityComponent(k)(ixb,hi_s))*xr(ixb,hi_s)*yr(ixb,hi_s) + 
			     (fieldValuesM[k](ixb,hi_s-2) - 3.*fieldValuesM[k](ixb,hi_s-1) + 3.*fieldValuesM[k](ixb,hi_s))*pow(yr(ixb,hi_s),2.)) + 
		     s_step*xr(ixb,hi_s)*((-fieldValuesM[k](ixb-1,hi_s) + fieldValuesM[k](ixb+1,hi_s))*xs(ixb,hi_s) + 
					  (otherVelocityComponent(k)(ixb-1,hi_s) - otherVelocityComponent(k)(ixb+1,hi_s))*ys(ixb,hi_s)))/
		    (r_step*(pow(xr(ixb,hi_s),2.) + pow(yr(ixb,hi_s),2.)));
		}
	    }
	}
    }
    otherVelocityComponent.updateGhostBoundaries();
}
#undef ys
#undef yr
#undef xs
#undef xr
#undef r_step
#undef s_step

const gridFunction operator+(const gridFunction & left, const gridFunction & right)
{
  gridFunction result = left;
  return result += right;
}

const gridFunction operator-(const gridFunction & left, const gridFunction & right)
{
  gridFunction result = left;
  return result -= right;
}

const gridFunction operator*(const gridFunction & left, const gridFunction & right)
{
  gridFunction result = left;
  return result *= right;
}

const gridFunction operator/(const gridFunction & left, const gridFunction & right)
{
  gridFunction result = left;
  return result /= right;
}

const gridFunction operator*(const double left, const gridFunction & right)
{
  gridFunction result = right;
  for (int k=0; k<right.nrGrids(); k++)
    {
      result(k) *= left;
    }
  return result;
}

const gridFunction operator*(const gridFunction & left, const double right)
{
  gridFunction result = left;
  for (int k=0; k<left.nrGrids(); k++)
    {
      result(k) *= right;
    }
  return result;
}

const gridFunction operator/(const gridFunction & left, const double right)
{
  gridFunction result = left;
  for (int k=0; k<left.nrGrids(); k++)
    {
      result(k) /= right;
    }
  return result;
}

double gridFunction::l_inf()
{
  double maxVal=0;
  double gridMax;

  for (int k=0; k<myGridM->nrGrids(); k++) 
    {
      //      where (myGridM->flagValuesM[k] == k+1)
      //	{
      gridMax = max(abs(fieldValuesM[k])*myGridM->maskM[k])/(k+1);
      //      std::cout << "gridmax at grid " << k << " is " << gridMax << endl;
      //	}
      maxVal = max(maxVal, gridMax);
    }

  return maxVal;
}

double gridFunction::l_2()
{
  double totSum=0;
  double gridSum;

  for (int k=0; k<myGridM->nrGrids(); k++) 
    {
      //      where (myGridM->flagValuesM[k] == k+1)
	//	{
	  gridSum = sum(pow(fieldValuesM[k],2.)*myGridM->maskM[k])/(k+1);
	  //	}
      totSum += gridSum;
    }

  totSum = sqrt(totSum);

  return totSum;
}

double gridFunction::l_h()
{
  double totSum=0;
  double gridSum;

  for (int k=0; k<myGridM->nrGrids(); k++) 
    {
      //      where (myGridM->flagValuesM[k] == k+1)
      //	{
	  gridSum = sum(pow(fieldValuesM[k],2.)*myGridM->maskM[k])/(k+1);
	  //	}
      totSum += gridSum*myGridM->rStep(k)*myGridM->sStep(k);
    }

  totSum = sqrt(totSum);

  return totSum;
}

doubleArray GridFunction::zero(const doubleArray& x, 
			       const doubleArray& y, 
			       const double t, 
			       int grid, 
			       int side)
{
  doubleArray retA = sin(x*0);
  return retA;
}
