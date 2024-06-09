#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>

#include "real.h"
#include "nrutil.h"
#include "hdf_stuff.h"

#define NR_END 1
#define FREE_ARG char*

int 
printf(const char *format, ...);

void 
nrerror(char error_text[]){
  printf("Numerical Recipies run-time error...\n");
  printf("%s\n",error_text);
  printf("...now exiting to system...\n");
  exit(1);
}

real *
vector(long nl, long nh){
  real *v;

  v = (real *)malloc( (nh-nl+1+NR_END)*sizeof(real) );
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

void 
free_vector(real *v, long nl, long nh){
  free((FREE_ARG) (v+nl-NR_END));
}

void 
read_vector( int fd, real *x, int low, int high ){
  int i;
  for ( i=low; i<=high; i++)
    read( fd, (char *) &(x[i]), sizeof(real) );
}

void 
write_vector( int fd, real *x, int low, int high ){
  int i;
  for ( i=low; i<=high; i++)
    write( fd, (char *) &(x[i]), sizeof(real) );
}

#ifndef HDF_DUMMY
/* define a macro for writing a standard int/real array of rank 1 into */
/* the vgroup `vgroup_id' */
int 
hput_vector( real x[], int low, int high, char *name, int vgroup_id ){
  const int32 rank=1; 
  int32 dims[1], start[1], edges[1], base[1], sds_id, istat, ref, file_id, sd_id;
 
  if (( file_id = get_file_id() )== -1){ 
    printf("ERROR: There is no open database file!n"); 
    return -1; 
  } 
  sd_id = get_sd_id();
 
  dims[0] = high - low + 1; 
  start[0] = 0; 
  edges[0] = dims[0]; 
  base[0] = low; 

/* create the array */ 
  sds_id = SDcreate(sd_id, name, DFNT_REAL, rank, dims); 
  istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*) &x[low]); 
/* Save array lower bounds */  \
  SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); \
/* Insert the sds into the current Vgroup.  */  
  ref = SDidtoref(sds_id);  
  Vaddtagref(vgroup_id, DFTAG_NDG, ref ); 
/* terminate access to the array */ 
  istat = SDendaccess(sds_id); 
  return TRUE; 
}

real * 
hget_vector( char * name, int *low, int *high, int vgroup_id ){ 
  int n, i, npairs, found=FALSE, status;
  int32 tag, ref, rank, nt, nattrs, sds_index, sds_id, attr_index, file_id, sd_id
    , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS], base;
  char sd_name[MAX_NC_NAME];
  real * x_ = NULL; 
 
  if (( file_id = get_file_id() )== -1){ 
    printf("ERROR: There is no open database file!\n"); 
    return NULL; 
  } 
  sd_id = get_sd_id();

/*  get the total number of tag/reference pairs in the current vgroup */  
  npairs = Vntagrefs(vgroup_id); 
  for( i=0; i<npairs; i++ ){  
/* get tag and ref */  
    status = Vgettagref(vgroup_id, i, &tag, &ref );  
    if( tag == DFTAG_NDG ){ /* this is a Numeric Data group */  
/* get index (which SDS in the file 0,1,2,...)  */  
      sds_index = SDreftoindex(sd_id, ref);  
/* select this SDS and get identifier  */ 
      sds_id = SDselect(sd_id, sds_index); 
      status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); 
/* The name, rank and number type must match */
      if( !strcmp(name, sd_name) && rank == 1 && nt == DFNT_REAL){ 
/* read the base attribute */
	if ((attr_index = SDfindattr(sds_id, "arrayBase")) >= 0){
	  status = SDreadattr(sds_id, attr_index, &base);
/*	  printf("hget_vector: %s: base = %i\n", name, base);*/
	}
	else{
	  printf("Warning: could not find the arrayBase attribute for the SDS `%s'\n",
		 sd_name);
	  base = 1;
	}
	*low = base;

        found=TRUE; 
	for( n=0; n<rank; n++ ){  
	  start[n]=0; 
	  edges[n]=dims[n]; 
	} 
	*high = *low + dims[0] - 1;
	x_ = vector(*low, *high);
 
/* now read in the array data */  
	status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) &x_[*low]); 
        status = SDendaccess(sds_id); 
        break; 
      } 
      status = SDendaccess(sds_id); 
    } 
  } 
  if( !found ){ 
    printf("get: ERROR searching for `%s'\n", name); 
    *low = *high = 0; 
  } 
  return x_; 
}

#else /* end if not HDF_DUMMY */

int 
hput_vector( real x[], int low, int high, char *name, int vgroup_id ){
  return 1;
}

real * 
hget_vector( char * name, int *low, int *high, int vgroup_id ){ 
  return NULL;
}

#endif /* end if HDF_DUMMY */
