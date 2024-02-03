 /*@@
   @file      ImportAHFinderFile.c
   @date      Wed 05 Sep 2001
   @author    Gerd Lanfermann, Thomas Radke
   @desc
              OpenDX Data Import module to read Cactus output data
              from an AHFinder HDF5 file
   @enddesc
   @version   $Header: /cactus/VizTools/OpenDXutils/src/ImportAHFinderFile.c,v 1.9 2003/03/28 13:01:42 tradke Exp $
 @@*/

#include <dx/dx.h>
#include <hdf5.h>
#include <math.h>

#ifndef M_PI
#define M_PI    3.14159265358979323846  /* pi */
#define M_PI_2  1.57079632679489661923  /* pi/2 */
#endif

/* define this to get some debugging messages */
/*#define DEBUG_ME 1*/

/* the shell environment variable to query
   for an optional AHFinder HDF5 input filename */
#define DX_IMPORT_HDF5  "DX_IMPORT_HDF5"

/* number of dataset IDs to allocate at once */
#define H5_BLOCK_INCREMENT  100

/* macro to check the return code of calls to the HDF5 library */
#define CHECK_HDF5_ERROR(fn_call)                                             \
    do                                                                        \
    {                                                                         \
      int error_code = fn_call;                                               \
                                                                              \
                                                                              \
      if (error_code < 0)                                                     \
      {                                                                       \
        DXMessage ("HDF5 call '%s' returned error code %d",                   \
                   #fn_call, error_code);                                     \
      }                                                                       \
    } while (0)


/* structure definitions describing a single dataset ID and containing
   iteration information */
/* FIXME: this assumes that there are 3 horizons per timestep */


typedef struct
{
  double time;
  /*Field field[3];*/
  float **field;
  float **curvature;
  float **center;
  int *npoints;
  int *cpoints;
  int **cdims;
} dataset_id_t;

typedef struct
{
  int nslices;
  int allocated;
  dataset_id_t *ids;
} iterate_info_t;


/* prototypes of routines defined in this source file */
Error m_ImportAHFinderFile (Object *in, Object *out);
static herr_t processDataset (hid_t file, const char *objectname, void *arg);
static int compare_fn (const void *_a, const void *_b);
Object *AlmToThetaPhiField(float *Alm, int Alm_count,
                           float *curv, int *cdims, float *center);

static float Ylm(float *Alm, int Alm_count, float theta, float phi);
static float faculty_div(int n, int m);
static float Y_factor(int l, int m);
static void A_to_lm(int A, int *l, int *m);
static int lm_to_A(int l, int m);


 /*@@
   @routine    m_ImportAHFinderFile
   @date       Wed 05 Sep 2001
   @author     Gerd Lanfermann, Thomas Radke
   @desc
               The main routine of the ImportAHFinderFile module
               Takes a filename and tries to open the identified file
               as a Cactus AHFinder output file in HDF5 file format.
               All datasets in the file are read as Array objects
               into a Field object which is then passed to connected
               modules.

               Note that the routine name prefix is necessary
               in order to link against the DX library routines.
   @enddesc

   @var        in
   @vdesc      pointer to array of input ports; a string containing the
               filename is expected on in[0]
   @vtype      Object *
   @vio        in
   @endvar
   @var        out
   @vdesc      pointer to array of output ports; the data field
               created from the Cactus file is output on out[0]
   @vtype      Object *
   @vio        in
   @endvar
@@*/
Error m_ImportAHFinderFile (Object *in, Object *out)
{
  Field dx_thetaphi;
  Group dx_allhorizons;
  int is_hdf5, i, h;
  hid_t file;
  char *filename, *filename_env;
  iterate_info_t info;
  H5E_auto_t print_error_fn;
  void      *print_error_fn_arg;

  /* if the file open fails
     also try the environment variable DX_IMPORT_HDF5 if set */
  filename = NULL;
  filename_env = getenv (DX_IMPORT_HDF5);
  DXExtractString (in[0], &filename);

  /* check to see that the file is accessible, and is an hdf file */
  CHECK_HDF5_ERROR (H5Eget_auto (&print_error_fn, &print_error_fn_arg, NULL));
  CHECK_HDF5_ERROR (H5Eset_auto (NULL, NULL, NULL));
  is_hdf5 = H5Fis_hdf5 (filename) > 0;
  if (! is_hdf5 && filename_env)
  {
    DXMessage ("Unable to read file from '%s', trying '%s'",
               filename, filename_env);
    filename = filename_env;
    is_hdf5 = H5Fis_hdf5 (filename) > 0;
  }
  if (! is_hdf5)
  {
    DXMessage ("File '%s' is not a valid AHFinder HDF5 file", filename);
    return (ERROR);
  }
  CHECK_HDF5_ERROR (H5Eset_auto (print_error_fn, print_error_fn_arg, NULL));

  /* initialize iteration info */
  info.allocated = info.nslices = 0;
  info.ids = NULL;

  /* open the file */
  CHECK_HDF5_ERROR (file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT));

  /* iterate over all datasets starting from "/" in the HDF5 file */
  CHECK_HDF5_ERROR (H5Giterate (file, "/", NULL, processDataset, &info));

  /* close the file */
  CHECK_HDF5_ERROR (H5Fclose (file));

  /* now sort the dataset ids by their time attribute */
  qsort (info.ids, info.nslices, sizeof (dataset_id_t), compare_fn);

  /* create a new series */
  out[0] = (Object) DXNewSeries ();
  if (! out[0])
  {
    DXMessage ("Couldn't allocate new time series");
    return (ERROR);
  }

  fprintf(stderr,"Number of timesteps: %d \n",info.nslices);

  /* fill in series members */
  for (i = 0; i < info.nslices; i++)
  {
    dx_allhorizons=DXNewGroup();
    for (h=0;h<3;h++)
    {
      if (info.ids[i].field[h]) {
        dx_thetaphi = (Field) AlmToThetaPhiField(info.ids[i].field[h], info.ids[i].npoints[h],
                                                 info.ids[i].curvature[h],
                                                 info.ids[i].cdims[h],
                                                 info.ids[i].center[h]);

        
        if (!DXSetMember(dx_allhorizons, NULL, (Object)dx_thetaphi)) {
          fprintf(stderr,"DXSetGroup failed: tstep %d horizon# %d",i,h);
          DXDelete((Object)(out[0]));
          DXDelete((Object)(dx_thetaphi));

          return(OK) ;
        }
      }
    }
    DXSetSeriesMember ((Series) out[0], i, info.ids[i].time, (Object)dx_allhorizons);
  }

  return (OK);
}


 /*@@
   @routine    processDataset
   @date       Wed 05 Sep 2001
   @author     Gerd Lanfermann, Thomas Radke
   @desc
               The worker routine which is called by H5Giterate().
               It reads a dataset from a file into an Array object
               and adds it to a Field object.
   @enddesc

   @var        file
   @vdesc      HDF5 object to start the iteration (toplevel hierarchy)
   @vtype      hid_t
   @vio        in
   @endvar
   @var        objectname
   @vdesc      name of the object at the current iteration
   @vtype      const char *
   @vio        in
   @endvar
   @var        arg
   @vdesc      pointer to the Field object where new datasets should be added to
   @vtype      void *
   @vio        in
   @endvar
@@*/
static herr_t processDataset (hid_t file,
                              const char *objectname,
                              void *_iterate_info)
{
  float *array;
  hid_t group, dataset, dataspace, attr, geo;
  hsize_t npoints=0, cpoints, *h5dims;
  H5G_stat_t object_info;
  int rank, idim, *sdims, new_size, horizon;
  float *cntr;

  char datasetname[50];
  H5E_auto_t print_error_fn;
  void      *print_error_fn_arg;
  iterate_info_t *info = (iterate_info_t *) _iterate_info;


#ifdef DEBUG_ME
  fprintf (stderr, "processDataset: processing object '%s'\n", objectname);
#endif

  /* we are interested in time groups only - skip anything else */
  CHECK_HDF5_ERROR (H5Gget_objinfo (file, objectname, 0, &object_info));
  if (object_info.type != H5G_GROUP ||
      strncmp (objectname, "Physical Time ", 14))
  {
    return (0);
  }

#ifdef DEBUG_ME
  fprintf (stderr, "processDataset: getting Physical time \n");
#endif

  /* check for a time attribute */
  CHECK_HDF5_ERROR (group = H5Gopen (file, objectname, H5P_DEFAULT));
  CHECK_HDF5_ERROR (attr = H5Aopen_name (group, "Physical Time"));
  if (attr < 0)
  {
    CHECK_HDF5_ERROR (H5Gclose (group));
    return (0);
  }

#ifdef DEBUG_ME
  fprintf (stderr, "processDataset: extending arrays \n");
#endif

  /* check whether the dataset id array needs to be extended */
  if (info->nslices >= info->allocated)
  {
    info->allocated += H5_BLOCK_INCREMENT;
    new_size = info->allocated * sizeof (dataset_id_t);
    if (info->ids)
    {
      info->ids = realloc (info->ids, new_size);
    }
    else
    {
      info->ids = malloc (new_size);
    }
  }
#ifdef DEBUG_ME
  fprintf (stderr, "processDataset: allocating fields \n");
#endif
  info->ids[info->nslices].field     = malloc(sizeof(float*)*3);
  info->ids[info->nslices].curvature = malloc(sizeof(float*)*3);
  info->ids[info->nslices].center    = malloc(sizeof(float*)*3);
  info->ids[info->nslices].npoints   = malloc(sizeof(int)*3);
  info->ids[info->nslices].cpoints   = malloc(sizeof(int)*3);
  info->ids[info->nslices].cdims     = malloc(sizeof(int*)*3);

  /* read the time attribute */
  CHECK_HDF5_ERROR (H5Aread (attr, H5T_NATIVE_DOUBLE,
                             &info->ids[info->nslices].time));
  CHECK_HDF5_ERROR (H5Aclose (attr));

  fprintf (stderr, "processDataset: getting Physical time %f\n",
           info->ids[info->nslices].time);

  /*********/
  /** YLM **/
  /*********/
  for (horizon = -1; horizon < 3; horizon++)
  {
    /* open the geometry dataset if it exists */
    if (horizon < 0)
    {
      strcpy (datasetname, "Apparent Horizon/Geometry/R");
    }
    else
    {
      sprintf (datasetname, "Apparent Horizon %d/Geometry/R", horizon);
    }
    CHECK_HDF5_ERROR (H5Eget_auto (&print_error_fn, &print_error_fn_arg, NULL));
    CHECK_HDF5_ERROR (H5Eset_auto (NULL, NULL, NULL));
    dataset = H5Dopen (group, datasetname, H5P_DEFAULT);
    CHECK_HDF5_ERROR (H5Eset_auto (print_error_fn, print_error_fn_arg, NULL));
    if (dataset >= 0)
    {

#ifdef DEBUG_ME
      fprintf (stderr, "processDataset YLM: reading dataset '%s' (horizon# %d)\n",
               datasetname, horizon);
#endif

      CHECK_HDF5_ERROR (dataspace = H5Dget_space (dataset));
      CHECK_HDF5_ERROR (rank      = H5Sget_simple_extent_ndims (dataspace));
      CHECK_HDF5_ERROR (npoints   = H5Sget_simple_extent_npoints (dataspace));

      /* allocate some temporary buffers */
      h5dims = malloc (rank * sizeof (hsize_t));
      sdims =  malloc (rank * sizeof (int));
      cntr  =  malloc (sizeof(float)*3);

      if (horizon < 0)
      {
        strcpy (datasetname, "Apparent Horizon/Geometry");
      }
      else
      {
        sprintf (datasetname, "Apparent Horizon %d/Geometry", horizon);
      }
      CHECK_HDF5_ERROR (geo  = H5Gopen (group, datasetname, H5P_DEFAULT));
      CHECK_HDF5_ERROR (attr = H5Aopen_name (geo, "Surface Center"));
      CHECK_HDF5_ERROR (H5Aread (attr, H5T_NATIVE_FLOAT, cntr));
      CHECK_HDF5_ERROR (H5Aclose (attr));
      CHECK_HDF5_ERROR (H5Gclose (geo));

      CHECK_HDF5_ERROR (H5Sget_simple_extent_dims (dataspace, h5dims, NULL));

      for (idim = 0; idim < rank; idim++)
      {
        sdims[idim] = (int) h5dims[idim];
      }

      array = malloc(sizeof(float)*npoints);

      if (! array)
      {
        DXMessage ("Failed to allocate new alm array: (npoints %d)",npoints);
        return (-1);
      }

      CHECK_HDF5_ERROR (H5Dread (dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT, array));
      CHECK_HDF5_ERROR (H5Dclose (dataset));
      CHECK_HDF5_ERROR (H5Sclose (dataspace));

      /* free temporary arrays */
      free (h5dims);
      free (sdims);
    }
    else
    {
#ifdef DEBUG_ME
      fprintf (stderr, "processDataset: no dataset '%s'\n", datasetname);
#endif

      array = NULL;
      cntr  = NULL;
    }

    info->ids[info->nslices].field[horizon < 0 ? 0 : horizon]  = array;
    info->ids[info->nslices].npoints[horizon < 0 ? 0 : horizon]= npoints;
    info->ids[info->nslices].center[horizon < 0 ? 0 : horizon] = cntr;

    /* leave the loop in single-horizon mode */
    if (horizon < 0 && dataset >= 0)
    {
      break;
    }
  }

  /************************/
  /** Gaussian curvature **/
  /************************/
  for (horizon = -1; horizon < 3; horizon++)
  {
    if (!info->ids[info->nslices].field[horizon < 0 ? 0 : horizon]) continue;

    /* open the geometry dataset if it exists */
    if (horizon < 0)
    {
      strcpy (datasetname, "Apparent Horizon/Points/Gaussian curvature");
    }
    else
    {
      sprintf (datasetname, "Apparent Horizon %d/Points/Gaussian curvature", horizon);
    }
    CHECK_HDF5_ERROR (H5Eget_auto (&print_error_fn, &print_error_fn_arg, NULL));
    CHECK_HDF5_ERROR (H5Eset_auto (NULL, NULL, NULL));
    dataset = H5Dopen (group, datasetname, H5P_DEFAULT);
    CHECK_HDF5_ERROR (H5Eset_auto (print_error_fn, print_error_fn_arg, NULL));
    if (dataset >= 0)
    {

      CHECK_HDF5_ERROR (dataspace = H5Dget_space (dataset));
      CHECK_HDF5_ERROR (rank      = H5Sget_simple_extent_ndims (dataspace));
      CHECK_HDF5_ERROR (cpoints   = H5Sget_simple_extent_npoints (dataspace));

      /* allocate some temporary buffers */
      h5dims = malloc (rank * sizeof (hsize_t));
      sdims  = malloc (rank * sizeof (int));

      CHECK_HDF5_ERROR (H5Sget_simple_extent_dims (dataspace, h5dims, NULL));

      for (idim = 0; idim < rank; idim++)
      {
        sdims[idim] = (int) h5dims[idim];
      }

#ifdef DEBUG_ME
      fprintf (stderr, "processDataset Curv: reading dataset '%s' (horizon# %d / %d x %d)\n",
               datasetname, horizon,sdims[0],sdims[1]);
#endif

      array = malloc(sizeof(float)*cpoints);

      if (! array)
      {
        DXMessage ("Failed to allocate new alm array");
        return (-1);
      }

      CHECK_HDF5_ERROR (H5Dread (dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT, array));
      CHECK_HDF5_ERROR (H5Dclose (dataset));
      CHECK_HDF5_ERROR (H5Sclose (dataspace));

      info->ids[info->nslices].cdims[horizon < 0 ? 0 : horizon]      = sdims;
      info->ids[info->nslices].curvature[horizon < 0 ? 0 : horizon]  = array;
      info->ids[info->nslices].cpoints[horizon < 0 ? 0 : horizon]    = cpoints;

      /* free temporary arrays */
      free (h5dims);

    }
    else
    {
#ifdef DEBUG_ME
      fprintf (stderr, "processDataset Curv: no dataset '%s'\n", datasetname);
#endif

      info->ids[info->nslices].cdims[horizon < 0 ? 0 : horizon]  = NULL;
      info->ids[info->nslices].curvature[horizon < 0 ? 0 : horizon] = NULL;
      info->ids[info->nslices].cpoints[horizon < 0 ? 0 : horizon]   = 0;
    }

    /* leave the loop in single-horizon mode */
    if (horizon < 0 && dataset >= 0)
    {
      break;
    }
  }

  /* close the time group */
  CHECK_HDF5_ERROR (H5Gclose (group));

  /* increment dataset counter */
  info->nslices++;

  return (0);
}


/* comparison function used by qsort(3) to sort datasets by time */
static int compare_fn (const void *_a, const void *_b)
{
  const dataset_id_t *a = (const dataset_id_t *) _a;
  const dataset_id_t *b = (const dataset_id_t *) _b;

  return (a->time > b->time);
}


Object *AlmToThetaPhiField(float *Alm, int Alm_count,
                           float *curv, int *cdims, float *center)

{
  Array ThetaPhiArray=NULL;
  Field ThetaPhiField=NULL;
  Point *pos_ptr;

  Array conn_array, pos_array;

  int   counts[2], phi_max,theta_max;
  int   bitant_off;
  float dtheta, dphi, rad;
  float theta, phi, x, y, z;
  float *tpdata_ptr;

  int pc,tc, numelements;
  int idx,idx2;

  phi_max   = cdims[0];
  theta_max = cdims[1];


  dtheta = M_PI / (float)(theta_max-1);
  dphi   = 2.0*M_PI / (float)(phi_max-1);

  numelements = phi_max * 2*theta_max;

  ThetaPhiField = DXNewField();
  if(!(ThetaPhiArray = DXNewArray (TYPE_FLOAT, CATEGORY_REAL, 0)))
    fprintf(stderr,"Failed to create ThetaPhi Array\n");

  if (!(DXAddArrayData(ThetaPhiArray, 0, numelements, NULL)))
    fprintf(stderr,"Failed to allocate memory for ThetaPhi Array: (elements %d)\n",numelements);
  tpdata_ptr = DXGetArrayData(ThetaPhiArray);

  DXSetStringAttribute((Object)ThetaPhiArray, "dep", "positions");
  DXSetComponentValue(ThetaPhiField, "data", (Object)ThetaPhiArray);

  /* create the grid position array */
  pos_array = DXNewArray(TYPE_FLOAT, CATEGORY_REAL, 1, 3);
  DXAddArrayData(pos_array, 0, numelements, NULL);
  pos_ptr = (Point*)DXGetArrayData(pos_array);

  idx=0;
  for (pc=0; pc<phi_max; pc++) {

    phi = pc*dphi;

    for (tc=0; tc<theta_max; tc++) {

      theta = tc*dtheta/2.0;

      idx = pc*phi_max+tc;

      rad = Ylm(Alm, Alm_count, theta, phi);

      x = rad*(sin(theta)*cos(phi))+center[0];
      y = rad*(sin(theta)*sin(phi))+center[1];
      z = rad*(cos(theta))         +center[2];

      pos_ptr[idx]     = DXPt(x,y,z);
      tpdata_ptr[idx]  = curv[idx];
    }
  }
  bitant_off=idx+1;

  for (pc=0; pc<phi_max; pc++) {

    phi = pc*dphi;

    for (tc=0; tc<theta_max; tc++) {

      theta = (M_PI_2)+(tc*dtheta/2.0);

      idx = pc*phi_max+tc;
      idx2= pc*phi_max+(theta_max-1-tc);

      rad = Ylm(Alm, Alm_count, theta, phi);

      x = rad*(sin(theta)*cos(phi))+center[0];
      y = rad*(sin(theta)*sin(phi))+center[1];
      z = rad*(cos(theta))         +center[2];

      pos_ptr[idx+bitant_off]     = DXPt(x,y,z);
      tpdata_ptr[idx+bitant_off]  = curv[idx2];
    }
  }
  counts[0] = 2*theta_max;
  counts[1] = phi_max;
  conn_array = DXMakeGridConnectionsV(2, counts);
  DXSetComponentValue(ThetaPhiField, "connections", (Object)conn_array);

  DXSetComponentValue(ThetaPhiField, "positions", (Object)pos_array);

  DXEndField(ThetaPhiField);

  return((Object*)ThetaPhiField);

}



static float Ylm(float *Alm, int Alm_count, float theta, float phi)
{
  float dU_p = 0, U_p = 0;
  float dU_n = 0, U_n = 0;
  float r_theta, c_theta;
  float r_phi, c_phi;
  float sin_theta, sin_phi;

  float A_p, A_n, Y_norm;
  float dV_p=0, V_p=0, dV_n=0, V_n=0, C_p=0, S_n=0, D=0;

  int l,m,L,M;

  A_to_lm(Alm_count, &L, &M);

  sin_theta = sin(theta);
  sin_phi   = sin(phi);

  if (theta > M_PI)
  {
    r_theta = -1;
    /*  c_theta =  2*sqr(cos(0.5*theta)); */
    c_theta =  2*cos(0.5*theta)*cos(0.5*theta);
  }
  else
  {
    r_theta = 1;
    /*  c_theta = -2*sqr(sin(0.5*theta)); */
    c_theta = -2*sin(0.5*theta)*sin(0.5*theta);
  }

  if (phi > (M_PI_2) && phi < 1.5*(M_PI_2))
  {
    r_phi = -1;
    c_phi =  2*cos(0.5*phi)*cos(0.5*phi);
  }
  else
  {
    r_phi =  1;
    c_phi = -2*sin(0.5*phi)*sin(0.5*phi);
  }

  for(m=L-1; m>=0; m--) {
    dU_p = 0, U_p = 0;
    dU_n = 0, U_n = 0;

    for(l=L-1; l>=m+1; l--)
    {
      A_p = Alm[lm_to_A(l,+m)];
      A_n = Alm[lm_to_A(l,-m)];

      Y_norm = Y_factor(l, m);

      A_p *= Y_norm;
      A_n *= Y_norm;

      dU_p = (A_p+(l+m+1)*r_theta*dU_p + (2*l+1)*c_theta*U_p) / (l-m);
      U_p  = r_theta*U_p + dU_p;

      dU_n = (A_n+(l+m+1)*r_theta*dU_n + (2*l+1)*c_theta*U_n) / (l-m);
      U_n  = r_theta*U_n + dU_n;
    }

    A_p =  Alm[lm_to_A(m,+m)];
    A_n =  Alm[lm_to_A(m,-m)];

    Y_norm = Y_factor(m, m);

    A_p *= Y_norm;
    A_n *= Y_norm;

    U_p = A_p+(2*m+1)*(dU_p*r_theta+c_theta*U_p);
    U_n = A_n+(2*m+1)*(dU_n*r_theta+c_theta*U_n);
    if (m>0) {
      dV_p = U_p-      (2*m+1)*sin_theta*(2*V_p*c_phi + dV_p*r_phi);
       V_p =    -r_phi*(2*m+1)*sin_theta* V_p + dV_p;

      dV_n = U_n-      (2*m+1)*sin_theta*(2*V_n*c_phi + dV_n*r_phi);
       V_n =    -r_phi*(2*m+1)*sin_theta* V_n + dV_n;
    }
  }



  C_p = U_p - sin_theta*( c_phi*V_p + r_phi*dV_p );
  S_n =     - sin_theta*sin_phi*V_n;

  D = S_n + C_p;

  return D;
}

static int lm_to_A(int l, int m)
{
  return l*l + l + m;
}

static void A_to_lm(int A, int *l, int *m)
{
/*   l = int( isqrt(A) ); */
  *l = sqrt((double) A);
  *m = A - (*l)*(*l) - (*l);
}

static float Y_factor(int l, int m)
{
  return (float) sqrt( (m==0? 1 : 2) * (2*l+1)  * faculty_div(l-m , l+m) );
}


static float faculty_div(int n, int m)
{
  float f = 1;
  int i;
  if (n>m) return 1.0/faculty_div(m,n);

  for(i=n+1; i<m+1; i++)
    f *= i;

  return (1./f);
}
