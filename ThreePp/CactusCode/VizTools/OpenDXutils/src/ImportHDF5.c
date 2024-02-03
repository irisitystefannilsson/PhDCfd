 /*@@
   @file      ImportHDF5.c
   @date      Tue 23 Jan 2001
   @author    Frank Herrmann, Gerd Lanfermann, Thomas Radke
   @desc
              Data Import module to read datasets from an HDF5 file
              as unigrid fields into OpenDX

              See the ../README file for general information.
              See the ../COPYING file for copyright and license information.
              See the ../doc/ subdirectory for documentation.
   @enddesc
   @hdate     12 Sep 2006
   @hauthor   Eh Tan, Computational Infrastructure for Geodynamics,
              http://www.geodynamics.org
   @hdesc     Added support for importing vector fields;
              also allow selection by name for the dataset to be imported
   @histor
   @endhistory

   @version   $Header: /cactus/VizTools/OpenDXutils/src/ImportHDF5.c,v 1.30 2006/09/12 08:52:31 tradke Exp $
 @@*/

#include <stdlib.h>
#include <string.h>

#include <dx/dx.h>
#include <hdf5.h>


/*****************************************************************************
 *************************     Macro Definitions   ***************************
 *****************************************************************************/

/* uncomment the following line to get some debug output */
/* #define DEBUG 1 */

/* number of dataset structures to allocate at once */
#define H5_BLOCK_INCREMENT  100

/* macro to check the return code of calls to the HDF5 library */
#define CHECK_ERROR(fn_call)                                                  \
        {                                                                     \
          if ((fn_call) < 0)                                                  \
          {                                                                   \
            fprintf (stderr, "HDF5 call in file %s line %d failed: %s\n",     \
                     __FILE__, __LINE__, #fn_call);                           \
          }                                                                   \
        }

/* the shell environment variable to query
   for an optional HDF5 input filename */
#define DX_IMPORT_HDF5  "DX_IMPORT_HDF5"

/*****************************************************************************
 **************************     Type Definitions   ***************************
 *****************************************************************************/

/* description of a dataset */
typedef struct
{
  char *name;             /* name of the dataset */
  double time;            /* value of the "time" attribute (if any) */
  int ndims;              /* number of dimensions */
  hsize_t *dims;          /* [ndims] */
  float *origin, *delta;  /* [ndims] */
  size_t typesize;
  H5T_class_t typeclass;

  /* hyperslab parameters to use */
  int single_precision;     /* whether to convert into single precision */
  hssize_t *start;          /* [ndims] */
  hsize_t *stride, *count;  /* [ndims] */
} dataset_t;

/* structure to keep state information as persistent data
   for each ImportCactusHDF5 module's instance */
typedef struct
{
  hid_t fid;              /* file descriptor for the input file */
  int is_gridftp_file;    /* flag indicating if input file is remote */
  char *last_filename;    /* name of the last opened input file */
  int last_index;         /* largest possible dataset index */
  int allocated;          /* number of allocated dataset structures */
  dataset_t *datasets;    /* list of datasets */
  char *objectname;       /* name of file object to be browsed */
  const char *file_env;   /* name of file passed in environment */
  int vector_import;      /* wheher to import as a vector array */
} file_t;


/*****************************************************************************
 *************************     Function Prototypes   *************************
 *****************************************************************************/

/* prototype of the exported worker routine of module ImportHDF5 */
Error m_ImportHDF5 (Object *in, Object *out);

/* prototypes of routines defined locally in this source file */
static int OpenFile (const char *filename, file_t *file, const Object *in);
static herr_t ExamineFile (hid_t group, const char *objectname, void *arg);
static Field ReadDataset (const file_t *file, dataset_t *this);
static void AddAttributes (const file_t *file, const dataset_t *this,
                           Object object);
static int CompareDatasetNames (const void *_a, const void *_b);


 /*@@
   @routine    m_ImportHDF5
   @date       Tue 23 Jan 2001
   @author     Gerd Lanfermann, Thomas Radke
   @desc
               The main worker routine of the ImportHDF5 module.

               It imports data from an HDF5 file <filename>.
               On the first time through, or if the <reopen> flag is set,
               the file will be browsed to get a list of all datasets contained
               in it. Individual datasets from that list can then be addressed
               by their index. If the datasets have a "time" attribute attached
               to them the list will be sorted by their values.
               The value of the largest possible dataset index is output on
               outtab[1].

               A specific dataset, selected by <index>, is then read from the
               input file - either completely as a full dataset, or - if any
               of the optional parameters <origin>, <thickness>, and <stride>
               were specified - as a hyperslab.
               According to the <single_precision> flag, double precision
               floating-point data will be read in single (default) or double
               precision.
               Each dataset is assumed to describe a regular grid. "origin"
               and "delta" attributes are used to create a positions and
               connections component for the imported data array.
               The imported field is output on outtab[0].
   @enddesc

   @var        in
   @vdesc      pointer to array of input ports;
               in[0] - string with the input filename
               in[1] - origin of slab to read
               in[2] - thickness of slab to read
               in[3] - stride for slab to read
               in[4] - index or name of the dataset to read
               in[5] - flag whether to reopen the file on each activation
               in[6] - flag whether to read REAL data in single precisison
               in[7] - user name for remote file access using standard Ftp
               in[8] - password for remote file access using standard Ftp
               in[9] - subject name for remote file access using GSI
               in[10] - number of parallel streams to use for remote file access
               in[11] - flag whether to import the dataset as a vector array
   @vtype      Object *
   @vio        in
   @endvar
   @var        out
   @vdesc      pointer to array of output ports;
               out[0] - the data as requested via in[1..4], returned as a Field
               out[1] - the largest possible dataset index for this input file
   @vtype      Object *
   @vio        inout
   @endvar
 @@*/
Error m_ImportHDF5 (Object *in, Object *out)
{
  Type type;
  int i, ndims, newndims, rank, dataset, reopen, retval;
  const int *vector;
  char *myID, *function;
  char *filename, *datasetname;
  Private private;
  file_t *file;
  dataset_t *this;


  /******************************************************************
   * get this module's instance's private data class from the cache *
   ******************************************************************/
  /* NOTE: DXGetModuleId() is assumend to return a string here.
           If this is changed in the OpenDX API at some time, it cannot be
           directly used as a key into the cache any longer. */
  myID = (char *) DXGetModuleId ();
  private = (Pointer) DXGetCacheEntryV (myID, 0, 0, NULL);
  if (private == NULL)
  {
    /* if it doesn't exist yet, create it and put into cache */
    file = calloc (1, sizeof (*file));
    file->fid = -1;
    private = DXNewPrivate (file, NULL);
    DXSetCacheEntryV ((Object) private, CACHE_PERMANENT, myID, 0, 0, NULL);
    file->file_env = getenv (DX_IMPORT_HDF5);
  }
  file = DXGetPrivateData (private);
  DXFreeModuleId (myID);

  /*****************************************************************
   * open the input file on the first time through, or on <reopen> *
   *****************************************************************/
  /* get the name of the HDF5 file to open */

  if (! in[0] && ! file->file_env)
  {
    DXSetError (ERROR_BAD_PARAMETER, "No filename given");
    return (ERROR);
  }

  /* FIXME:
   * Quick hack which means it will always take from
   * environment if present.
   * I.e. in[0] will be ignored.
   */
  if(in[0] && ! file->file_env)
  {
    DXExtractString (in[0], &filename);
  }
  else
  {
    filename = strdup(file->file_env);
  }

  /* get the reopen flag */
  reopen = 0;
  if (in[5])
  {
    DXExtractInteger (in[5], &reopen);
  }

  /* open the file on the first trip through or if the filename has changed */
  if (file->last_filename == NULL ||
      strcmp (filename, file->last_filename) || reopen)
  {
    if (OpenFile (filename, file, &in[7]) == ERROR)
    {
      return (ERROR);
    }

    /* set output port out[2] to max_index */
    out[1] = (Object) DXMakeInteger (file->last_index);
  }

  if (file->last_index < 0)
  {
    DXSetError (ERROR_BAD_PARAMETER,
                "There are no datasets in file '%s'", filename);
    return (ERROR);
  }

  /*****************************************************************
   * read the other input parameters and check for sensible values *
   *****************************************************************/
  /* read the dataset index/name and check that it is valid */
  dataset = 0;
  if (in[4])
  {
    /* is it an integer? */
    if (DXExtractInteger (in[4], &dataset))
    {
      if (dataset < 0 || dataset > file->last_index)
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "Invalid dataset index %d (must be in range [0, %d)",
                    dataset, file->last_index);
      }
    }
    /* if not then it must be a string */
    else if (DXExtractString (in[4], &datasetname))
    {
      /* find the dataset index whose name is datasetname */
      for (dataset = 0; dataset <= file->last_index; dataset++)
      {
        if (strcmp (datasetname, file->datasets[dataset].name) == 0)
        {
          break;
        }
      }
      if (dataset > file->last_index)
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "There is no dataset '%s' in file '%s'",
                    datasetname, filename);
        return (ERROR);
      }
    }
    else
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "Couldn't read dataset name from input tab");
      return (ERROR);
    }

  }
  this = &file->datasets[dataset];

  /* read the vectorimport flag */
  file->vector_import = 0;
  newndims = this->ndims;
  if (in[11])
  {
    DXExtractInteger (in[11], &file->vector_import);
    if (file->vector_import)
    {
      /* reduce the dimensionality by 1 since the last dimension is
         implicitly imported as a vector */
      newndims = this->ndims - 1;
    }
  }

  /* read the hyperslab origin */
  memset (this->start, 0, this->ndims * sizeof (*this->start));
  if (in[1])
  {
    DXGetArrayInfo ((Array) in[1], &ndims, &type, NULL, &rank, NULL);
    if (type != TYPE_INT || (rank != 0 && rank != 1))
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "origin must be given as integer list or integer vector");
      return (ERROR);
    }
    if (rank == 1)
    {
      DXGetArrayInfo ((Array) in[1], NULL, NULL, NULL, NULL, &ndims);
    }
    if (ndims > newndims)
    {
      DXWarning ("origin has %d elements but selected dataset %d has only "
                 "%d dimensions (exceeding elements will be ignored)",
                 ndims, dataset, newndims);
      ndims = newndims;
    }
    vector = DXGetArrayData ((Array) in[1]);
    for (i = 0; i < ndims; i++)
    {
      if (vector[i] < 0 || vector[i] >= (int) this->dims[newndims-i-1])
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "origin[%d] = %d is out of range (must be in range [0, ",
                    "%d])", i, vector[i], (int) this->dims[newndims-i-1]-1);
        return (ERROR);
      }
      this->start[newndims-i-1] = vector[i];
    }
  }

  /* read the hyperslab extent (thickness) */
  for (i = 0; i < this->ndims; i++)
  {
    this->count[i] = this->dims[i] - this->start[i];
  }
  if (in[2])
  {
    DXGetArrayInfo ((Array) in[2], &ndims, &type, NULL, &rank, NULL);
    if (type != TYPE_INT || (rank != 0 && rank != 1))
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "thickness must be given as integer list or integer vector");
      return (ERROR);
    }
    if (rank == 1)
    {
      DXGetArrayInfo ((Array) in[2], NULL, NULL, NULL, NULL, &ndims);
    }
    if (ndims > newndims)
    {
      DXWarning ("thickness vector has %d elements but selected dataset %d has "
                 "only %d dimensions (exceeding elements will be ignored)",
                 ndims, dataset, newndims);
      ndims = newndims;
    }
    vector = DXGetArrayData ((Array) in[2]);
    for (i = 0; i < ndims; i++)
    {
      if (vector[i] < 0 || vector[i] > (int) this->dims[newndims-i-1])
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "thickness[%d] = %d is out of range (must be in range [0, "
                    "%d])", i, vector[i], (int) this->dims[newndims-i-1]);
        return (ERROR);
      }
      if (vector[i])
      {
        this->count[newndims-i-1] = vector[i];
      }
    }
  }

  /* verify that start + extent are still within the dataset dimensions */
  for (i = 0; i < this->ndims; i++)
  {
    if (this->start[i] + this->count[i] > this->dims[i])
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "origin (%d) + thickness (%d) in dimension %d is out of "
                  "range (must be in range [1, %d])", (int) this->start[i],
                  (int) this->count[i], this->ndims-i-1, (int) this->dims[i]);
      return (ERROR);
    }
  }

  /* read the hyperslab stride */
  for (i = 0; i < this->ndims; i++)
  {
    this->stride[i] = 1;
  }
  if (in[3])
  {
    DXGetArrayInfo ((Array) in[3], &ndims, &type, NULL, &rank, NULL);
    if (type != TYPE_INT || (rank != 0 && rank != 1))
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "stride must be given as integer list or integer vector");
      return (ERROR);
    }
    if (rank == 1)
    {
      DXGetArrayInfo ((Array) in[3], NULL, NULL, NULL, NULL, &ndims);
    }
    if (ndims > newndims)
    {
      DXWarning ("stride vector has %d elements but selected dataset %d has "
                 "only %d dimensions (exceeding elements will be ignored)",
                 ndims, dataset, newndims);
      ndims = newndims;
    }
    vector = DXGetArrayData ((Array) in[3]);
    for (i = 0; i < ndims; i++)
    {
      if (vector[i] < 1)
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "stride[%d] = %d is out of range (must be positive)",
                    i, vector[i]);
        return (ERROR);
      }
      this->stride[newndims-i-1] = vector[i];
    }
  }

  /* calculate the real counts according to the strides */
  for (i = 0; i < this->ndims; i++)
  {
    this->count[i] = (this->count[i] + this->stride[i] - 1) / this->stride[i];
  }

  /* read the precision flag */
  this->single_precision = 1;
  if (in[6])
  {
    DXExtractInteger (in[6], &this->single_precision);
  }

  /*****************************************************************
   * now read the requested dataset, either from file or cache     *
   *****************************************************************/
  /* build a unique name for this request (also across multiple files) */
  function = malloc (strlen (filename) + strlen (this->name) + 1);
  sprintf (function, "%s%s", filename, this->name);

  /* check whether the dataset is available from the cache
     if not, read it from the file and put it in the cache */
  out[0] = DXGetCacheEntryV (function, 0, 4, &in[1]);
#ifdef DEBUG
  DXMessage ("requested dataset '%s' (index %d out of %d from file '%s')",
             this->name, dataset, file->last_index, filename);
  DXMessage ("reading dataset from %s...", out[0] ? "cache" : "file");
#endif
  if (! out[0])
  {
    out[0] = (Object) ReadDataset (file, this);
    DXSetCacheEntryV (out[0], (double) 0, function, 0, 4, &in[1]);
  }
  free (function);

  retval = out[0] ? OK : ERROR;
  if (retval == OK)
  {
    if (reopen)
    {
      DXReadyToRun (DXGetModuleId ());
    }
  }

  return (retval);
}


 /*@@
   @routine    OpenFile
   @date       Tue 3 December 2002
   @author     Thomas Radke
   @desc
               Opens a given HDF5 file and browses through all the datasets
               therein.
   @enddesc

   @var        filename
   @vdesc      filename of the file to open
   @vtype      const char *
   @vio        in
   @endvar
   @var        file
   @vdesc      structure to fill out with file and dataset information
   @vtype      file_t *
   @vio        inout
   @endvar
   @var        in
   @vdesc      input ports containing
                - the user name (in[0]) and password (in[1]) for standard Ftp
                  authentification
                - the subject name (in[2]) for GSI authentification
                - the number of parallel streams to use for remote file access
                  (in[3])
   @vtype      const Object *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 for success, or negative if file could not be opened
   @endreturndesc
@@*/
static int OpenFile (const char *filename, file_t *file, const Object *in)
{
  int i;
  hid_t plist;
#ifdef H5_HAVE_GRIDFTP
  H5FD_gridftp_fapl_t gridftp_fapl;
#endif


#ifndef H5_HAVE_GRIDFTP
  /* prevent compiler warning about unused parameter */
  in = in;
#endif

  /*****************************************************************
   * close a previously opened file and free its resources         *
   *****************************************************************/
  if (file->last_filename)
  {
#ifdef DEBUG
    DXMessage ("closing file '%s'", file->last_filename);
#endif
    free (file->last_filename);
    file->last_filename = NULL;
  }
  if (file->fid >= 0)
  {
    CHECK_ERROR (H5Fclose (file->fid));
  }
  if (file->datasets)
  {
    for (i = 0; i <= file->last_index; i++)
    {
      free (file->datasets[i].dims);
      free (file->datasets[i].start);
      free (file->datasets[i].count);
      free (file->datasets[i].origin);
      free (file->datasets[i].name);
    }
    free (file->datasets);
    file->datasets = NULL;
  }
  if (file->objectname)
  {
    free (file->objectname);
    file->objectname = NULL;
  }

  /*****************************************************************
   * open the new file                                             *
   *****************************************************************/
  /*
   * open the file as one of the following until one succeeds:
   *  - a remote HDF5 file using the GridFtp driver
   *  - a streamed HDF5 file using the Stream driver
   *  - a regular HDF5 file on a UNIX filesystem using the standard driver
   */
#ifdef DEBUG
  DXMessage ("opening file '%s'", filename);
#endif
  CHECK_ERROR (plist = H5Pcreate (H5P_FILE_ACCESS));
  file->fid = -1;
  H5E_BEGIN_TRY
  {
#ifdef H5_HAVE_GRIDFTP
    if (strncmp (filename, "ftp://", 6) == 0 ||
        strncmp (filename, "gsiftp://", 9) == 0)
    {
      /* set the property list to use the GridFtp VFD with default properties */
      H5Pset_fapl_gridftp (plist, NULL);

      /* get the GridFtp VFD default properties and modify them
         according to the options from the input tabs */
      H5Pget_fapl_gridftp (plist, &gridftp_fapl);
      if (in[0])
      {
        DXExtractString (in[0], &gridftp_fapl.user);
      }
      if (in[1])
      {
        DXExtractString (in[1], &gridftp_fapl.password);
      }
      if (in[2])
      {
        DXExtractString (in[2], &gridftp_fapl.subject);
      }
      if (in[3])
      {
        DXExtractInteger (in[3], &i);
        gridftp_fapl.num_streams = (unsigned int) i;
      }
      H5Pset_fapl_gridftp (plist, &gridftp_fapl);
      file->fid = H5Fopen (filename, H5F_ACC_RDONLY, plist);
      file->is_gridftp_file = file->fid >= 0;
    }
#endif
    if (file->fid < 0)
    {
#ifdef H5_HAVE_STREAM
      H5Pset_fapl_stream (plist, NULL);
      file->fid = H5Fopen (filename, H5F_ACC_RDONLY, plist);
#endif
      if (file->fid < 0)
      {
        H5Pset_fapl_sec2 (plist);
        file->fid = H5Fopen (filename, H5F_ACC_RDONLY, plist);
      }
    }
  } H5E_END_TRY;
  CHECK_ERROR (H5Pclose (plist));

  if (file->fid < 0)
  {
    DXSetError (ERROR_BAD_PARAMETER, "Could not open HDF5 file '%s'", filename);
    return (ERROR);
  }

  /* initialize file info structure */
  file->last_filename = strdup (filename);

  /*****************************************************************
   * browse through all the datasets contained therein             *
   *****************************************************************/
  file->objectname = strdup ("");
  file->allocated = 0; file->last_index = -1;
  CHECK_ERROR (H5Giterate (file->fid, "/", NULL, ExamineFile, file));

#ifdef DEBUG
  DXMessage ("file '%s' contains %d datasets", filename, file->last_index+1);
#endif

  /* now sort the datasets by their time attribute */
  qsort (file->datasets, (size_t) file->last_index + 1, sizeof (dataset_t),
         CompareDatasetNames);

  return (OK);
}


 /*@@
   @routine    ExamineFile
   @date       Tue 23 Jan 2001
   @author     Frank Herrmann, Gerd Lanfermann, Thomas Radke
   @desc
               The worker routine which is called by H5Giterate().
               It reads a dataset from a file into an Array object
               and adds it to a Field object.
   @enddesc

   @var        group
   @vdesc      HDF5 object to start the iteration
   @vtype      hid_t
   @vio        in
   @endvar
   @var        name
   @vdesc      name of the object at the current iteration
   @vtype      const char *
   @vio        in
   @endvar
   @var        _file
   @vdesc      pointer to the info structure describing the file
   @vtype      void *
   @vio        in
   @endvar
@@*/
static herr_t ExamineFile (hid_t group, const char *name, void *_file)
{
  int nmembers = 0, is_good;
  size_t i, ndims, typesize;
  hsize_t *dims;
  hid_t dataset, datatype, dataspace;
  hid_t attr_origin, attr_delta, attr_time;
  H5T_class_t typeclass, member1typeclass = -1, member0typeclass = -1;
  H5G_stat_t object_info;
  size_t new_size;
  char *objectname;
  dataset_t *this;
  file_t *file = _file;


  /* build the full name for the current object to process */
  objectname = file->objectname;
  file->objectname = malloc (strlen (objectname) + strlen (name) + 2);
  sprintf (file->objectname, "%s/%s", objectname, name);

  /* we are interested in datasets only - skip anything else */
  CHECK_ERROR (H5Gget_objinfo (group, file->objectname, 0, &object_info));
  if (object_info.type != H5G_DATASET)
  {
    if (object_info.type == H5G_GROUP)
    {
      /* skip groups named "Cactus Parameters" */
      if (strcmp (name, "Cactus Parameters"))
      {
        /* iterate over all datasets in this group */
        CHECK_ERROR (H5Giterate (group, file->objectname, NULL, ExamineFile,
                                 file));
      }
    }
    free (file->objectname);
    file->objectname = objectname;

    return (0);
  }

  /* check the dataset's datatype - make sure it is either integer or real */
  CHECK_ERROR (dataset   = H5Dopen (group, name, H5P_DEFAULT));
  CHECK_ERROR (datatype  = H5Dget_type (dataset));
  CHECK_ERROR (typeclass = H5Tget_class (datatype));
  typesize = H5Tget_size (datatype);
  if (typeclass == H5T_COMPOUND)
  {
    hid_t membertype;
    nmembers = H5Tget_nmembers (datatype);
    CHECK_ERROR (membertype = H5Tget_member_type (datatype, 0));
    CHECK_ERROR (member0typeclass = H5Tget_class (membertype));
    CHECK_ERROR (H5Tclose (membertype));
    CHECK_ERROR (membertype = H5Tget_member_type (datatype, 1));
    CHECK_ERROR (member1typeclass = H5Tget_class (membertype));
    CHECK_ERROR (H5Tclose (membertype));
  }
  CHECK_ERROR (H5Tclose (datatype));
  is_good = typeclass == H5T_INTEGER ||
            typeclass == H5T_FLOAT   ||
            (typeclass == H5T_COMPOUND && nmembers == 2 &&
             member0typeclass == H5T_FLOAT && member1typeclass == H5T_FLOAT);
  if (! is_good)
  {
    CHECK_ERROR (H5Dclose (dataset));
    free (file->objectname);
    file->objectname = objectname;

    return (0);
  }

  /* increment dataset counter */
  file->last_index++;

  /* check whether the dataset id array needs to be extended */
  if (file->last_index >= file->allocated)
  {
    file->allocated += H5_BLOCK_INCREMENT;
    new_size = file->allocated * sizeof (dataset_t);
    if (file->datasets)
    {
      file->datasets = realloc (file->datasets, new_size);
    }
    else
    {
      file->datasets = malloc (new_size);
    }
  }

  /* read the dimensions */
  CHECK_ERROR (dataspace = H5Dget_space (dataset));
  ndims = H5Sget_simple_extent_ndims (dataspace);
  dims = malloc (ndims * sizeof (hsize_t));
  CHECK_ERROR (H5Sget_simple_extent_dims (dataspace, dims, NULL));
  CHECK_ERROR (H5Sclose (dataspace));

  this = &file->datasets[file->last_index];
  memset (this, 0, sizeof (*this));
  this->ndims = ndims;
  this->dims = dims;
  this->name = strdup (file->objectname);
  this->typesize = typesize;
  this->typeclass = typeclass;
  this->start = malloc (ndims * sizeof (hssize_t));
  this->count = malloc (2 * ndims * sizeof (hsize_t));
  this->stride = this->count + ndims;

  /* initialize the dataset's attributes to some sensible defaults
     and then try to read them */
  this->origin = malloc (2 * ndims * sizeof (float));
  this->delta = this->origin + ndims;
  for (i = 0; i < ndims; i++)
  {
    this->origin[i] = 0;
    this->delta[i] = 1;
  }
  this->time = file->last_index;
  H5E_BEGIN_TRY
  {
    attr_origin = H5Aopen_name (dataset, "origin");
    attr_delta  = H5Aopen_name (dataset, "delta");
    attr_time   = H5Aopen_name (dataset, "time");
  } H5E_END_TRY;
  if (attr_origin >= 0)
  {
    CHECK_ERROR (H5Aread (attr_origin, H5T_NATIVE_FLOAT, this->origin));
    CHECK_ERROR (H5Aclose (attr_origin));
  }
  if (attr_delta >= 0)
  {
    CHECK_ERROR (H5Aread (attr_delta, H5T_NATIVE_FLOAT, this->delta));
    CHECK_ERROR (H5Aclose (attr_delta));
  }
  if (attr_time >= 0)
  {
    CHECK_ERROR (H5Aread (attr_time, H5T_NATIVE_DOUBLE, &this->time));
    CHECK_ERROR (H5Aclose (attr_time));
  }

  CHECK_ERROR (H5Dclose (dataset));
  free (file->objectname);
  file->objectname = objectname;

  return (0);
}


 /*@@
   @routine    ReadDataset
   @date       Tue 11 June 2002
   @author     Frank Herrmann, Gerd Lanfermann, Thomas Radke
   @desc
               Reads a single dataset from the given HDF5 file.
   @enddesc

   @var        file
   @vdesc      info structure describing the HDF5 input file
   @vtype      const file_t *
   @vio        in
   @endvar
   @var        this
   @vdesc      structure to identify the dataset to read
   @vtype      dataset_t *
   @vio        inout
   @endvar

   @returntype int
   @returndesc
               a new Field object containing the requested dataset,
               or NULL if the Field object coulnd't be created
   @endreturndesc
@@*/
static Field ReadDataset (const file_t *file, dataset_t *this)
{
  Type dxtype;
  Category dxcategory;
  Array array;
  Field retval;
  hid_t dataset, datatype, slabspace, filespace, plist;
  int i, j, newndims, vector_length;
  int *dims;
  float *origin, *delta, *bbox;
  Pointer data;


  /* determine the datatype to use */
  if (this->typeclass == H5T_FLOAT)
  {
    if (this->typesize == sizeof (float) || this->single_precision)
    {
      datatype = H5T_NATIVE_FLOAT; dxtype = TYPE_FLOAT;
    }
    else
    {
      datatype = H5T_NATIVE_DOUBLE; dxtype = TYPE_DOUBLE;
    }
    dxcategory = CATEGORY_REAL;
  }
  else if (this->typeclass == H5T_INTEGER)
  {
    datatype = H5T_NATIVE_INT; dxtype = TYPE_INT; dxcategory = CATEGORY_REAL;
  }
  else /* if (this->typeclass == H5T_COMPOUND) */
  {
    if (this->typesize == 2*sizeof (float) || this->single_precision)
    {
      CHECK_ERROR (datatype = H5Tcreate (H5T_COMPOUND, 2*sizeof (float)));
      CHECK_ERROR (H5Tinsert (datatype, "real", 0, H5T_NATIVE_FLOAT));
      CHECK_ERROR (H5Tinsert (datatype, "imag", sizeof (float),
                              H5T_NATIVE_FLOAT));
      dxtype = TYPE_FLOAT;
    }
    else
    {
      CHECK_ERROR (datatype = H5Tcreate (H5T_COMPOUND, 2*sizeof (double)));
      CHECK_ERROR (H5Tinsert (datatype, "real", 0, H5T_NATIVE_DOUBLE));
      CHECK_ERROR (H5Tinsert (datatype, "imag", sizeof (double),
                              H5T_NATIVE_DOUBLE));
      dxtype = TYPE_DOUBLE;
    }
    dxcategory = CATEGORY_COMPLEX;
  }

  if (file->vector_import)
  {
    /* check the length of the vector */
    newndims = this->ndims - 1;
    vector_length = this->dims[this->ndims - 1];
    array = DXNewArray (dxtype, dxcategory, 1, vector_length);
  }
  else
  {
    newndims = this->ndims;
    vector_length = 1;
    array = DXNewArray (dxtype, dxcategory, 0);
  }
  if (array == NULL)
  {
    DXSetError (ERROR_NO_MEMORY, "Couldn't allocate new DX array");
    return (NULL);
  }

  /* allocate some temporary buffers */
  dims   = calloc ((size_t) this->ndims, sizeof (int));
  origin = calloc ((size_t) (this->ndims+1) * this->ndims, sizeof (float));
  delta  = origin + this->ndims;

  /* by default: select full HDF5 dataspace, copy DX dims, origin, delta */
  for (i = 0; i < this->ndims; i++)
  {
    dims[i] = this->count[i];
    delta[(this->ndims-i-1)*this->ndims + i] =
      this->delta[i] * this->stride[this->ndims-i-1];
    origin[i] = this->origin[i] +
                this->start[this->ndims-i-1] * this->delta[i];
  }

  CHECK_ERROR (slabspace = H5Screate_simple (this->ndims, this->count, NULL));
  i = H5Sget_simple_extent_npoints (slabspace);
  data = DXGetArrayData (DXAddArrayData (array, 0, i/vector_length, NULL));
  if (data)
  {
    CHECK_ERROR (filespace = H5Screate_simple (this->ndims, this->dims, NULL));
    if (i > 1)
    {
      CHECK_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET, this->start,
                                        this->stride, this->count, NULL));
    }

    /* set the buffer size to read the full hyperslab */
    CHECK_ERROR (plist = H5Pcreate (H5P_DATASET_XFER));
    if (! file->is_gridftp_file)
    {
      CHECK_ERROR (H5Pset_buffer (plist, i * this->typesize, NULL, NULL));
    }

    CHECK_ERROR (dataset = H5Dopen (file->fid, this->name, H5P_DEFAULT));
    CHECK_ERROR (H5Dread (dataset, datatype, slabspace, filespace, plist, data));
    CHECK_ERROR (H5Dclose (dataset));
    CHECK_ERROR (H5Pclose (plist));
    CHECK_ERROR (H5Sclose (filespace));

    retval = DXNewField ();
  }
  else
  {
    DXSetError (ERROR_NO_MEMORY, "Couldn't allocate memory for new DX array");
    retval = NULL;
  }

  /* close the dataspace and datatype */
  CHECK_ERROR (H5Sclose (slabspace));
  if (datatype != H5T_NATIVE_INT && datatype != H5T_NATIVE_FLOAT &&
      datatype != H5T_NATIVE_DOUBLE)
  {
    CHECK_ERROR (H5Tclose (datatype));
  }

  if (retval)
  {
    /* add all attributes of this dataset to the new field */
    AddAttributes (file, this, (Object) retval);

    /* set the individual components in the field */
    DXSetComponentValue (retval, "data", (Object) array);
    if (newndims > 0)
    {
      array = DXMakeGridPositionsV (newndims, dims, origin, delta);
      DXSetComponentValue (retval, "positions", (Object) array);

      /* eliminate dimensions along which to slab */
      for (i = j = 0; i < newndims; i++)
      {
        if (this->count[i] > 1)
        {
          dims[j++] = this->count[i];
        }
      }
      array = DXMakeGridConnectionsV (j, dims);
      DXSetComponentValue (retval, "connections", (Object) array);
    }
    DXEndField (retval);

    /* set bounding box according to the full dataset */
    bbox = DXGetArrayData ((Array) DXGetComponentValue (retval, "box"));
    for (i = 0; i < (1 << newndims); i++)
    {
      for (j = 0; j < newndims; j++)
      {
        bbox[i*newndims + j] = this->origin[j];
        if (i & (1 << j))
        {
          bbox[i*newndims + j] += (this->dims[newndims-j-1]-1) *
                                     this->delta[j];
        }
      }
    }
  }
  else
  {
    DXDelete ((Object) array);
  }

  /* free temporary arrays */
  free (origin);
  free (dims);

  return (retval);
}


 /*@@
   @routine    AddAttributes
   @date       Sun 1 June 2003
   @author     Thomas Radke
   @desc
               Reads all attributes from the given dataset and adds them
               to the given DX object.
   @enddesc

   @var        file
   @vdesc      structure with file and dataset information
   @vtype      const file_t *
   @vio        inout
   @endvar
   @var        object
   @vdesc      DX object to attach attributes
   @vtype      Object
   @vio        inout
   @endvar
@@*/
static void AddAttributes (const file_t *file, const dataset_t *this,
                           Object object)
{
  int i, nelems;
  void *value;
  Object dxattr;
  char attrname[256];
  H5T_class_t typeclass;
  hid_t dataset, datatype, dataspace, attr;


  CHECK_ERROR (dataset = H5Dopen (file->fid, this->name, H5P_DEFAULT));
  CHECK_ERROR (i = H5Aget_num_attrs (dataset));
  while (--i >= 0)
  {
    dxattr = NULL;

    CHECK_ERROR (attr = H5Aopen_idx (dataset, (unsigned int) i));
    CHECK_ERROR (H5Aget_name (attr, sizeof (attrname) - 1, attrname));
    CHECK_ERROR (datatype = H5Aget_type (attr));
    CHECK_ERROR (dataspace = H5Aget_space (attr));
    nelems = H5Sget_simple_extent_npoints (dataspace);

    /* determine the datatype to use */
    CHECK_ERROR (typeclass = H5Tget_class (datatype));
    if (typeclass == H5T_INTEGER)
    {
      dxattr = (Object) DXNewArray (TYPE_INT, CATEGORY_REAL, 0);
      value = DXGetArrayData (DXAddArrayData ((Array) dxattr, 0, nelems, NULL));
      CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, value));
    }
    else if (typeclass == H5T_FLOAT)
    {
      dxattr = (Object) DXNewArray (TYPE_DOUBLE, CATEGORY_REAL, 0);
      value = DXGetArrayData (DXAddArrayData ((Array) dxattr, 0, nelems, NULL));
      CHECK_ERROR (H5Aread (attr, H5T_NATIVE_DOUBLE, value));
    }
    else if (typeclass == H5T_STRING)
    {
      nelems = H5Tget_size (datatype);
      value = malloc ((size_t) nelems + 1);
      CHECK_ERROR (H5Aread (attr, datatype, value));
      ((char *) value)[nelems] = 0;
      dxattr = (Object) DXMakeString (value);
      free (value);
    }
    else
    {
      DXWarning ("dataset attribute '%s' is not of type INTEGER or FLOAT "
                 "(attribute will be ignored)", attrname);
    }
    if (dxattr)
    {
      DXSetAttribute (object, attrname, dxattr);
    }

    CHECK_ERROR (H5Sclose (dataspace));
    CHECK_ERROR (H5Tclose (datatype));
    CHECK_ERROR (H5Aclose (attr));
  }
  CHECK_ERROR (H5Dclose (dataset));

  /* also add the name of the dataset itself as an attribute "name"
     if no such attribute already exists */
  if (! DXGetAttribute (object, "name"))
  {
    DXSetStringAttribute (object, "name", this->name);
  }
}


/* comparison function used by qsort(3) to sort datasets by time */
static int CompareDatasetNames (const void *_a, const void *_b)
{
  const dataset_t *a = (const dataset_t *) _a;
  const dataset_t *b = (const dataset_t *) _b;


  return (a->time - b->time > 0 ? +1 : -1);
}
