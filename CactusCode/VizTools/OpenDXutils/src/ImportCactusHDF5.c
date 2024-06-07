 /*@@
   @file      ImportCactusHDF5.c
   @date      Wed 2 April 2003
   @author    Thomas Radke
   @desc
              Data Import module to read datasets from a Cactus HDF5 file
              as unigrid fields into OpenDX.
              In addition to the functionality provided by ImportHDF5,
              this module is also able to read chunked and unchunked Cactus
              datafiles directly.

              See the ../README file for general information.
              See the ../COPYING file for copyright and license information.
              See the ../doc/ subdirectory for documentation.
   @enddesc
   @version   $Header: /cactus/VizTools/OpenDXutils/src/ImportCactusHDF5.c,v 1.11 2006/07/18 09:54:34 tradke Exp $
 @@*/

#include <stdlib.h>
#include <string.h>
#include <unistd.h>   /* sysconf(_SC_OPEN_MAX) */

#include <dx/dx.h>
#include <hdf5.h>


/*****************************************************************************
 *************************     Macro Definitions   ***************************
 *****************************************************************************/

/* uncomment the following line to get some debug output */
/* #define DEBUG 1 */

/* the shell environment variable to query
   for an optional Cactus HDF5 input filename */
#define DX_IMPORT_CACTUS_HDF5  "DX_IMPORT_CACTUS_HDF5"

/* number of file descriptors that should be reserved for system usage
   (these are at least stdin, stdout, stderr; some MPI implementations
    may need additional descriptors for internal use) */
#define RESERVED_FILE_DESCRIPTORS  5

/* number of dataset structures to allocate at once */
#define H5_BLOCK_INCREMENT  100

/* name of the special table-of-contents dataset */
#define TOC_DATASETNAME "/Table of Contents"

/* name of the group which holds global attributes
   for a Cactus HDF5 chunked datafile */
#define GLOBAL_ATTRIBUTES_GROUP "Global Attributes"

/* macro to check the return code of calls to the HDF5 library */
#define CHECK_ERROR(fn_call)                                                  \
        {                                                                     \
          if ((fn_call) < 0)                                                  \
          {                                                                   \
            fprintf (stderr, "HDF5 call in file %s line %d failed: %s\n",     \
                     __FILE__, __LINE__, #fn_call);                           \
          }                                                                   \
        }


/*****************************************************************************
 **************************     Type Definitions   ***************************
 *****************************************************************************/

/* description of a dataset */
typedef struct
{
  char *name;               /* name of the dataset */
  double time;              /* value of the "time" attribute (if any) */
  int ndims;                /* number of dimensions */
  hsize_t *dims;            /* [ndims] */
  float *origin, *delta;    /* [ndims] */
  size_t typesize;
  H5T_class_t typeclass;

  /* hyperslab parameters to use */
  int single_precision;     /* whether to convert into single precision */
  hssize_t *start;          /* [ndims] */
  hsize_t *thickness;       /* [ndims] */
  hsize_t *stride, *count;  /* [ndims] */

  /* the following describe a chunked dataset */
  hsize_t *chunkdims;       /* [ndims * nprocs] */
  hssize_t *chunkoffset;    /* [ndims * nprocs] */
} dataset_t;

/* structure to keep state information as persistent data
   for each ImportCactusHDF5 module's instance */
typedef struct
{
  int last_index;           /* largest possible dataset index */
  int allocated;            /* number of allocated dataset structures */
  dataset_t *datasets;      /* list of datasets */
  char *objectname;         /* name of file object to be browsed */

  /* the following keep information for chunked datafiles */
  int unchunked, nprocs, ioproc_every, nfiles;
  hid_t *fid;
  char **filename;
  int max_filedescriptors;
  int is_gridftp_file;      /* flag indicating if this is a remote HDF5 file */
  char *last_filename;
} file_t;


/*****************************************************************************
 *************************     Function Prototypes   *************************
 *****************************************************************************/

/* prototype of the exported worker routine of module ImportCactusHDF5 */
Error m_ImportCactusHDF5 (Object *in, Object *out);

/* prototypes of routines defined locally in this source file */
static int OpenFile (const char *filename, file_t *file, const Object *in);
static int ReadTOC (file_t *file);
static herr_t ExamineFile (hid_t group, const char *objectname, void *_file);
static int GetChunkLayout (file_t *file, dataset_t *this);
static Field ReadDataset (const file_t *file, dataset_t *this);
static void ReadThisDataset (const file_t *file, const dataset_t *this,
                             hid_t datatype, hid_t slabspace, void *data);
static void AddAttributes (const file_t *file, const dataset_t *this,
                           Object field);
static int CompareDatasetNames (const void *_a, const void *_b);


 /*@@
   @routine    m_ImportCactusHDF5
   @date       Tue 23 Jan 2001
   @author     Gerd Lanfermann, Thomas Radke
   @desc
               The main worker routine of the ImportCactusHDF5 module.

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
               in[4] - index of the dataset to read
               in[5] - flag whether to reopen the file on each activation
               in[6] - flag whether to read REAL data in single precisison
               in[7] - user name for remote file access using standard Ftp
               in[8] - password for remote file access using standard Ftp
               in[9] - subject name for remote file access using GSI
               in[10] - number of parallel streams to use for remote file access
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
Error m_ImportCactusHDF5 (Object *in, Object *out)
{
  Type type;
  int i, ndims, rank, dataset, reopen, retval;
  const int *vector;
  char *myID, *function;
  char *filename;
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
    private = DXNewPrivate (file, NULL);
    DXSetCacheEntryV ((Object) private, CACHE_PERMANENT, myID, 0, 0, NULL);
  }
  file = DXGetPrivateData (private);
  DXFreeModuleId (myID);

  /*****************************************************************
   * open the input file on the first time through, or on <reopen> *
   *****************************************************************/
  /* get the filename to open, either from the environment
     (if DX_IMPORT_CACTUS_HDF5 is set) or from in[0] */
  filename = getenv (DX_IMPORT_CACTUS_HDF5);
  if (filename == NULL)
  {
    if (! in[0])
    {
      DXSetError (ERROR_BAD_PARAMETER, "No filename given");
      return (ERROR);
    }
    DXExtractString (in[0], &filename);
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

  /*****************************************************************
   * read the other input parameters and check for sensible values *
   *****************************************************************/
  /* read the dataset index and check that it is valid */
  dataset = 0;
  if (in[4])
  {
    DXExtractInteger (in[4], &dataset);
  }
  this = dataset >= 0 && dataset <= file->last_index ?
         &file->datasets[dataset] : NULL;
  if (! this)
  {
    if (file->last_index >= 0)
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "Invalid dataset index %d (must be in range [0, %d]",
                  dataset, file->last_index);
    }
    else
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "There are no datasets in file '%s'", filename);
    }
    return (ERROR);
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
    if (ndims > this->ndims)
    {
      DXWarning ("origin has %d elements but selected dataset %d has only "
                 "%d dimensions (exceeding elements will be ignored)",
                 ndims, dataset, this->ndims);
      ndims = this->ndims;
    }
    vector = DXGetArrayData ((Array) in[1]);
    for (i = 0; i < ndims; i++)
    {
      if (vector[i] < 0 || vector[i] >= (int) this->dims[this->ndims-i-1])
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "origin[%d] = %d is out of range (must be in range [0, "
                    "%d])", i, vector[i], (int) this->dims[this->ndims-i-1]-1);
        return (ERROR);
      }
      this->start[this->ndims-i-1] = vector[i];
    }
  }

  /* read the hyperslab thickness */
  for (i = 0; i < this->ndims; i++)
  {
    this->thickness[i] = this->dims[i] - this->start[i];
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
    if (ndims > this->ndims)
    {
      DXWarning ("thickness vector has %d elements but selected dataset %d has "
                 "only %d dimensions (exceeding elements will be ignored)",
                 ndims, dataset, this->ndims);
      ndims = this->ndims;
    }
    vector = DXGetArrayData ((Array) in[2]);
    for (i = 0; i < ndims; i++)
    {
      if (vector[i] < 0 || vector[i] > (int) this->dims[this->ndims-i-1])
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "thickness[%d] = %d is out of range (must be in range [0, "
                    "%d])", i, vector[i], (int) this->dims[this->ndims-i-1]);
        return (ERROR);
      }
      if (vector[i])
      {
        this->thickness[this->ndims-i-1] = vector[i];
      }
    }
  }

  /* verify that start + extent are still within the dataset dimensions */
  for (i = 0; i < this->ndims; i++)
  {
    if (this->start[i] + this->thickness[i] > this->dims[i])
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "origin (%d) + thickness (%d) in dimension %d is out of "
                  "range (must be in range [1, %d])", (int) this->start[i],
                  (int) this->thickness[i], this->ndims-i-1, (int) this->dims[i]);
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
    if (ndims > this->ndims)
    {
      DXWarning ("stride vector has %d elements but selected dataset %d has "
                 "only %d dimensions (exceeding elements will be ignored)",
                 ndims, dataset, this->ndims);
      ndims = this->ndims;
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
      this->stride[this->ndims-i-1] = vector[i];
    }
  }

  /* calculate the real counts from the thickness and strides */
  for (i = 0; i < this->ndims; i++)
  {
    this->count[i] = (this->thickness[i] + this->stride[i]-1) / this->stride[i];
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
  DXMessage ("requested dataset '%s' (index %d out of %d from file '%s'",
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
  char *tmp, *template;
  hid_t file_id, attr, plist, group;
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
  for (i = 0; i < file->nfiles; i++)
  {
    if (file->filename[i])
    {
#ifdef DEBUG
      DXMessage ("closing file '%s'", file->filename[i]);
#endif
      free (file->filename[i]);
    }
    if (file->fid[i] >= 0)
    {
      CHECK_ERROR (H5Fclose (file->fid[i]));
    }
  }
  free (file->filename);
  free (file->fid);
  if (file->datasets)
  {
    for (i = 0; i <= file->last_index; i++)
    {
      free (file->datasets[i].dims);
      free (file->datasets[i].start);
      free (file->datasets[i].count);
      free (file->datasets[i].origin);
      free (file->datasets[i].name);
      if (! file->unchunked)
      {
        free (file->datasets[i].chunkdims);
        free (file->datasets[i].chunkoffset);
      }
    }
    free (file->datasets);
    file->datasets = NULL;
  }
  if (file->objectname)
  {
    free (file->objectname);
    file->objectname = NULL;
  }
  if (file->last_filename)
  {
    free (file->last_filename);
    file->last_filename = NULL;
  }

  /* mark the file structure as being empty */
  file->nfiles = 0;

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
  file_id = -1;
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
      file_id = H5Fopen (filename, H5F_ACC_RDONLY, plist);
      file->is_gridftp_file = file_id >= 0;
    }
#endif
    if (file_id < 0)
    {
#ifdef H5_HAVE_STREAM
      H5Pset_fapl_stream (plist, NULL);
      file_id = H5Fopen (filename, H5F_ACC_RDONLY, plist);
#endif
      if (file_id < 0)
      {
        H5Pset_fapl_sec2 (plist);
        file_id = H5Fopen (filename, H5F_ACC_RDONLY, plist);
      }
    }
  } H5E_END_TRY;
  CHECK_ERROR (H5Pclose (plist));

  if (file_id < 0)
  {
    DXSetError (ERROR_BAD_PARAMETER, "Could not open HDF5 file '%s'", filename);
    return (ERROR);
  }

  /*****************************************************************
   * check whether we deal with chunked or unchuned datafiles      *
   *****************************************************************/
  /* read the 'nprocs' and 'ioproc_every' attributes
     from the GLOBAL_ATTRIBUTES_GROUP group */
  file->nprocs = file->ioproc_every = file->unchunked = 0;
  H5E_BEGIN_TRY
  {
    group = H5Gopen (file_id, GLOBAL_ATTRIBUTES_GROUP, H5P_DEFAULT);
  } H5E_END_TRY;
  if (group >= 0)
  {
    CHECK_ERROR (attr = H5Aopen_name (group, "nprocs"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &file->nprocs));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (group, "ioproc_every"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &file->ioproc_every));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (attr = H5Aopen_name (group, "unchunked"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &file->unchunked));
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (H5Gclose (group));
  }
  else
  {
    CHECK_ERROR (H5Fclose (file_id));
    DXSetError (ERROR_BAD_PARAMETER, "Could not find chunk attributes in HDF5 "
                "file '%s' ! Is this really a Cactus HDF5 file ?", filename);
    return (ERROR);
  }

  /* get the number of chunked input files to process */
  file->nfiles = (file->nprocs + file->ioproc_every - 1) / file->ioproc_every;

  /* query maximum number of files which can be opened at the same time */
  file->max_filedescriptors = sysconf (_SC_OPEN_MAX);
  if (file->max_filedescriptors < 0)
  {
    DXWarning ("Cannot query file descriptor limit ! Assuming none...");
  }
  /* subtract by number of reserved file descriptors */
  file->max_filedescriptors -= RESERVED_FILE_DESCRIPTORS;
  if (file->max_filedescriptors < 0)
  {
    file->max_filedescriptors = file->nfiles;
  }

  /* initialize file info structure */
  file->fid = malloc (file->nfiles * sizeof (hid_t));
  file->filename = malloc (file->nfiles * sizeof (char *));
  file->last_filename = strdup (filename);

  /* now get all chunked input filenames and check that they can be opened */
  if (file->nfiles == 1)
  {
    /* not much to be done here */
    file->fid[0] = file_id;
    file->filename[0] = strdup (filename);
  }
  else
  {
    /* close the first file (it might not be the one written by processor 0) */
    CHECK_ERROR (H5Fclose (file_id));

    /* get the basename of input file(s) */
    tmp = strstr (filename, ".file_");
    if (tmp == NULL)
    {
      DXSetError (ERROR_BAD_PARAMETER, "Cannot parse HDF5 filename '%s' ! "
                  "Is this really a chunked Cactus HDF5 file ?", filename);
      return (ERROR);
    }

    /* build the input filename template */
    template = malloc (strlen (filename) + 2);
    strcpy (template, filename);
    template[tmp - filename + 6] = 0;
    strcat (template, "%d.h5");

    /* now loop through all the files */
    for (i = 0; i < file->nfiles; i++)
    {
      /* build the input filename */
      file->filename[i] = malloc (strlen (template) + 10);
      sprintf (file->filename[i], template, i);
      file->fid[i] = H5Fopen (file->filename[i], H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file->fid[i] < 0)
      {
        DXSetError (ERROR_BAD_PARAMETER, "Cannot open chunked HDF5 file '%s'",
                    file->filename[i]);
        return (ERROR);
      }

      /* close file if file descriptor limit would be exceeded */
      if (i > file->max_filedescriptors)
      {
        CHECK_ERROR (H5Fclose (file->fid[i]));
        file->fid[i] = -1;
      }
    }

    free (template);
  }

  /* try to read metadata from a TOC dataset;
     if that fails then do it the hard way by iterating over all datasets
     starting from "/" in the HDF5 file */
  file->objectname = strdup ("");
  file->allocated = 0; file->last_index = -1;
  if (ReadTOC (file) < 0)
  {
    CHECK_ERROR (H5Giterate (file->fid[0], "/", NULL, ExamineFile, file));
  }

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
  int nmembers = 0;
  size_t i, ndims, is_valid_dataset;
  hsize_t tmp;
  hsize_t *dims;
  hid_t object, datatype, dataspace;
  hid_t member0typeclass = -1, member1typeclass = -1;
  hid_t attr_global_size, attr_origin, attr_delta, attr_time;
  H5G_stat_t object_info;
  size_t new_size;
  char *objectname;
  dataset_t *this;
  file_t *file = _file;


#ifdef DEBUG
  DXMessage ("Examining object '%s'", name);
#endif

  ndims = 0;
  dims = NULL;
  object = -1;

  /* build the full name for the current object to process */
  objectname = file->objectname;
  file->objectname = malloc (strlen (objectname) + strlen (name) + 2);
  sprintf (file->objectname, "%s/%s", objectname, name);

  /* check the type of the current object */
  CHECK_ERROR (H5Gget_objinfo (group, file->objectname, 0, &object_info));
  if (object_info.type == H5G_GROUP)
  {
    /* check whether this is a group of chunked datasets */
    CHECK_ERROR (object = H5Gopen (group, name, H5P_DEFAULT));

    /* a group with chunked data must have a 'global_size' attribute */
    H5E_BEGIN_TRY
    {
      attr_global_size = H5Aopen_name (object, "global_size");
      if (attr_global_size >= 0)
      {
        /* read the global size of the unchunked dataset */
        datatype = H5Aget_type (attr_global_size);
        dataspace = H5Aget_space (attr_global_size);
        ndims = H5Sget_simple_extent_npoints (dataspace);
        if (H5Tget_class (datatype) != H5T_INTEGER)
        {
          ndims = 0;
        }
        if (ndims > 0)
        {
          dims = calloc (ndims, sizeof (hsize_t));
          H5Aread (attr_global_size, H5T_NATIVE_HSIZE, dims);

          /* convert dims[] from fortran order (HDF5) into C (DX) */
          for (i = 0; i < ndims / 2; i++)
          {
            tmp = dims[i];
            dims[i] = dims[ndims - i - 1];
            dims[ndims - i - 1] = tmp;
          }
        }
        H5Sclose (dataspace);
        H5Tclose (datatype);
        H5Aclose (attr_global_size);
      }
    } H5E_END_TRY;

    if (! dims)
    {
      /* if a 'global_size' attribute is missing, we treat it as an ordinary
         group and (recursively) iterate over its datasets
         Note that the special group named "Cactus Parameters" is not
         traversed. */
      CHECK_ERROR (H5Gclose (object));
      if (strcmp (name, "Cactus Parameters"))
      {
        CHECK_ERROR (H5Giterate (group, file->objectname, NULL, ExamineFile,
                                 file));
      }

      free (file->objectname);
      file->objectname = objectname;

      return (0);
    }
  }
  else if (object_info.type != H5G_DATASET)
  {
    /* object is neither a group nor a dataset - silently ignore it */
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

  /* read the dimensions for an unchunked dataset */
  if (object_info.type == H5G_DATASET)
  {
    CHECK_ERROR (object   = H5Dopen (group, name, H5P_DEFAULT));
    CHECK_ERROR (dataspace = H5Dget_space (object));
    ndims = H5Sget_simple_extent_ndims (dataspace);
    dims = malloc (ndims * sizeof (hsize_t));
    CHECK_ERROR (H5Sget_simple_extent_dims (dataspace, dims, NULL));
    CHECK_ERROR (H5Sclose (dataspace));
  }

  this = &file->datasets[file->last_index];
  memset (this, 0, sizeof (*this));
  this->ndims = ndims;
  this->dims = dims;
  this->name = strdup (file->objectname);
  this->start = malloc (ndims * sizeof (hssize_t));
  this->count = malloc (3 * ndims * sizeof (hsize_t));
  this->stride = this->count + ndims;
  this->thickness = this->count + 2*ndims;

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
    attr_origin = H5Aopen_name (object, "origin");
    attr_delta  = H5Aopen_name (object, "delta");
    attr_time   = H5Aopen_name (object, "time");
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

  /* check the object's datatype - make sure it is either integer or real */
  if (object_info.type == H5G_GROUP)
  {
    /* also get the layout of all chunks */
    is_valid_dataset = GetChunkLayout (file, this);
  }
  else
  {
    CHECK_ERROR (datatype = H5Dget_type (object));
    CHECK_ERROR (this->typeclass = H5Tget_class (datatype));
    this->typesize = H5Tget_size (datatype);
    if (this->typeclass == H5T_COMPOUND)
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
    is_valid_dataset = this->typeclass == H5T_INTEGER ||
                       this->typeclass == H5T_FLOAT ||
                       (this->typeclass == H5T_COMPOUND && nmembers == 2 &&
                        member0typeclass == H5T_FLOAT &&
                        member1typeclass == H5T_FLOAT);
  }

  /* close the object */
  if (object_info.type == H5G_GROUP)
  {
    CHECK_ERROR (H5Gclose (object));
  }
  else
  {
    CHECK_ERROR (H5Dclose (object));
  }

  if (! is_valid_dataset)
  {
    free (this->start);
    free (this->count);
    free (this->origin);
    file->last_index--;
  }
#ifdef DEBUG
  else if (file->unchunked)
  {
    DXMessage ("Dataset '%s' is unchunked", file->objectname);
  }
  else
  {
    DXMessage ("Dataset '%s' contains %d chunks in %d files",
               file->objectname, file->nprocs, file->nfiles);
  }
#endif

  free (file->objectname);
  file->objectname = objectname;

  return (0);
}


 /*@@
   @routine    GetChunkLayout
   @date       Sun 1 June 2003
   @author     Thomas Radke
   @desc
               Determines the file layout for a chunked dataset.
   @enddesc

   @var        file
   @vdesc      info structure describing the HDF5 input file
   @vtype      file_t *
   @vio        inout
   @endvar
   @var        this
   @vdesc      structure to identify the dataset to read
   @vtype      dataset_t
   @vio        inout
   @endvar

   @returntype int
   @returndesc
               TRUE(1) for success, or FALSE(0) if chunks have invalid datatype
   @endreturndesc
@@*/
static int GetChunkLayout (file_t *file, dataset_t *this)
{
  int i, j, idx, is_valid_dataset, nmembers = 0, is_good;
  H5T_class_t member0typeclass = -1, member1typeclass = -1;
  unsigned int nchunks, chunk;
  char chunkname[30];
  hid_t group, dataset, attr, dataspace, datatype;
  hssize_t tmp;
  const dataset_t *first;


  first = &file->datasets[0];

  nchunks = this->ndims * file->nprocs;
  this->chunkdims = calloc (nchunks, sizeof (hsize_t));
  this->chunkoffset = calloc (nchunks, sizeof (hssize_t));

  /*
   * In order to save time for browsing through each chunked file for every
   * dataset, we first check whether the current dataset has the same
   * dimensionality as the very first one.
   * If so it will be assumed that the current dataset has the same layout
   * for all its chunks as the first dataset. Browsing can be skipped in
   * this case.
   */
  if (file->last_index > 0 && this->ndims == first->ndims &&
      memcmp (this->dims, first->dims, this->ndims * sizeof (hsize_t)) == 0)
  {
#ifdef DEBUG
    DXMessage ("Dataset '%s' has same topology as first one", this->name);
#endif

    this->typeclass = first->typeclass;
    this->typesize = first->typesize;
    memcpy (this->chunkdims, first->chunkdims, nchunks * sizeof (hsize_t));
    memcpy (this->chunkoffset, first->chunkoffset, nchunks * sizeof (hssize_t));

    return (1);
  }

#ifdef DEBUG
  DXMessage ("Getting chunk info of dataset '%s'", this->name);
#endif

  for (idx = j = 0, is_valid_dataset = 1;
       idx < file->nfiles && is_valid_dataset; idx++)
  {
    /* re-open the file if it was closed before */
    if (idx > file->max_filedescriptors)
    {
#ifdef DEBUG
      DXMessage ("reopening input file '%s'", file->filename[idx]);
#endif
      CHECK_ERROR (file->fid[idx] =
                   H5Fopen (file->filename[idx], H5F_ACC_RDONLY, H5P_DEFAULT));
    }

    /* get the number of chunks in this file */
    nchunks = file->ioproc_every;
    if (idx == file->nfiles - 1 && file->nprocs % file->ioproc_every)
    {
      nchunks = file->nprocs % file->nfiles;
    }

    /* open the chunk group and iterate over all individual chunks */
    CHECK_ERROR (group = H5Gopen (file->fid[idx], file->objectname, H5P_DEFAULT));
    for (chunk = 0; chunk < nchunks; chunk++, j += this->ndims)
    {
      /* check that such a chunk actually exists (it will not if the processor
         didn't contribute to a hyperslab selection during the output) */
      sprintf (chunkname, "chunk%d", chunk);
      H5E_BEGIN_TRY
      {
        dataset = H5Dopen (group, chunkname, H5P_DEFAULT);
      }
      H5E_END_TRY
      if (dataset < 0)
      {
        continue;
      }

      /* check the chunk's datatype - make sure it is either integer or real */
      CHECK_ERROR (datatype = H5Dget_type (dataset));
      CHECK_ERROR (this->typeclass = H5Tget_class (datatype));
      this->typesize = H5Tget_size (datatype);
      if (this->typeclass == H5T_COMPOUND)
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
      is_good = this->typeclass == H5T_INTEGER ||
                this->typeclass == H5T_FLOAT   ||
                (this->typeclass == H5T_COMPOUND && nmembers == 2 &&
                 member0typeclass == H5T_FLOAT && member1typeclass ==H5T_FLOAT);
      if (! is_good)
      {
        CHECK_ERROR (H5Dclose (dataset));
        is_valid_dataset = 0;
        break;
      }

      /* read the 'chunk_origin' attribute of this chunk */
      CHECK_ERROR (attr = H5Aopen_name (dataset, "chunk_origin"));
      CHECK_ERROR (H5Aread (attr, H5T_NATIVE_HSSIZE, this->chunkoffset + j));
      CHECK_ERROR (H5Aclose (attr));

      /* get the chunk dims */
      CHECK_ERROR (dataspace = H5Dget_space (dataset));
      CHECK_ERROR (H5Sget_simple_extent_dims (dataspace, this->chunkdims + j,
                                              NULL));
      CHECK_ERROR (H5Sclose (dataspace));
      CHECK_ERROR (H5Dclose (dataset));

      /* convert chunkoffset[] from fortran order (HDF5) into C (DX) */
      for (i = 0; i < this->ndims / 2; i++)
      {
        tmp = this->chunkoffset[j + i];
        this->chunkoffset[j + i] = this->chunkoffset[j + this->ndims-i-1];
        this->chunkoffset[j + this->ndims-i-1] = tmp;
      }
    }

    /* close the group again */
    CHECK_ERROR (H5Gclose (group));

    /* close input file if file descriptor limit would be exceeded */
    if (idx > file->max_filedescriptors)
    {
#ifdef DEBUG
      DXMessage ("Temporarily closing input file '%s'", file->filename[idx]);
#endif
      CHECK_ERROR (H5Fclose (file->fid[idx]));
    }
  }

  if (! is_valid_dataset)
  {
    free (this->chunkdims);
    free (this->chunkoffset);
  }

  return (is_valid_dataset);
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
  hid_t datatype, slabspace;
  int i, j;
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
    datatype = H5T_NATIVE_INT; dxtype = TYPE_INT;
    dxcategory = CATEGORY_REAL;
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

  array = DXNewArray (dxtype, dxcategory, 0);
  if (array == NULL)
  {
    DXSetError (ERROR_NO_MEMORY, "Couldn't allocate new DX array");
    return (NULL);
  }

  /* allocate some temporary buffers */
  dims   = malloc (this->ndims * sizeof (int));
  origin = calloc ((size_t) (this->ndims+1) * this->ndims, sizeof (float));
  delta = origin + this->ndims;

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
  data = DXGetArrayData (DXAddArrayData (array, 0, i, NULL));
  if (data)
  {
    ReadThisDataset (file, this, datatype, slabspace, data);
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
    array = DXMakeGridPositionsV (this->ndims, dims, origin, delta);
    DXSetComponentValue (retval, "positions", (Object) array);

    /* eliminate dimensions along which to slab */
    for (i = j = 0; i < this->ndims; i++)
    {
      if (this->count[i] > 1)
      {
        dims[j++] = this->count[i];
      }
    }
    array = DXMakeGridConnectionsV (j, dims);
    DXSetComponentValue (retval, "connections", (Object) array);
    DXEndField (retval);

    /* set bounding box according to the full dataset */
    bbox = DXGetArrayData ((Array) DXGetComponentValue (retval, "box"));
    for (i = 0; i < (1 << this->ndims); i++)
    {
      for (j = 0; j < this->ndims; j++)
      {
        bbox[i*this->ndims + j] = this->origin[j];
        if (i & (1 << j))
        {
          bbox[i*this->ndims + j] += (this->dims[this->ndims-j-1]-1) *
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
   @routine    ReadThisDataset
   @date       Sun 1 June 2003
   @author     Thomas Radke
   @desc
               Reads the given dataset from the input file as a hyperslab.
   @enddesc

   @var        file
   @vdesc      structure with file and dataset information
   @vtype      const file_t *
   @vio        inout
   @endvar
   @var        this
   @vdesc      structure to identify the dataset to read
   @vtype      const dataset_t *
   @vio        in
   @endvar
   @var        datatype
   @vdesc      memory data type to use for the hyperslab read
   @vtype      hid_t
   @vio        in
   @endvar
   @var        slabspace
   @vdesc      memory data space to use for the hyperslab read
   @vtype      hid_t
   @vio        in
   @endvar
   @var        data
   @vdesc      pointer to memory buffer for the hyperslab to read
   @vtype      void *
   @vio        out
   @endvar
@@*/
static void ReadThisDataset (const file_t *file, const dataset_t *this,
                             hid_t datatype, hid_t slabspace, void *data)
{
  int i, j, idx;
  unsigned int nchunks, chunk;
  hid_t dataset, group, filespace, plist;
  char chunkname[30];
  hsize_t *chunkcount;
  hssize_t *chunkstart, *chunkoffset;


  if (file->unchunked)
  {
    CHECK_ERROR (filespace = H5Screate_simple (this->ndims, this->dims, NULL));
    CHECK_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET, this->start,
                                      this->stride, this->count, NULL));

    /* set the buffer size to read the full hyperslab */
    CHECK_ERROR (plist = H5Pcreate (H5P_DATASET_XFER));
    if (! file->is_gridftp_file)
    {
      i = H5Sget_simple_extent_npoints (slabspace);
      CHECK_ERROR (H5Pset_buffer (plist, i * this->typesize, NULL, NULL));
    }

    CHECK_ERROR (dataset = H5Dopen (file->fid[0], this->name, H5P_DEFAULT));
    CHECK_ERROR (H5Dread (dataset, datatype, slabspace, filespace, plist,data));
    CHECK_ERROR (H5Dclose (dataset));
    CHECK_ERROR (H5Pclose (plist));
    CHECK_ERROR (H5Sclose (filespace));

    return;
  }

  chunkstart = malloc (2 * this->ndims * sizeof (hssize_t));
  chunkoffset = chunkstart + this->ndims;
  chunkcount = malloc (this->ndims * sizeof (hsize_t));

  /* now read all the chunks from all input files */
  for (idx = j = 0; idx < file->nfiles; idx++)
  {
    /* re-open the file if it was closed before */
    if (idx > file->max_filedescriptors)
    {
#ifdef DEBUG
      DXMessage ("Reopening input file '%s'", file->filename[idx]);
#endif
      CHECK_ERROR (file->fid[idx] =
                   H5Fopen (file->filename[idx], H5F_ACC_RDONLY, H5P_DEFAULT));
    }

    /* get the number of chunks in this file */
    nchunks = file->ioproc_every;
    if (idx == file->nfiles - 1 && file->nprocs % file->ioproc_every)
    {
      nchunks = file->nprocs % file->nfiles;
    }

    for (chunk = 0; chunk < nchunks; chunk++, j += this->ndims)
    {
      memcpy (chunkcount, this->chunkdims+j, this->ndims * sizeof (hsize_t));
      memcpy (chunkoffset, this->chunkoffset+j, this->ndims * sizeof(hssize_t));

      for (i = 0; i < this->ndims; i++)
      {
        if (chunkoffset[i] >= this->start[i] + (hssize_t) this->thickness[i] ||
            chunkoffset[i] + (hssize_t) chunkcount[i] <= this->start[i])
        {
          break;
        }

        /* get the first point on this chunk */
        if (this->start[i] >= chunkoffset[i])
        {
          chunkstart[i] = this->start[i] - chunkoffset[i];
        }
        else
        {
          chunkstart[i] = (chunkoffset[i] - this->start[i] +
                           this->stride[i] - 1) / this->stride[i];
          chunkstart[i] = chunkstart[i]*this->stride[i] -
                          (chunkoffset[i] - this->start[i]);
          /* make sure the point is within this chunk */
          if (chunkstart[i] >= (hssize_t) chunkcount[i] ||
              chunkstart[i] + (hssize_t) chunkoffset[i] >=
              this->start[i] + (hssize_t) this->thickness[i])
          {
            break;
          }
        }

        /* get the number of points from this chunk */
        chunkcount[i] -= chunkstart[i];
        if (chunkcount[i] > this->thickness[i])
        {
          chunkcount[i] = this->thickness[i];
        }
        if (chunkcount[i] + chunkoffset[i] + chunkstart[i] >
            this->start[i] + this->thickness[i])
        {
          chunkcount[i] = this->start[i] + this->thickness[i] -
                          (chunkoffset[i] + chunkstart[i]);
        }
        /* calculate the real chunk counts according to the strides */
        chunkcount[i] = (chunkcount[i] + this->stride[i] - 1) / this->stride[i];

        /* calculate the chunk offsets according to the strides */
        if (this->start[i] < chunkoffset[i])
        {
          chunkoffset[i] = (chunkoffset[i] - this->start[i] +
                            this->stride[i] - 1) / this->stride[i];
        }
        else
        {
          chunkoffset[i] = 0;
        }
      }
      if (i != this->ndims)
      {
        continue;
      }

      /* define the filespace/memspace mapping as hyperslab selections */
      CHECK_ERROR (filespace = H5Screate_simple (this->ndims, this->chunkdims+j,
                                                 NULL));
      CHECK_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET, chunkstart,
                                        this->stride, chunkcount, NULL));

      CHECK_ERROR (H5Sselect_hyperslab (slabspace, H5S_SELECT_SET, chunkoffset,
                                        NULL, chunkcount, NULL));

      /* build the object name of this chunk */
      sprintf (chunkname, "chunk%d", chunk);
      CHECK_ERROR (group = H5Gopen (file->fid[idx], this->name, H5P_DEFAULT));
      CHECK_ERROR (dataset = H5Dopen (group, chunkname, H5P_DEFAULT));
      CHECK_ERROR (H5Dread (dataset, datatype, slabspace, filespace,
                            H5P_DEFAULT, data));
      CHECK_ERROR (H5Dclose (dataset));
      CHECK_ERROR (H5Gclose (group));
      CHECK_ERROR (H5Sclose (filespace));
    }

    /* close input file if file descriptor limit would be exceeded */
    if (idx > file->max_filedescriptors)
    {
#ifdef DEBUG
      DXMessage ("Temporarily closing input file '%s'", file->filename[idx]);
#endif
      CHECK_ERROR (H5Fclose (file->fid[idx]));
    }
  }

  free (chunkstart);
  free (chunkcount);
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
   @var        this
   @vdesc      structure to identify the dataset to read
   @vtype      const dataset_t *
   @vio        in
   @endvar
   @var        field
   @vdesc      DX field to attach attributes
   @vtype      Object
   @vio        inout
   @endvar
@@*/
static void AddAttributes (const file_t *file, const dataset_t *this,
                           Object field)
{
  int i, nelems;
  void *value;
  Object dxattr;
  char attrname[256];
  H5T_class_t typeclass;
  hid_t object, datatype, dataspace, attr;


  if (file->unchunked)
  {
    CHECK_ERROR (object = H5Dopen (file->fid[0], this->name, H5P_DEFAULT));
  }
  else
  {
    CHECK_ERROR (object = H5Gopen (file->fid[0], this->name, H5P_DEFAULT));
  }

  CHECK_ERROR (i = H5Aget_num_attrs (object));
  while (--i >= 0)
  {
    dxattr = NULL;

    CHECK_ERROR (attr = H5Aopen_idx (object, (unsigned int) i));
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
      DXWarning ("object attribute '%s' is not of type INTEGER or FLOAT "
                 "(attribute will be ignored)", attrname);
    }
    if (dxattr)
    {
      DXSetAttribute (field, attrname, dxattr);
    }

    CHECK_ERROR (H5Sclose (dataspace));
    CHECK_ERROR (H5Tclose (datatype));
    CHECK_ERROR (H5Aclose (attr));
  }

  if (file->unchunked)
  {
    CHECK_ERROR (H5Dclose (object));
  }
  else
  {
    CHECK_ERROR (H5Gclose (object));
  }

  /* also add the name of the object itself as an attribute "name"
     if no such attribute already exists */
  if (! DXGetAttribute (field, "name"))
  {
    DXSetStringAttribute (field, "name", this->name);
  }
}


/* comparison function used by qsort(3) to sort datasets by time */
static int CompareDatasetNames (const void *_a, const void *_b)
{
  const dataset_t *a = _a, *b = _b;


  return (a->time - b->time > 0 ? +1 : -1);
}


 /*@@
   @routine    ReadTOC
   @date       Tue 3 December 2002
   @author     Thomas Radke
   @desc
               Tries to read the file and dataset information from a special
               table-of-contents dataset - to shortcut the "official" but
               slow way of using H5Giterate().
   @enddesc

   @var        file
   @vdesc      structure to fill out with file and dataset information
   @vtype      file_t *
   @vio        inout
   @endvar

   @returntype int
   @returndesc
               0 for success, or negative if TOC could not be read
   @endreturndesc
@@*/
static int ReadTOC (file_t *file)
{
  hid_t dataset, dataspace, datatype, plist;
  hid_t string_datatype, hsize_array_datatype, float_array_datatype;
  int i, ndims, retval;
  hsize_t ndatasets, maxdim;
  size_t size, max_namelen;
  size_t offset[6];
  char *buffer;
  dataset_t *this;


  /* initialize HDF5 objects and return value for cleanup */
  dataset = dataspace = datatype = plist = -1;
  string_datatype = hsize_array_datatype = float_array_datatype = -1;
  retval = -1;

  /* try to open a TOC dataset */
  H5E_BEGIN_TRY
  {
    dataset = H5Dopen (file->fid[0], TOC_DATASETNAME, H5P_DEFAULT);
  } H5E_END_TRY;
  if (dataset < 0)
  {
    goto cleanup;
  }

  /* TOC dataspace must be one-dimensional */
  CHECK_ERROR (dataspace = H5Dget_space (dataset));
  ndims = H5Sget_simple_extent_ndims (dataspace);
  if (ndims != 1)
  {
    goto cleanup;
  }
  CHECK_ERROR (H5Sget_simple_extent_dims (dataspace, &ndatasets, NULL));

  /* get the TOC compound datatype */
  CHECK_ERROR (datatype = H5Dget_type (dataset));

  /* check that all structure elements we need are contained
     in this compound datatype */
  if (H5Tget_member_index (datatype, "name") < 0 ||
      H5Tget_member_index (datatype, "time") < 0 ||
      H5Tget_member_index (datatype, "ndims") < 0 ||
      H5Tget_member_index (datatype, "dims") < 0 ||
      H5Tget_member_index (datatype, "origin") < 0 ||
      H5Tget_member_index (datatype, "delta") < 0)
  {
    goto cleanup;
  }

  /* get the string datatype for the dataset names */
  i = H5Tget_member_index (datatype, "name");
  CHECK_ERROR (string_datatype = H5Tget_member_type (datatype, i));
  max_namelen = H5Tget_size (string_datatype);

  /* get maxdim as the dimensionality of the array structure elements */
  i = H5Tget_member_index (datatype, "origin");
  CHECK_ERROR (float_array_datatype = H5Tget_member_type (datatype, i));
  if (H5Tget_array_ndims (float_array_datatype) != 1)
  {
    goto cleanup;
  }
  CHECK_ERROR (H5Tget_array_dims (float_array_datatype, &maxdim));
  CHECK_ERROR (H5Tclose (float_array_datatype));

  CHECK_ERROR (H5Tclose (datatype));

#ifdef DEBUG
  DXMessage ("found TOC dataset with %u entries", (unsigned int) ndatasets);
#endif

  /* create an array datatype of appropriate dimensionality to read the
     origin and delta metadata information */
  CHECK_ERROR (hsize_array_datatype = H5Tarray_create (H5T_NATIVE_HSIZE, 1,
                                                       &maxdim));
  CHECK_ERROR (float_array_datatype = H5Tarray_create (H5T_NATIVE_FLOAT, 1,
                                                       &maxdim));

  /* set the buffer size to read the full TOC dataset at once */
  CHECK_ERROR (plist = H5Pcreate (H5P_DATASET_XFER));
  size = H5Dget_storage_size (dataset);
  CHECK_ERROR (H5Pset_buffer (plist, size, NULL, NULL));

  /* create a corresponding compound datatype (with native datatypes)
     to read the TOC into memory */
  size = max_namelen + sizeof (int) + sizeof (double) +
         maxdim * (2 * sizeof (float) + sizeof (hsize_t));
  CHECK_ERROR (datatype = H5Tcreate (H5T_COMPOUND, size));
  offset[0] = 0;
  CHECK_ERROR (H5Tinsert (datatype, "name", offset[0], string_datatype));
  offset[1] = offset[0] + max_namelen;
  CHECK_ERROR (H5Tinsert (datatype, "time", offset[1], H5T_NATIVE_DOUBLE));
  offset[2] = offset[1] + sizeof (double);
  CHECK_ERROR (H5Tinsert (datatype, "ndims", offset[2], H5T_NATIVE_INT));
  offset[3] = offset[2] + sizeof (int);
  CHECK_ERROR (H5Tinsert (datatype, "dims", offset[3], hsize_array_datatype));
  offset[4] = offset[3] + maxdim * sizeof (hsize_t);
  CHECK_ERROR (H5Tinsert (datatype, "origin", offset[4], float_array_datatype));
  offset[5] = offset[4] + maxdim * sizeof (float);
  CHECK_ERROR (H5Tinsert (datatype, "delta", offset[5], float_array_datatype));

  /* read the TOC dataset all at once */
  buffer = malloc ((size_t) ndatasets * size);
  CHECK_ERROR (H5Dread (dataset, datatype, H5S_ALL, H5S_ALL, plist, buffer));

  /* fill out the dataset[] structure */
  this = file->datasets = malloc ((size_t) ndatasets * sizeof (dataset_t));
  for (i = 0; i < (int) ndatasets; i++, this++)
  {
    this->name = strdup (buffer + offset[0]);
    memcpy (&this->time, buffer + offset[1], sizeof (double));
    memcpy (&this->ndims, buffer + offset[2], sizeof (int));

    this->dims = malloc (this->ndims * sizeof (hsize_t));
    this->origin = malloc (2 * this->ndims * sizeof (float));
    this->delta = this->origin + this->ndims;
    memcpy (this->dims, buffer + offset[3], this->ndims * sizeof (hsize_t));
    memcpy (this->origin, buffer + offset[4], this->ndims * sizeof (float));
    memcpy (this->delta, buffer + offset[5], this->ndims * sizeof (float));

    buffer += size;
  }
  file->last_index = ndatasets - 1;

  /* free read buffer */
  free (buffer -= ndatasets * size);

  /* set the return value to success */
  retval = 0;

cleanup:
  if (plist >= 0) H5Pclose (plist);
  if (dataset >= 0) H5Dclose (dataset);
  if (datatype >= 0) H5Tclose (datatype);
  if (dataspace >= 0) H5Sclose (dataspace);
  if (float_array_datatype >= 0) H5Tclose (float_array_datatype);
  if (hsize_array_datatype >= 0) H5Tclose (hsize_array_datatype);
  if (string_datatype >= 0) H5Tclose (string_datatype);

  return (retval);
}
