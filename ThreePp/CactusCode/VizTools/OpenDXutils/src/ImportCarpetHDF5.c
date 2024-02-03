 /*@@
   @file      ImportCarpetHDF5.c
   @date      Sun 1 June 2003
   @author    Thomas Radke
   @desc
              Data Import module to read datasets from a Carpet HDF5 file
              as a multigrid group of FMR patches into OpenDX

              See the ../README file for general information.
              See the ../COPYING file for copyright and license information.
              See the ../doc/ subdirectory for documentation.
   @enddesc
   @version   $Header: /cactus/VizTools/OpenDXutils/src/ImportCarpetHDF5.c,v 1.23 2007/06/14 16:02:00 tradke Exp $
 @@*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <regex.h>

#include <dx/dx.h>
#include <hdf5.h>


/*****************************************************************************
 *************************     Macro Definitions   ***************************
 *****************************************************************************/

/* uncomment the following line to get some debug output */
/* #define DEBUG 1 */

/* name of the group which holds metadata describing a Carpet datafile
   for a Cactus HDF5 chunked datafile */
#define FILE_METADATA_GROUP "Parameters and Global Attributes"

/* number of patch_t structures to allocate at once */
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


/*****************************************************************************
 **************************     Type Definitions   ***************************
 *****************************************************************************/

/* description of a single FMR patch (a plain HDF5 dataset in the input file) */
typedef struct
{
  char *name;             /* name of the patch */
  hid_t fid;              /* file descriptor of corresponding input file */
  double time;            /* value of the "time" attribute */
  int timestep;           /* value of the "timestep" attribute */
  int level;              /* value of the "level" attribute */
  size_t typesize;
  H5T_class_t typeclass;

  int factor;             /* refinement factor */
  int ndims;              /* number of dimensions */
  hsize_t *dims;          /* [ndims] */
  float *origin, *delta;  /* [ndims] */
  int *iorigin;           /* [ndims] */

  /* hyperslab parameters to use */
  hssize_t *start;          /* [ndims] */
  hsize_t *stride, *count;  /* [ndims] */
  
  /* valid range of points on this patch */
  int *range;    /* [2*ndims] */    

} patch_t;

/* description of a single FMR timestep */
typedef struct
{
  double time;            /* value of the "time" attribute */
  int timestep;           /* value of the "timestep" attribute */
  size_t npatches;        /* number of patches */
  patch_t *patch;         /* [npatches] */
  size_t typesize;
  H5T_class_t typeclass;
  int ndims;              /* number of dimensions */
  int active_levels;

  /* hyperslab parameters to use */
  int single_precision;     /* whether to convert into single precision */
  int *origin;              /* [ndims] */
  int *thickness;           /* [ndims] */
  hsize_t *stride;          /* [ndims] */
  int invalidate_fine_grid; 
} timestep_t;

/* structure to keep state information as persistent data
   for each ImportCactusHDF5 module's instance */
typedef struct
{
  size_t nfiles;          /* number of input files (for parallel output) */
  hid_t *fid;             /* file descriptors for the input file(s) */
  hid_t current_fid;      /* currently open file (when browsing contents) */
  char *last_filename;    /* name of the last opened input file */
  char *objectname;       /* name of file object to be browsed */
  size_t npatches;        /* total number of patches */
  patch_t *patch;         /* list of patches */
  size_t ntimesteps;      /* total number of timesteps */
  timestep_t *timesteps;  /* list of timesteps */
  size_t allocated;       /* number of allocated timestep/patch structures */
  int is_gridftp_file;    /* flag indicating if this is a remote HDF5 file */
  int maxndims;           /* maximum number of dimensions over all timesteps */
  int *maxthickness;      /* [maxndims] */
  int *bboxes;            /* [2 * maxndims * max_level] */
  int max_level;
} file_t;

/* structure containing the data and bbox fields for an individual patch */
typedef struct
{
  Field data;
  Field bbox;
} DataBBox_t;


/*****************************************************************************
 *************************     Function Prototypes   *************************
 *****************************************************************************/

/* prototype of the exported worker routine of module ImportCarpetHDF5 */
Error m_ImportCarpetHDF5 (Object *in, Object *out);

/* prototypes of routines defined locally in this source file */
static int OpenFile (const char *filename,
                     const char *regex, const regex_t *preg,
                     file_t *file, const Object *in);
static herr_t ExamineFile (hid_t group, const char *objectname, void *arg);
static timestep_t *ReadParameters (const Object *in,
                                   const char *regex, const regex_t *preg,
                                   const file_t *file);
static void ReadPatch (timestep_t *timestep, patch_t *this, DataBBox_t *patch);
static void AddAttributes (const patch_t *this, Object object);
static void InvalidateOverlappingBoundaries (const timestep_t *this,
                                             const file_t *file,
                                             Group group);
static int ComparePatches (const void *_a, const void *_b);


 /*@@
   @routine    m_ImportCarpetHDF5
   @date       Sun 1 June 2003
   @author     Thomas Radke
   @desc
               The main worker routine of the ImportCarpetHDF5 module.

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
               in[4] - list of active refinement levels
               in[5] - index of the dataset to read
               in[6] - flag whether to reopen the file on each activation
               in[7] - flag whether to read REAL data in single precisison
               in[8] - user name for remote file access using standard Ftp
               in[9] - password for remote file access using standard Ftp
               in[10] - subject name for remote file access using GSI
               in[11] - number of parallel streams to use for remote file access
               in[12] - string with the names of variables to import
               in[13] - flag whether points are invalidated on fine grid
   @vtype      Object *
   @vio        in
   @endvar
   @var        out
   @vdesc      pointer to array of output ports;
               out[0] - the data as requested via in[1..5], returned as a Group
               out[1] - the bboxes for all patches of active refinement levels
               out[2] - the largest possible dataset index for this input file
   @vtype      Object *
   @vio        inout
   @endvar
 @@*/
Error m_ImportCarpetHDF5 (Object *in, Object *out)
{
  int i, reopen, set_ntimesteps_tab, retval;
  size_t len;
  DataBBox_t patch;
  char patchname[32];
  char *myID, *filename, *function;
  Private private;
  file_t *file;
  char *regex;
  regex_t preg;
  timestep_t *this;


  /* turn off automatic error printing */
  H5Eset_auto (NULL, NULL, NULL);

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
  /* get the name of the HDF5 file to open */
  if (! in[0])
  {
    DXSetError (ERROR_BAD_PARAMETER, "No filename given");
    return (ERROR);
  }
  DXExtractString (in[0], &filename);

  /* get the reopen flag */
  reopen = 0;
  if (in[6])
  {
    DXExtractInteger (in[6], &reopen);
  }

  /* get the varnames parameter */
  regex = ".";
  if (in[12]) {
    DXExtractString (in[12], &regex);
  }
  if (regcomp (&preg, regex, REG_EXTENDED | REG_ICASE)) {
    DXSetError (ERROR_BAD_PARAMETER, "Invalid regular expression '%s' given "
                "for 'varnames' input parameter", regex);
    return (ERROR);
  }

  /* open the file on the first trip through or if the filename has changed */
  set_ntimesteps_tab = 0;
  if (file->last_filename == NULL || strcmp (filename, file->last_filename) ||
      reopen)
  {
    if (OpenFile (filename, regex, &preg, file, &in[8]) == ERROR)
    {
      return (ERROR);
    }

    /* defer setting the 'ntimesteps' output tab until ReadParameters()
       returned successfully */
    set_ntimesteps_tab = 1;
  }

  /*****************************************************************
   * read the other input parameters and check for sensible values *
   *****************************************************************/
  this = ReadParameters (in, regex, &preg, file);
  if (! this)
  {
    return (ERROR);
  }

  if (set_ntimesteps_tab) {
    /* set output port out[2] to max_index */
    out[2] = (Object) DXMakeInteger ((int) file->ntimesteps - 1);
  }

  /*****************************************************************
   * now read all the patches of the requested timestep            *
   *****************************************************************/
  /* build a unique name for this request (also across multiple files) */
  len = strlen (filename) + 1;
  for (i = 0; i < (int) this->npatches; i++) {
    if (regexec (&preg, this->patch[i].name, 0, NULL, 0)) {
      len += strlen (this->patch[i].name);
    }
  }
  function = malloc (len);
  sprintf (function, "%s", filename);
  for (i = 0; i < (int) this->npatches; i++) {
    if (regexec (&preg, this->patch[i].name, 0, NULL, 0)) {
      strcat (function, this->patch[i].name);
    }
  }

  /* check whether the dataset is available from the cache
     if not, read it from the file and put it in the cache */
  out[0] = DXGetCacheEntryV (function, 0, 5, &in[1]);
  out[1] = DXGetCacheEntryV (function, 1, 5, &in[1]);
#ifdef DEBUG
  DXMessage ("requested dataset(s) at time %f (timestep %d) from file '%s'",
             this->time, this->timestep, filename);
  DXMessage ("reading dataset(s) from %s...", out[0] ? "cache" : "file");
#endif
  if (! (out[0] && out[1])) {
    out[0] = (Object) DXNewMultiGrid ();
    out[1] = (Object) DXNewMultiGrid ();
    for (i = 0; i < (int) this->npatches; i++) {

      if (! (this->active_levels & (1 << this->patch[i].level))) continue;

      if (regexec (&preg, this->patch[i].name, 0, NULL, 0)) continue;

      sprintf (patchname, "patch %d", i);
      ReadPatch (this, this->patch + i, &patch);
      if (patch.data)
      {
        DXSetMember ((Group) out[0], patchname, (Object) patch.data);
      }
      if (patch.bbox)
      {
        DXSetMember ((Group) out[1], patchname, (Object) patch.bbox);
      }
    }

    DXGetMemberCount ((Group) out[0], &i);
    if (i <= 0)
    {
      DXWarning ("no patches found for selected refinement level(s) and chosen "
                 "hyperslab parameters");
      DXSetMember ((Group) out[0], "dummy patch", (Object) DXNewField ());
    }
    else
    {
      /* invalidate overlapping boundaries in nested patches */
      InvalidateOverlappingBoundaries (this, file, (Group) out[0]);

#if 0
      /* restore the original bbox (discard the one which depends on invalid
         positions) */
      /*** FIXME: removing a component named "valid box" doesn't seem to work
                  so we have to copy the "box" component onto "valid box" ***/
      DXRename (out[0], "box", "valid box");
#endif
    }
    DXGetMemberCount ((Group) out[1], &i);
    if (i <= 0)
    {
      DXSetMember ((Group) out[1], "dummy patch", (Object) DXNewField ());
    }

    DXSetFloatAttribute ((Object) out[0], "time", this->time);

    DXSetCacheEntryV (out[0], (double) 0, function, 0, 5, &in[1]);
    DXSetCacheEntryV (out[1], (double) 0, function, 1, 5, &in[1]);
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
   @date       Sun 1 June 2003
   @author     Thomas Radke
   @desc
               Opens a given HDF5 file and browses through all the timesteps
               therein.
   @enddesc

   @var        filename
   @vdesc      filename of the file to open
   @vtype      const char *
   @vio        in
   @endvar
   @var        regex
   @vdesc      regular expression string to match against variable names
   @vtype      const char *
   @vio        in
   @endvar
   @var        preg
   @vdesc      compiled regular expression for matching variable names
   @vtype      const regex_t *
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
static int OpenFile (const char *filename,
                     const char *regex, const regex_t *preg,
                     file_t *file, const Object *in)
{
  size_t i, j, k, nmatches;
  int l, maxthickness;
  hid_t plist, group, attr;
  char *msg, *tmp, *template, *chunked_filename;
  timestep_t *timestep;
#ifdef H5_HAVE_GRIDFTP
  H5FD_gridftp_fapl_t gridftp_fapl;
#endif


  /* prevent compiler warnings about unused parameter */
  in = in;

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
  for (i = 0; i < file->nfiles; i++)
  {
    CHECK_ERROR (H5Fclose (file->fid[i]));
  }
  if (file->timesteps)
  {
    for (i = 0; i < file->ntimesteps; i++)
    {
      for (j= 0; j < file->timesteps[i].npatches; j++)
      {
        free (file->timesteps[i].patch[j].start);
        free (file->timesteps[i].patch[j].dims);
        free (file->timesteps[i].patch[j].origin);
        free (file->timesteps[i].patch[j].iorigin);
        free (file->timesteps[i].patch[j].name);
        free (file->timesteps[i].patch[j].range);
      }
      free (file->timesteps[i].stride);
    }
    free (file->bboxes);
    free (file->maxthickness);
    free (file->timesteps);
    free (file->patch);
    file->timesteps = NULL; file->patch = NULL;
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
  file->current_fid = -1;
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
        DXExtractInteger (in[3], &l);
        gridftp_fapl.num_streams = (unsigned int) l;
      }
      H5Pset_fapl_gridftp (plist, &gridftp_fapl);
      file->current_fid = H5Fopen (filename, H5F_ACC_RDONLY, plist);
      file->is_gridftp_file = file->current_fid >= 0;
    }
#endif
    if (file->current_fid < 0)
    {
#ifdef H5_HAVE_STREAM
      H5Pset_fapl_stream (plist, NULL);
      file->current_fid = H5Fopen (filename, H5F_ACC_RDONLY, plist);
#endif
      if (file->current_fid < 0)
      {
        H5Pset_fapl_sec2 (plist);
        file->current_fid = H5Fopen (filename, H5F_ACC_RDONLY, plist);
      }
    }
  } H5E_END_TRY;
  CHECK_ERROR (H5Pclose (plist));

  if (file->current_fid < 0)
  {
    DXSetError (ERROR_BAD_PARAMETER, "Could not open HDF5 file '%s'", filename);
    return (ERROR);
  }

  /*****************************************************************
   * check whether we deal with chunked or unchunked datafiles     *
   *****************************************************************/
  /* read the 'nioprocs' attribute from the file's metadata group
     Old-style Carpet output files may not contain such a group.
     In this case, a single output file is assumed. */
  file->nfiles = 1;
  H5E_BEGIN_TRY
  {
    group = H5Gopen (file->current_fid, FILE_METADATA_GROUP, H5P_DEFAULT);
  } H5E_END_TRY;
  if (group >= 0)
  {
    CHECK_ERROR (attr = H5Aopen_name (group, "nioprocs"));
    CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &file->nfiles));
    assert (file->nfiles > 0);
    CHECK_ERROR (H5Aclose (attr));
    CHECK_ERROR (H5Gclose (group));
  }
  file->fid = malloc (file->nfiles * sizeof (hid_t));

  /* fill in the physical file information */
  if (file->nfiles == 1)
  {
    file->fid[0] = file->current_fid;
  }
  else
  {
    /* close the current file (it might not be first in the set) */
    CHECK_ERROR (H5Fclose (file->current_fid));

    /* get the basename of chunked input file(s) */
    tmp = strstr (filename, ".file_");
    if (tmp == NULL)
    {
      DXSetError (ERROR_BAD_PARAMETER, "Cannot parse HDF5 filename '%s' ! "
                  "Is this really a chunked Carpet HDF5 file ?", filename);
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
      chunked_filename = malloc (strlen (template) + 10);
      sprintf (chunked_filename, template, i);
      file->fid[i] = H5Fopen (chunked_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file->fid[i] < 0)
      {
        DXSetError (ERROR_BAD_PARAMETER, "Cannot open chunked HDF5 file '%s'",
                    chunked_filename);
        free (chunked_filename);
        return (ERROR);
      }
      free (chunked_filename);
    }

    free (template);
  }

  /* initialize file info structure */
  file->last_filename = strdup (filename);

  /*****************************************************************
   * browse through all the datasets contained therein             *
   *****************************************************************/
  file->allocated = file->npatches = file->maxndims = 0;
  for (i = 0; i < file->nfiles; i++)
  {
    file->objectname = strdup ("");
    file->current_fid = file->fid[i];
    CHECK_ERROR (H5Giterate (file->fid[i], "/", NULL, ExamineFile, file));
    free (file->objectname);
  }

  if (file->npatches == 0)
  {
    DXSetError (ERROR_BAD_PARAMETER, "No valid datasets found in HDF5 file "
                "'%s'", filename);
    return (ERROR);
  }

  /* check that this file contains any datasets matching the 'varnames'
     input string parameter */
  nmatches = 0;
  for (i = 0; i < file->npatches; i++) {
    if (regexec (preg, file->patch[i].name, 0, NULL, 0) == 0) {
      nmatches++;
    }
  }
  if (nmatches == 0) {
    DXWarning ("Input file '%s' contains no datasets matching '%s'",
               file->last_filename, regex);
  }

  /* now sort the patches by their time attribute */
  qsort (file->patch, file->npatches, sizeof (patch_t), ComparePatches);

  /* get the maximum refinement level (over all patches in the file) */
  file->max_level = -1;
  for (i = 0; i < file->npatches; i++)
  {
    if (file->max_level < file->patch[i].level)
    {
      file->max_level = file->patch[i].level;
    }
  }

  file->timesteps = malloc (file->npatches * sizeof (timestep_t));
  for (i = 0, timestep = file->timesteps; i < file->npatches; timestep++)
  {
    timestep->patch = file->patch + i;
    timestep->time = timestep->patch[0].time;
    timestep->timestep = timestep->patch[0].timestep;
    timestep->typesize = timestep->patch[0].typesize;
    timestep->typeclass = timestep->patch[0].typeclass;
    timestep->ndims = timestep->patch[0].ndims;

    /* count the number of patches for this timestep */
    timestep->npatches = 0;
    do
    {
      i++;
      timestep->npatches++;
    } while (i < file->npatches &&
             timestep->timestep == file->patch[i].timestep);

    /* compute the refinement factors for each patch */
    /* NOTE: there is a restriction in Carpet that a new level will be refined
       by the same refinement factor in all directions */
    for (j = 0; j < timestep->npatches; j++)
    {
      timestep->patch[j].factor =
        1 << (file->max_level - timestep->patch[j].level);
      for (k = 0; k < (size_t) (timestep->ndims+1) / 2; k++)
      {
        l = timestep->patch[j].iorigin[timestep->ndims - k - 1];
        timestep->patch[j].iorigin[timestep->ndims - k - 1] =
          timestep->patch[j].iorigin[k] * timestep->patch[j].factor;
        timestep->patch[j].iorigin[k] = l * timestep->patch[j].factor;
      }
    }

    /* allocate vectors */
    timestep->stride = malloc (timestep->ndims * sizeof (hsize_t));
    timestep->origin = malloc (2 * timestep->ndims * sizeof (hssize_t));
    timestep->thickness = timestep->origin + 1*timestep->ndims;
  }
  file->ntimesteps = timestep - file->timesteps;

  /* get the dimensions for the underlying grid covering all patches */
  file->maxthickness = calloc ((size_t) file->maxndims, sizeof (int));
  file->bboxes = calloc ((size_t) (2 * file->maxndims * file->max_level),
                         sizeof (int));
  for (i = 0; i < file->npatches; i++)
  {
    const patch_t *patch = &file->patch[i];
    for (j = 0; j < (size_t) patch->ndims; j++)
    {
      maxthickness = patch->iorigin[j] + patch->dims[j] * patch->factor;
      if (file->maxthickness[j] < maxthickness)
      {
        file->maxthickness[j] = maxthickness;
      }
    }
  }

  DXMessage ("file '%s' contains %d timestep(s) with %d patch(es) in total",
             filename, file->ntimesteps, file->npatches);
  if (file->ntimesteps > 0)
  {
    msg = malloc (50 + 7 * (size_t) file->timesteps[0].ndims);
    sprintf (msg, "grid dimensions are (%d",
             file->maxthickness[file->maxndims-1]);
    for (l = file->maxndims-2; l >= 0; l--)
    {
      sprintf (msg, "%s x %d", msg, file->maxthickness[l]);
    }
    DXMessage ("%s) with respect to finest refinement level %d ",
               msg, file->max_level);
#if 0
    for (l = 0; l < file->max_level; l++)
    {
      printf ("level %d: ", l);
      for (m = 0; m < file->maxndims; m++)
      {
        printf ("[%d, %d]", file->bboxes[2*(l*file->maxndims+m) + 0],
                            file->bboxes[2*(l*file->maxndims+m) + 1]);
      }
    }
#endif
    free (msg);
  }

  return (OK);
}


 /*@@
   @routine    ExamineFile
   @date       Sun 1 June 2003
   @author     Thomas Radke
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
  size_t typesize;
  hid_t dataset, datatype, dataspace, attr;
  H5T_class_t typeclass, member1typeclass = -1, member0typeclass = -1;
  H5G_stat_t object_info;
  size_t new_size;
  char *objectname;
  patch_t *this;
  file_t *file = _file;


#ifdef DEBUG
  DXMessage ("Examining object '%s'", name);
#endif

  /* build the full name for the current object to process */
  objectname = file->objectname;
  file->objectname = malloc (strlen (objectname) + strlen (name) + 2);
  sprintf (file->objectname, "%s/%s", objectname, name);

  /* we are interested in datasets only - skip anything else */
  CHECK_ERROR (H5Gget_objinfo (group, name, 0, &object_info));
  if (object_info.type != H5G_DATASET)
  {
    if (object_info.type == H5G_GROUP)
    {
      /* iterate over all datasets in this group (if it isn't file metadata) */
      if (strcmp (name, FILE_METADATA_GROUP))
      {
        CHECK_ERROR (H5Giterate (group, file->objectname,NULL,ExamineFile,file))
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

  /* check whether the patch id array needs to be extended */
  if (file->npatches >= file->allocated)
  {
    file->allocated += H5_BLOCK_INCREMENT;
    new_size = file->allocated * sizeof (patch_t);
    if (file->patch)
    {
      file->patch = realloc (file->patch, new_size);
    }
    else
    {
      file->patch = malloc (new_size);
    }
  }

  this = &file->patch[file->npatches];
  memset (this, 0, sizeof (*this));
  this->name = strdup (file->objectname);
  this->typesize = typesize;
  this->typeclass = typeclass;
  this->fid = file->current_fid;

  /* read the dimensions */
  CHECK_ERROR (dataspace = H5Dget_space (dataset));
  this->ndims = H5Sget_simple_extent_ndims (dataspace);
  this->iorigin = malloc (this->ndims * sizeof (int));
  this->origin = malloc (2 * this->ndims * sizeof (float));
  this->delta = this->origin + this->ndims;
  this->start = malloc (this->ndims * sizeof (hssize_t));
  this->dims = malloc (3 * this->ndims * sizeof (hsize_t));
  this->count  = this->dims + 1*this->ndims;
  this->stride = this->count + 2*this->ndims;
  CHECK_ERROR (H5Sget_simple_extent_dims (dataspace, this->dims, NULL));
  CHECK_ERROR (H5Sclose (dataspace));

  /* read attributes */
  CHECK_ERROR (attr = H5Aopen_name (dataset, "origin"));
  CHECK_ERROR (H5Aread (attr, H5T_NATIVE_FLOAT, this->origin));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (attr = H5Aopen_name (dataset, "delta"));
  CHECK_ERROR (H5Aread (attr, H5T_NATIVE_FLOAT, this->delta));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (attr = H5Aopen_name (dataset, "time"));
  CHECK_ERROR (H5Aread (attr, H5T_NATIVE_DOUBLE, &this->time));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (attr = H5Aopen_name (dataset, "timestep"));
  CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &this->timestep));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (attr = H5Aopen_name (dataset, "level"));
  CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, &this->level));
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (attr = H5Aopen_name (dataset, "iorigin"));
  CHECK_ERROR (H5Aread (attr, H5T_NATIVE_INT, this->iorigin));
  CHECK_ERROR (H5Aclose (attr));

  CHECK_ERROR (H5Dclose (dataset));
  free (file->objectname);
  file->objectname = objectname;

  /* get the maximum number of dimensions over all datasets */
  if (file->maxndims < this->ndims)
  {
    file->maxndims = this->ndims;
  }

  /* increment patch counter */
  file->npatches++;

  return (0);
}


 /*@@
   @routine    ReadParameters
   @date       Sun 1 June 2003
   @author     Thomas Radke
   @desc
               Reads the parameters from the input tabs
   @enddesc

   @var        in
   @vdesc      list of input tabs
   @vtype      const Object *
   @vio        in
   @endvar
   @var        regex
   @vdesc      regular expression string to match against variable names
   @vtype      const char *
   @vio        in
   @endvar
   @var        preg
   @vdesc      compiled regular expression for matching variable names
   @vtype      const regex_t *
   @vio        in
   @endvar
   @var        file
   @vdesc      info structure describing the HDF5 input file
   @vtype      const file_t *
   @vio        in
   @endvar

   @returntype timestep_t *
   @returndesc
               pointer to the structure identifying the timestep to read,
               or NULL if there was an invalid input parameter setting
   @endreturndesc
@@*/
static timestep_t *ReadParameters (const Object *in,
                                   const char *regex, const regex_t *preg,
                                   const file_t *file)
{
  Type type;
  int i, j, level, active_levels, ndims, nmatches, rank, timestep;
  const int *vector;
  timestep_t *this;


  /* read the timestep index and check that it is valid */
  timestep = 0;
  if (in[5])
  {
    DXExtractInteger (in[5], &timestep);
  }
  this = timestep >= 0 && timestep < (int) file->ntimesteps ?
         &file->timesteps[timestep] : NULL;
  if (! this)
  {
    if (file->ntimesteps > 0)
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "Invalid timestep index %d (must be in range [0, %d)",
                  timestep, file->ntimesteps);
    }
    else
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "There are no datasets in file '%s'", file->last_filename);
    }
    return (NULL);
  }

  /* check the 'varnames' parameter */
  nmatches = 0;
  for (i = 0; i < (int) this->npatches; i++) {
    if (regexec (preg, this->patch[i].name, 0, NULL, 0) == 0) nmatches++;
  }
  if (nmatches == 0) {
    DXSetError (ERROR_BAD_PARAMETER,
                "This timestep doesn't contain any datasets with names "
                "matching '%s'", regex);
    return (NULL);
  }

  /* read the active levels */
  this->active_levels = -1;
  if (in[4])
  {
    DXGetArrayInfo ((Array) in[4], &ndims, &type, NULL, &rank, NULL);
    if (type != TYPE_INT || (rank != 0 && rank != 1))
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "levels must be given as integer list or integer vector");
      return (NULL);
    }
    if (rank == 1)
    {
      DXGetArrayInfo ((Array) in[4], NULL, NULL, NULL, NULL, &ndims);
    }
    vector = DXGetArrayData ((Array) in[4]);
    this->active_levels = 0;
    for (i = 0; i < ndims; i++)
    {
      if (vector[i] < 0)
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "levels[%d] = %d is out of range (must be in range [0, ",
                    "%d])", i, vector[i], file->max_level);
        return (NULL);
      }
      if (vector[i] > file->max_level)
      {
        DXWarning ("levels[%d] = %d specifies a non-existing level (finest "
                   "level is %d) - specified level will be ignored",
                   i, vector[i], file->max_level);
      }
      this->active_levels |= 1 << vector[i];
    }
    if (! this->active_levels)
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "levels parameter didn't specify any valid level");
      return (NULL);
    }
  }

  /* read the hyperslab origin */
  memset (this->origin, 0, this->ndims * sizeof (*this->origin));
  if (in[1])
  {
    DXGetArrayInfo ((Array) in[1], &ndims, &type, NULL, &rank, NULL);
    if (type != TYPE_INT || (rank != 0 && rank != 1))
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "origin must be given as integer list or integer vector");
      return (NULL);
    }
    if (rank == 1)
    {
      DXGetArrayInfo ((Array) in[1], NULL, NULL, NULL, NULL, &ndims);
    }
    if (ndims > this->ndims)
    {
      DXWarning ("origin has %d elements but selected timestep %d has only "
                 "%d dimensions (exceeding elements will be ignored)",
                 ndims, timestep, this->ndims);
      ndims = this->ndims;
    }
    vector = DXGetArrayData ((Array) in[1]);
    for (i = 0; i < ndims; i++)
    {
      if (vector[i] < 0 || vector[i] >= file->maxthickness[this->ndims-i-1])
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "origin[%d] = %d is out of range (must be in range [0, "
                    "%d])", i, vector[i], file->maxthickness[this->ndims-i-1] - 1);
        return (NULL);
      }
      this->origin[this->ndims-i-1] = vector[i];
    }
  }

  /* read the hyperslab extent (thickness) */
  for (i = 0; i < this->ndims; i++)
  {
    this->thickness[i] = file->maxthickness[i] - this->origin[i];
  }
  if (in[2])
  {
    DXGetArrayInfo ((Array) in[2], &ndims, &type, NULL, &rank, NULL);
    if (type != TYPE_INT || (rank != 0 && rank != 1))
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "thickness must be given as integer list or integer vector");
      return (NULL);
    }
    if (rank == 1)
    {
      DXGetArrayInfo ((Array) in[2], NULL, NULL, NULL, NULL, &ndims);
    }
    if (ndims > this->ndims)
    {
      DXWarning ("thickness vector has %d elements but selected timestep %d "
                 "has only %d dimensions (exceeding elements will be ignored)",
                 ndims, timestep, this->ndims);
      ndims = this->ndims;
    }
    vector = DXGetArrayData ((Array) in[2]);
    for (i = 0; i < ndims; i++)
    {
      if (vector[i] < 0 || vector[i] > file->maxthickness[this->ndims-i-1])
      {
        DXSetError (ERROR_BAD_PARAMETER,
                    "thickness[%d] = %d is out of range (must be in range [0, "
                    "%d])", i, vector[i], file->maxthickness[this->ndims-i-1]);
        return (NULL);
      }
      if (vector[i])
      {
        this->thickness[this->ndims-i-1] = vector[i];
      }
    }
  }

  /* verify that start + extent are still within the timestep dimensions */
  for (i = 0; i < this->ndims; i++)
  {
    if (this->origin[i] + this->thickness[i] > file->maxthickness[i])
    {
      DXSetError (ERROR_BAD_PARAMETER,
                  "origin (%d) + thickness (%d) in dimension %d is out of "
                  "range (must be in range [1, %d])", this->origin[i],
                  this->thickness[i], this->ndims-i-1, file->maxthickness[i]);
      return (NULL);
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
      return (NULL);
    }
    if (rank == 1)
    {
      DXGetArrayInfo ((Array) in[3], NULL, NULL, NULL, NULL, &ndims);
    }
    if (ndims > this->ndims)
    {
      DXWarning ("stride vector has %d elements but selected timestep %d has "
                 "only %d dimensions (exceeding elements will be ignored)",
                 ndims, timestep, this->ndims);
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
        return (NULL);
      }
      this->stride[this->ndims-i-1] = vector[i];
      if (vector[i] != 1)
      {
        /* count the number of active levels */
        for (j = active_levels = 0, level = -1; j < (int) this->npatches; j++)
        {
          if (this->active_levels & (1 << this->patch[j].level))
          {
            if (level < 0)
            {
              level = this->patch[j].level;
            }
            if (level != this->patch[j].level)
            {
              active_levels++;
            }
          }
        }
        if (active_levels != 0)
        {
          DXWarning ("Striding for more than one active level isn't "
                     "implemented yet, setting of stride[%d] = %d will be "
                     "ignored", i, vector[i]);
          this->stride[this->ndims-i-1] = 1;
        }
      }
    }
  }

  /* read the precision flag */
  this->single_precision = 1;
  if (in[7])
  {
    DXExtractInteger (in[7], &this->single_precision);
  }

   /* read the invalidation flag */
  this->invalidate_fine_grid = 1;
  if (in[13])
  {
    DXExtractInteger (in[13], &this->invalidate_fine_grid);
  }

  return (this);
}


 /*@@
   @routine    ReadPatch
   @date       Sun 1 June 2003
   @author     Thomas Radke
   @desc
               Reads a single timestep from the given HDF5 file.
   @enddesc

   @var        timestep
   @vdesc      structure to identify the timestep to read
   @vtype      timestep_t *
   @vio        in
   @endvar
   @var        this
   @vdesc      structure to identify the patch to read
   @vtype      patch_t *
   @vio        inout
   @endvar
   @var        patch
   @vdesc      return structure with data and bbox fields
               The bbox field object should always exist, the data field object
               is created only if patch is overlapping with the requested slab
               (otherwise it will be NULL).
   @vtype      DataBBox_t *
   @vio        out
   @endvar
@@*/
static void ReadPatch (timestep_t *timestep, patch_t *this, DataBBox_t *patch)
{
  Type dxtype;
  Category dxcategory;
  Array array;
  hid_t dataset, datatype, slabspace, filespace, plist;
  int i, j;
  int *dims;
  float *forigin, *fdelta;
  Pointer data;
  float positions[(1 << 3) * 3];
  RGBColor color = DXRGB (0.7, 0.7, 0.0), delta = DXRGB (0.0, 0.0, 0.0);
  const Line connections[] = {{0, 1}, {0, 2}, {1, 3}, {2, 3},
                              {4, 5}, {4, 6}, {5, 7}, {6, 7},
                              {0, 4}, {1, 5}, {2, 6}, {3, 7}};


  patch->data = patch->bbox = NULL;

  /* calculate bbox */
  for (i = j = 0; i < (1 << this->ndims); i++)
  {
    for (j = 0; j < this->ndims; j++)
    {
      positions[i*this->ndims + j] = this->origin[j];
      if ((i / (1 << (this->ndims-1-j))) & 1)
      {
        positions[i*this->ndims + j] += (this->dims[this->ndims-1-j]-1) *
                                        this->delta[j];
      }
    }
  }

  /* create bbox field object */
  patch->bbox = DXNewField ();
  if (! patch->bbox)
  {
    return;
  }
  /* set bbox positions */
  array = DXNewArrayV (TYPE_FLOAT, CATEGORY_REAL, 1, &j);
  DXAddArrayData (array, 0, i, NULL);
  memcpy (DXGetArrayData (array), positions, i*j * sizeof (float));
  DXSetComponentValue (patch->bbox, "positions", (Object) array);

  /* set bbox connections */
  array = DXNewArray (TYPE_INT, CATEGORY_REAL, 1, 2);
  DXAddArrayData (array, 0, i*j/2, NULL);
  memcpy (DXGetArrayData (array), connections, sizeof (connections));
  DXSetStringAttribute ((Object) array, "element type", "lines");
  DXSetComponentValue (patch->bbox, "connections", (Object) array);

  /* set bbox default color, using regular (compact) array */
  array = (Array) DXNewRegularArray (TYPE_FLOAT, j, i, &color, &delta);
  DXSetStringAttribute ((Object) array, "dep", "positions");
  DXSetComponentValue (patch->bbox, "colors", (Object) array);

  /* finish bbox field object */
  DXEndField (patch->bbox);

  for (i = 0; i < this->ndims; i++)
  {
    if (this->iorigin[i] >= timestep->origin[i] +
                            (hssize_t) timestep->thickness[i] ||
        this->iorigin[i] + (hssize_t) this->dims[i]*this->factor <=
        timestep->origin[i])
    {
#if 0
fprintf (stderr, "timestep %d patch %d outside: %d >= %d + %d  ||  %d + %d*%d <= %d\n",
timestep->timestep, this - timestep->patch,
(int) this->iorigin[i], (int) timestep->origin[i], (int) timestep->thickness[i],
(int) this->iorigin[i], (int) this->dims[i], (int) this->factor, (int) timestep->origin[i]);
#endif
      return;
    }

    this->start[i] = 0;
    if (timestep->origin[i] > this->iorigin[i])
    {
      this->start[i] = (timestep->origin[i] - this->iorigin[i] +
                        (this->factor >> 1)) / this->factor;
      if (this->start[i] >= (hssize_t) this->dims[i])
      {
        this->start[i] = this->dims[i] - 1;
      }
    }

    this->count[i] = this->dims[i] - this->start[i];
    if (timestep->origin[i] + timestep->thickness[i] <
        this->iorigin[i] + (int) this->dims[i]*this->factor)
    {
      this->count[i] = (timestep->origin[i] + timestep->thickness[i] -
                        this->iorigin[i] - this->start[i]*this->factor +
                        this->factor-1) / this->factor;
#if 0
DXMessage ("start[%d] %d = (%d - %d + %d) / %d\n", i, (int) this->start[i], (int) timestep->origin[i], (int) this->iorigin[i], (int) (this->factor >> 1), (int) this->factor);
DXMessage ("count[%d] %d = (%d + %d - %d - %d*%d + %d) / %d\n", i, (int) this->count[i], (int) timestep->origin[i], (int) timestep->thickness[i], (int) this->iorigin[i], (int) this->start[i], (int) this->factor, (int) this->factor-1, (int) this->factor);
#endif
    }
    if (this->count[i] == 0)
    {
      this->count[i] = 1;
    }
    else if (this->count[i] > (hsize_t) timestep->thickness[i])
    {
      this->count[i] = timestep->thickness[i];
    }
  }

  /* determine the datatype to use */
  if (this->typeclass == H5T_FLOAT)
  {
    if (this->typesize == sizeof (float) || timestep->single_precision)
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
    if (this->typesize == 2*sizeof (float) || timestep->single_precision)
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
    return;
  }

  /* allocate some temporary buffers */
  dims    = malloc (this->ndims * sizeof (int));
  forigin = calloc ((size_t) (this->ndims+1) * this->ndims, sizeof (float));
  fdelta  = forigin + timestep->ndims;

  /* by default: select full HDF5 dataspace, copy DX dims, origin, delta */
  for (i = 0; i < this->ndims; i++)
  {
    this->count[i] = (this->count[i] + timestep->stride[i] - 1) /
                     timestep->stride[i];
    dims[i] = this->count[i];
    fdelta[(this->ndims-i-1)*this->ndims + i] =
      this->delta[i] * timestep->stride[this->ndims-i-1];
    forigin[i] = this->origin[i] +
                 this->start[this->ndims-i-1] * this->delta[i];
  }

  CHECK_ERROR (slabspace = H5Screate_simple (this->ndims, this->count, NULL));
  i = H5Sget_simple_extent_npoints (slabspace);
  data = DXGetArrayData (DXAddArrayData (array, 0, i, NULL));
  if (data)
  {
    CHECK_ERROR (filespace = H5Screate_simple (this->ndims, this->dims, NULL));
    CHECK_ERROR (H5Sselect_hyperslab (filespace, H5S_SELECT_SET, this->start,
                                      timestep->stride, this->count, NULL));

    /* set the buffer size to read the full hyperslab */
    CHECK_ERROR (plist = H5Pcreate (H5P_DATASET_XFER));
    CHECK_ERROR (H5Pset_buffer (plist, i * this->typesize, NULL, NULL));

    CHECK_ERROR (dataset = H5Dopen (this->fid, this->name, H5P_DEFAULT));
    CHECK_ERROR (H5Dread (dataset, datatype, slabspace, filespace, plist,data));
    CHECK_ERROR (H5Dclose (dataset));
    CHECK_ERROR (H5Pclose (plist));
    CHECK_ERROR (H5Sclose (filespace));

    /* create data field object */
    patch->data = DXNewField ();
  }
  else
  {
    DXSetError (ERROR_NO_MEMORY, "Couldn't allocate memory for new DX array");
    patch->data = NULL;
  }

  /* close the dataspace and datatype */
  CHECK_ERROR (H5Sclose (slabspace));
  if (datatype != H5T_NATIVE_INT && datatype != H5T_NATIVE_FLOAT &&
      datatype != H5T_NATIVE_DOUBLE)
  {
    CHECK_ERROR (H5Tclose (datatype));
  }

  if (patch->data)
  {
    /* add all attributes of this dataset to the new field */
    AddAttributes (this, (Object) patch->data);

    /* set the individual components in the field */
    DXSetComponentValue (patch->data, "data", (Object) array);
    array = DXMakeGridPositionsV (this->ndims, dims, forigin, fdelta);
    DXSetComponentValue (patch->data, "positions", (Object) array);

    /* eliminate dimensions along which to slab */
    for (i = j = 0; i < this->ndims; i++)
    {
      if (this->count[i] > 1)
      {
        dims[j++] = this->count[i];
      }
    }
    array = DXMakeGridConnectionsV (j, dims);
    DXSetComponentValue (patch->data, "connections", (Object) array);
    DXEndField (patch->data);
  }
  else
  {
    DXDelete ((Object) array);
  }

  /* free temporary arrays */
  free (forigin);
  free (dims);
}


 /*@@
   @routine    AddAttributes
   @date       Sun 1 June 2003
   @author     Thomas Radke
   @desc
               Reads all attributes from the given dataset and adds them
               to the given DX object.
   @enddesc

   @var        this
   @vdesc      structure to identify the timestep to read
   @vtype      const patch_t *
   @vio        in
   @endvar
   @var        object
   @vdesc      DX object to attach attributes
   @vtype      Object
   @vio        inout
   @endvar
@@*/
static void AddAttributes (const patch_t *this, Object object)
{
  int i, nelems;
  void *value;
  Object dxattr;
  char attrname[256];
  H5T_class_t typeclass;
  hid_t dataset, datatype, dataspace, attr;


  CHECK_ERROR (dataset = H5Dopen (this->fid, this->name, H5P_DEFAULT));
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


 /*@@
   @routine    InvalidateOverlappingBoundaries
   @date       Sun 1 June 2003
   @author     Thomas Radke
   @desc
               Loops through all patches in the given multigrid group
               and invalidates points on coarser patches which overlap
               which refined patches.
   @enddesc

   @var        this
   @vdesc      structure to identify the timestep to read
   @vtype      const patch_t *
   @vio        in
   @endvar
   @var        file
   @vdesc      structure describing the datafile in general
   @vtype      const file_t *
   @vio        in
   @endvar
   @var        group
   @vdesc      DX multigrid group with individual patches
   @vtype      Group
   @vio        inout
   @endvar
@@*/
static void InvalidateOverlappingBoundaries (const timestep_t *this,
                                             const file_t *file,
                                             Group group)
{
  size_t i, j;
  int idx, nextlevel, num_invalidated;
  int point[3], frange[6], crange[6];
  hsize_t ii, jj, kk;
  Object cfield, ffield;
  patch_t *c, *f;
  char patchname[32];
  InvalidComponentHandle handle;


  cfield = NULL;

  /*** FIXME: limited to 2- and 3-dimensional patches ***/
  if (this->ndims < 2 || this->ndims > 3)
  {
    DXWarning ("removal of coarse-level grids overlapping with finer grids "
               "has been implemented for 2- and 3-dimensional grids only");
    return;
  }

  /* loop over all patches */
  int first_patch = 1;
  for (i = 0; i < this->npatches; i++)
  {
    /* skip inactive levels */
    if (! (this->active_levels & (1 << this->patch[i].level)))
    {
      continue;
    }

    for (nextlevel = this->patch[i].level + 1;
         ! ((1 << nextlevel) & this->active_levels) &&
            nextlevel <= file->max_level;
         nextlevel++);
    
    

    /* loop over all following patches */
    for (j = i + 1; j < this->npatches; j++)
    {
      /* skip patches which aren't on the next or on this level */
      if (this->patch[j].level != nextlevel && this->patch[j].level != this->patch[i].level)
      {
        continue;
      }
      
      /* check that both patches really are members of the group */
      sprintf (patchname, "patch %d", i);
      cfield = DXGetMember (group, patchname);
      sprintf (patchname, "patch %d", j);
      ffield = DXGetMember (group, patchname);
      if (! cfield || ! ffield)
      {
        continue;
      }

      c = this->patch + i;
      f = this->patch + j;

#ifdef DEBUG
      DXMessage ("patch %d of level %d potentially overlaps with patch %d "
                 "of level %d\n", i, c->level, j, f->level);
#endif

      /* set the start and end vectors to exclude the first and last point */
      if (!f->range)
      {
         f->range = malloc (2 * this->ndims * sizeof (int));
      }
      
      for (ii = 0; ii < (hsize_t) this->ndims; ii++)
      {
      
         /* ranges on the finer grid (boundary points inclusive) */
         if (f->level != c->level)
         {
            frange[2*ii+0] = f->start[ii] * f->factor + f->iorigin[ii];
            frange[2*ii+1] = frange[2*ii+0] + (f->count[ii]-1) * f->factor;
#if 0
DXMessage ("dim %d: starting with lower range %d = %d*%d + %d and upper range %d = %d + %d*%d", (int) ii,
           frange[2*ii+0], (int) f->start[ii], (int) f->factor, (int) f->iorigin[ii], (int) frange[2*ii+1], (int) frange[2*ii+0], (int) f->count[ii], (int) f->factor);
#endif

            /* ranges on the coarser grid (boundary points exclusive) */
            crange[2*ii+0] = frange[2*ii+1] + 1;
            crange[2*ii+1] = frange[2*ii+0] - 1;
            if (frange[2*ii+0] % c->factor)
            {
               if (this->thickness[ii] > 1)
               {
                  crange[2*ii+0] = frange[2*ii+0];
               }
               frange[2*ii+0] += c->factor - frange[2*ii+0] % c->factor;
            }
            if (frange[2*ii+1] % c->factor)
            {
               if (this->thickness[ii] > 1)
               {
                  crange[2*ii+1] = frange[2*ii+1];
               }
               frange[2*ii+1] -= frange[2*ii+1] % c->factor;
            }
            
            f->range[2*ii+0] = frange[2*ii+0]; 
            f->range[2*ii+1] = frange[2*ii+1]; 
         }
         else
	 {
	    if (first_patch)
            {
	       f->range[2*ii+0] = f->start[ii] * f->factor + f->iorigin[ii];
	       f->range[2*ii+1] = f->range[2*ii+0] + (f->count[ii]-1) * f->factor;
	    }
	
            frange[2*ii+0] = f->range[2*ii+0]; 
            frange[2*ii+1] = f->range[2*ii+1];
            
            crange[2*ii+0] = frange[2*ii+1] + 1;
            crange[2*ii+1] = frange[2*ii+0] - 1;
         }
      }

      /* invalidate all positions on the coarser patch
         which are already covered by positions on the finer grid */
      handle = DXCreateInvalidComponentHandle (cfield, NULL, "positions");
      idx = num_invalidated = 0;
      point[0] = c->start[0]*c->factor + c->iorigin[0];
#if 0
DXMessage ("level %d: disabling points %d-%d and %d-%d and %d-%d", c->level,
           frange[4], frange[5], frange[2], frange[3], frange[0], frange[1]);
#endif
      for (kk = 0; kk < c->count[0]; kk++, point[0] += c->factor)
      {
        point[1] = c->start[1]*c->factor + c->iorigin[1];
        for (jj = 0; jj < c->count[1]; jj++, point[1] += c->factor)
        {
          if (this->ndims == 2)
          {
            if ((this->thickness[0] == 1 ||
                 (frange[0] < point[0] && point[0] < frange[1])) &&
                (this->thickness[1] == 1 ||
                 (frange[2] < point[1] && point[1] < frange[3])))
            {
              DXSetElementInvalid (handle, idx);
              num_invalidated++;
            }
            idx++;
          }
          else
          {
            point[2] = c->start[2]*c->factor + c->iorigin[2];
            for (ii = 0; ii < c->count[2]; ii++, point[2] += c->factor, idx++)
            {
              if ((this->thickness[0] == 1 ||
                   (frange[0] < point[0] && point[0] < frange[1])) &&
                  (this->thickness[1] == 1 ||
                   (frange[2] < point[1] && point[1] < frange[3])) &&
                  (this->thickness[2] == 1 ||
                   (frange[4] < point[2] && point[2] < frange[5])))
              {
                DXSetElementInvalid (handle, idx);
                num_invalidated++;
              }
            }
          }
        }
      }
      if (num_invalidated)
      {
        DXSaveInvalidComponent ((Field) cfield, handle);
      }
      DXFreeInvalidComponentHandle (handle);

      

      /* invalidate overlapping boundaries on the finer patch */
      handle = DXCreateInvalidComponentHandle (ffield, NULL, "positions");
      idx = num_invalidated = 0;
      point[0] = f->start[0]*f->factor + f->iorigin[0];
      for (kk = 0; kk < f->count[0]; kk++, point[0] += f->factor)
      {
        point[1] = f->start[1]*f->factor + f->iorigin[1];
        for (jj = 0; jj < f->count[1]; jj++, point[1] += f->factor)
        {
          if (this->ndims == 2)
          {
            if ((crange[0] <= point[0] && point[0] <  frange[0]) ||
                (frange[1] <  point[0] && point[0] <= crange[1]) ||
                (crange[2] <= point[1] && point[1] <  frange[2]) ||
                (frange[3] <  point[1] && point[1] <= crange[3]))
            {
              DXSetElementInvalid (handle, idx);
              num_invalidated++;
            }
            idx++;
          }
          else
          {
            point[2] = f->start[2]*f->factor + f->iorigin[2];
            for (ii = 0; ii < f->count[2]; ii++, point[2] += f->factor, idx++)
            {
              if ((crange[0] <= point[0] && point[0] <  frange[0]) ||
                  (frange[1] <  point[0] && point[0] <= crange[1]) ||
                  (crange[2] <= point[1] && point[1] <  frange[2]) ||
                  (frange[3] <  point[1] && point[1] <= crange[3]) ||
                  (crange[4] <= point[2] && point[2] <  frange[4]) ||
                  (frange[5] <  point[2] && point[2] <= crange[5]))
              {
                DXSetElementInvalid (handle, idx);
                num_invalidated++;
              }
            }
          }
        }
      }
      if (num_invalidated && this->invalidate_fine_grid)
      {
        DXSaveInvalidComponent ((Field) ffield, handle);
      }
      DXFreeInvalidComponentHandle (handle);
    }

    /* also invalidate connections now (which depend on the positions) */
    DXInvalidateConnections (cfield);
#if 0
    DXCull (cfield);
    DXChangedComponentStructure ((Field) cfield, "positions");
#endif

    first_patch = 0;
  }
}


/* comparison function used by qsort(3) to sort patches by time and level */
static int ComparePatches (const void *_a, const void *_b)
{
  const patch_t *a = _a;
  const patch_t *b = _b;


  if (a->timestep == b->timestep)
  {
    return (a->level > b->level ? +1 : -1);
  }
  else
  {
    return (a->time > b->time ? +1 : -1);
  }
}
