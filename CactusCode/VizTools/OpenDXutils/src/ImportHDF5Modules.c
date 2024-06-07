/* Automatically generated - may need to edit! */

#include <dx/dx.h>
#include <dx/modflags.h>

#if defined(intelnt) || defined(WIN32)
#include <windows.h>
#endif

#if defined(__cplusplus)
extern "C" Error DXAddModule (char *, ...);
#else
extern Error DXAddModule (char *, ...);
#endif

#if defined(__cplusplus)
extern "C" Error m_ImportHDF5(Object*, Object*);
#endif
#if defined(intelnt) || defined(WIN32)
void FAR WINAPI DXEntry()
#else
  #if defined(__cplusplus)
    extern "C" void DXEntry()
  #else
    void DXEntry()
  #endif
#endif
{
    {
#if defined(__cplusplus)
        extern "C" Error m_ImportHDF5(Object *, Object *);
#else
        extern Error m_ImportHDF5(Object *, Object *);
#endif
        DXAddModule("ImportHDF5", m_ImportHDF5, 
            MODULE_ASYNC,
            12, "filename", "origin", "thickness", "stride", "index", "reopen[visible:0]", "single_precision[visible:0]", "user[visible:0]", "password[visible:0]", "subject[visible:0]", "num_streams[visible:0]", "vectorimport[visible:0]",
            2, "result", "max_index");
    }
    {
#if defined(__cplusplus)
        extern "C" Error m_ImportCactusHDF5(Object *, Object *);
#else
        extern Error m_ImportCactusHDF5(Object *, Object *);
#endif
        DXAddModule("ImportCactusHDF5", m_ImportCactusHDF5, 
            MODULE_ASYNC,
            11, "filename", "origin", "thickness", "stride", "index", "reopen[visible:0]", "single_precision[visible:0]", "user[visible:0]", "password[visible:0]", "subject[visible:0]", "num_streams[visible:0]",
            2, "result", "max_index");
    }
    {
#if defined(__cplusplus)
        extern "C" Error m_ImportCarpetHDF5(Object *, Object *);
#else
        extern Error m_ImportCarpetHDF5(Object *, Object *);
#endif
        DXAddModule("ImportCarpetHDF5", m_ImportCarpetHDF5, 
            MODULE_ASYNC,
            14, "filename", "origin", "thickness", "stride", "levels", "timestep", "reopen[visible:0]", "single_precision[visible:0]", "user[visible:0]", "password[visible:0]", "subject[visible:0]", "num_streams[visible:0]", "varnames", "invalidate_fine_grid[visible:0]",
            3, "result", "bboxes", "max_timestep");
    }
    {
#ifndef __cplusplus
        extern Error m_ImportAHFinderFile(Object *, Object *);
#endif
        DXAddModule("ImportAHFinderFile", m_ImportAHFinderFile, 0,
            1, "name",
            1, "result");
    }
}
