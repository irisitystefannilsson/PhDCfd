# -*-makefile-*-
# Compile using gmake ARCH=<arch> BOPT=<opt> <target>

ARCH             = LINUX
BOPT             = O_c++
HOME=/home/stefan/LLNLfiles
# LINKER FLAGS GIVEN BELOW
LDFLAGS_LINUX = 
CC = g++
CLINKER = g++
LDFLAGS = $(LDFLAGS_$(ARCH))

# PPlusPlus is an environment variable that
# points to the directory above the lib/include directories
PPP_INCLUDE = -I$(HOME)/A++-P++/P++/install/include
PPPLIB_DIR = $(HOME)/A++-P++/P++/install/lib
PPP_INCLUDE_G = -I$(HOME)/A++-P++-0.8.3-Debug/P++/install/include
PPPLIB_DIR_G = $(HOME)/A++-P++-0.8.3-Debug/P++/install/lib

AZTEC_INCLUDE = -I$(HOME)/Aztec/lib
AZTEC_LIB = $(HOME)/Aztec/lib

# HDF5_DIR is an environment variable that
# points to the directory above the lib/include directories
HDF5_DIR = /usr/include/hdf5/openmpi
HDF5_INCLUDE = -I/usr/include/hdf5/openmpi/
HDF5_LIB = -L/usr/lib/x86_64-linux-gnu

LIBS_LINUX = -L$(PPPLIB_DIR) -L$(AZTEC_LIB) \
	  -lPpp -lPpp_static -laztec \
	  -L$(HDF5_LIB) -lhdf5_openmpi -lpthread /usr/lib/x86_64-linux-gnu/openblas-serial/libopenblas.so.0 -lm -lmpi -lmpi_cxx /usr/lib/x86_64-linux-gnu/openblas-serial/libopenblas.so.0
LIBS_LINUX_G = -L$(PPPLIB_DIR_G) -L$(AZTEC_LIB) \
	  -lPpp -lPpp_static -laztec \
	  -L$(HDF5_LIB) -lhdf5_openmpi -lpthread /usr/lib/x86_64-linux-gnu/openblas-serial/libopenblas.so.0 -lm -lmpi -lmpi_cxx /usr/lib/x86_64-linux-gnu/openblas-serial/libopenblas.so.0
# INCLUDE SEARCH PATH GIVEN BELOW

INCLUDES= $(PPP_INCLUDE_G) $(AZTEC_INCLUDE) $(HDF5_INCLUDE) -I/usr/lib/petscdir/petsc3.18/x86_64-linux-gnu-real/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/

CPPFLAGS_LINUX = ${INCLUDES} -DHAVE_CONFIG_H -I. -DHAVE_CONFIG_H -I. -DAZ_MPI -fpermissive -std=c++1z
CFLAGS_LINUX = -DHAVE_CONFIG_H -I. -DAZ_MPI

CPPFLAGS = $(CPPFLAGS_$(ARCH)) -D$(ARCH)
CFLAGS = $(CFLAGS_$(ARCH)) 

#include ${PETSC_DIR}/bmake/common/base

LIBRARIES= $(LIBS_$(ARCH)) -L/usr/lib/petscdir/petsc3.18/x86_64-linux-gnu-real/lib -lpetsc_real
LIBRARIES_G= $(LIBS_$(ARCH)_G) -L/usr/lib/petscdir/petsc3.18/x86_64-linux-gnu-real/lib -lpetsc_real

MAKEDEPEND_LINUX=makedepend
MAKEDEPEND = $(MAKEDEPEND_$(ARCH))

# SOURCES GIVEN BELOW

# These files are needed by all applications
srcs= GridFunc.cc CompGrid.cc OGEquation.cc

# OBJECTS GIVEN BELOW
objs= $(srcs:.cc=.o)

objs_g= $(srcs:.cc=_g.o)

my_auxiliary_routines.o: my_auxiliary_routines.cc my_auxiliary_routines.hh

ins_sharp.o: ins_sharp.cc ins_sharp.hh

heatEq.o: heatEq.cc ins_sharp.hh

ins.o: ins.cc my_auxiliary_routines.hh

GridFunc.o: GridFunc.cc GridFunc.hh

OGEquation.o: OGEquation.cc OGEquation.hh

CompGrid.o: CompGrid.cc CompGrid.hh

my_auxiliary_routines_g.o: my_auxiliary_routines.cc my_auxiliary_routines.hh
	${CC} -g -c -o $@ $(CPPFLAGS) $(CFLAGS) ${PETSC_CCPPFLAGS} my_auxiliary_routines.cc
ins_sharp_g.o: ins_sharp.cc ins_sharp.hh
	${CC} -g -c -o $@ $(CPPFLAGS) $(CFLAGS) ${PETSC_CCPPFLAGS} ins_sharp.cc
heatEq_g.o: heatEq.cc ins_sharp.hh
	${CC} -g -c -o $@ $(CPPFLAGS) $(CFLAGS) ${PETSC_CCPPFLAGS} heatEq.cc
ins_g.o: ins.cc my_auxiliary_routines.hh
	${CC} -g -c -o $@ $(CPPFLAGS) $(CFLAGS) ${PETSC_CCPPFLAGS} ins.cc
GridFunc_g.o: GridFunc.cc GridFunc.hh
	${CC} -g -c -o $@ $(CPPFLAGS) $(CFLAGS) ${PETSC_CCPPFLAGS} GridFunc.cc
OGEquation_g.o: OGEquation.cc OGEquation.hh
	${CC} -g -c -o $@ $(CPPFLAGS) $(CFLAGS) ${PETSC_CCPPFLAGS} OGEquation.cc
CompGrid_g.o: CompGrid.cc CompGrid.hh
	${CC} -g -c -o $@ $(CPPFLAGS) $(CFLAGS) ${PETSC_CCPPFLAGS} CompGrid.cc
ins: ins.o my_auxiliary_routines.o $(objs)
	${CLINKER} -o $@ $(LDFLAGS) ins.o my_auxiliary_routines.o $(objs) $(LIBRARIES)

ins_sharp: ins_sharp.o $(objs)
	${CLINKER} -o $@ $(LDFLAGS) ins_sharp.o $(objs) $(LIBRARIES_G)

ins_g: ins_g.o my_auxiliary_routines_g.o $(objs_g)
	${CLINKER} -o $@ $(LDFLAGS) ins_g.o my_auxiliary_routines_g.o $(objs_g) $(LIBRARIES_G)

ins_sharp_g: ins_sharp_g.o $(objs_g)
	${CLINKER} -o $@ $(LDFLAGS) ins_sharp_g.o $(objs_g) $(LIBRARIES_G)

heatEq: heatEq.o $(objs)
	${CLINKER} -o $@ $(LDFLAGS) heatEq.o $(objs) $(LIBRARIES_G)

heatEq_g: heatEq_g.o $(objs_g)
	${CLINKER} -o $@ $(LDFLAGS) heatEq_g.o $(objs_g) $(LIBRARIES_G)
