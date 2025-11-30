# -*-makefile-*-

ARCH=LINUX
HOME=/home/stefan/LLNLfiles
# LINKER FLAGS GIVEN BELOW
LDFLAGS_LINUX= 
CC=g++
CLINKER=g++
LDFLAGS=$(LDFLAGS_$(ARCH))

# PPlusPlus is an environment variable that
# points to the directory above the lib/include directories
PPP_INCLUDE=-I/tmp/A++-P++-0.8.3/P++/install/include
PPPLIB_DIR=/tmp/A++-P++-0.8.3/P++/install/lib
PPP_INCLUDE_G=-I$(HOME)/A++-P++-0.8.3-Debug/P++/install/include
PPPLIB_DIR_G=/tmp/A++-P++-0.8.3-Debug/P++/install/lib

AZTEC_INCLUDE=-I$(HOME)/Aztec/lib
AZTEC_LIB=$(HOME)/Aztec/lib

# HDF5_DIR is an environment variable that
# points to the directory above the lib/include directories
HDF5_DIR=/usr/include/hdf5/openmpi
HDF5_INCLUDE=-I/usr/include/hdf5/openmpi/
HDF5_LIB=-L/usr/lib/x86_64-linux-gnu

LIBS_LINUX=-L$(PPPLIB_DIR) -L$(AZTEC_LIB) \
	  -lPpp -lPpp_static -laztec \
	  -L$(HDF5_LIB) -lhdf5_openmpi -lpthread /usr/lib/x86_64-linux-gnu/openblas-serial/libopenblas.so.0 -lm -lmpi /usr/lib/x86_64-linux-gnu/openblas-serial/libopenblas.so.0
LIBS_LINUX_G=-L$(PPPLIB_DIR_G) -L$(AZTEC_LIB) \
	  -lPpp -lPpp_static -laztec \
	  -L$(HDF5_LIB) -lhdf5_openmpi -lpthread /usr/lib/x86_64-linux-gnu/openblas-serial/libopenblas.so.0 -lm -lmpi /usr/lib/x86_64-linux-gnu/openblas-serial/libopenblas.so.0
# INCLUDE SEARCH PATH GIVEN BELOW

INCLUDES= $(PPP_INCLUDE_G) $(AZTEC_INCLUDE) $(HDF5_INCLUDE) -I/usr/lib/petscdir/petsc3.22/x86_64-linux-gnu-real/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/

CPPFLAGS_LINUX=${INCLUDES} -DHAVE_CONFIG_H -I. -DAZ_MPI
CFLAGS_LINUX=-fpermissive -std=c++17 

CPPFLAGS=$(CPPFLAGS_$(ARCH))
CFLAGS=$(CFLAGS_$(ARCH)) 

LIBRARIES= $(LIBS_$(ARCH)) -L/usr/lib/petscdir/petsc3.122/x86_64-linux-gnu-real/lib -lpetsc_real
LIBRARIES_G= $(LIBS_$(ARCH)_G) -L/usr/lib/petscdir/petsc3.22/x86_64-linux-gnu-real/lib -lpetsc_real

# SOURCES GIVEN BELOW

# These files are needed by all applications
srcs= GridFunc.cc CompGrid.cc OGEquation.cc

# OBJECTS GIVEN BELOW
objs= $(srcs:.cc=.o)

%.o: %.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(CPPFLAGS) -O2

objs_g= $(srcs:.cc=_g.o)

%_g.o: %.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(CPPFLAGS) -g

my_auxiliary_routines.o: my_auxiliary_routines.cc my_auxiliary_routines.hh

ins_sharp.o: ins_sharp.cc ins_sharp.hh

heatEq.o: heatEq.cc ins_sharp.hh

GridFunc.o: GridFunc.cc GridFunc.hh

OGEquation.o: OGEquation.cc OGEquation.hh

CompGrid.o: CompGrid.cc CompGrid.hh

my_auxiliary_routines_g.o: my_auxiliary_routines.cc my_auxiliary_routines.hh

ins_sharp_g.o: ins_sharp.cc ins_sharp.hh

poissonEq_g.o: poissonEq.cc

heatEq_g.o: heatEq.cc ins_sharp.hh

GridFunc_g.o: GridFunc.cc GridFunc.hh

OGEquation_g.o: OGEquation.cc OGEquation.hh

CompGrid_g.o: CompGrid.cc CompGrid.hh

ins_sharp: ins_sharp.o $(objs)
	${CLINKER} -o $@ $(LDFLAGS) ins_sharp.o $(objs) $(LIBRARIES)

ins_g: ins_g.o my_auxiliary_routines_g.o $(objs_g)
	${CLINKER} -o $@ $(LDFLAGS) ins_g.o my_auxiliary_routines_g.o $(objs_g) $(LIBRARIES_G)

ins_sharp_g: ins_sharp_g.o $(objs_g)
	${CLINKER} -o $@ $(LDFLAGS) ins_sharp_g.o $(objs_g) $(LIBRARIES_G)

poissonEq_g: poissonEq_g.o $(objs_g)
	${CLINKER} -o $@ $(LDFLAGS) poissonEq_g.o $(objs_g) $(LIBRARIES_G)

poissonEq: poissonEq.o $(objs)
	${CLINKER} -o $@ $(LDFLAGS) poissonEq.o $(objs) $(LIBRARIES)

heatEq: heatEq.o $(objs)
	${CLINKER} -o $@ $(LDFLAGS) heatEq.o $(objs) $(LIBRARIES)

heatEq_g: heatEq_g.o $(objs_g)
	${CLINKER} -o $@ $(LDFLAGS) heatEq_g.o $(objs_g) $(LIBRARIES_G)

clean:
	rm *.o heatEq heatEq_g ins ins_g ins_sharp ins_sharp_g
