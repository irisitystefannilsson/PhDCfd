# GNU Makefile to build the ImportHDF5 and ImportAHFinderHDF5 modules
# provided in this directory as OpenDX loadable modules
#
# $Header: /cactus/VizTools/OpenDXutils/src/Makefile,v 1.13 2005/04/26 16:32:15 tradke Exp $

# check that the environment points to the OpenDX and HDF5 installations
ifeq ($(DXROOT),)
  $(error "You must set DXROOT to point to your OpenDX installation !")
endif
ifeq ($(HDF5_DIR),)
  $(error "You must set HDF5_DIR to point to your HDF5 installation !")
endif

# find out the architecture by looking for a $(DXROOT)/lib_$(DXARCH) directory
ifeq ($(DXARCH),)
  DXARCH := $(shell ls $(DXROOT) | grep lib_ | cut -c5-)
endif

# include the DXARCH makefile
include $(DXROOT)/lib_$(DXARCH)/arch.mak

# find out what compile and link flags HDF5 itself needs
HDF5_CFLAGS := $(shell grep -s 'Compiler switches:' $(HDF5_DIR)/lib/libhdf5.settings | cut -f2-)
HDF5_LIBS   := $(shell grep -s 'Extra libraries:' $(HDF5_DIR)/lib/libhdf5.settings | cut -d ':' -f2-)

# now set the makefile variables for compiling and linking
build_loadable_modules = 'yes'
CFLAGS += $(DX_RTL_CFLAGS) $(DXABI) -I$(DXROOT)/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi 
CFLAGS += -I$(HDF5_DIR)/include $(HDF5_CFLAGS) -fPIC

CC	= $(DXCC)
LD	= $(SHARED_LINK)
# the DX_RTL_LDFLAGS for Darwin are broken;
# hardcode them here in the hope they are working
ifeq ($(build_loadable_modules), 'yes')
ifeq ($(DXARCH), macos)
  LDLIBS += -bundle -flat_namespace -bundle_loader $(DXROOT)/bin_$(DXARCH)/dxexec -Wl,-iDXEntry:_DXEntry -e DXEntry -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi
else
  LDLIBS += $(DX_RTL_LDFLAGS) $(DXABI)
endif
else
  LDLIBS += -L$(DXROOT)/lib_$(DXARCH) -lDX $(DX_X_LINK_LIBS)
endif
LDLIBS += -L$(HDF5_DIR)/lib -lhdf5 $(HDF5_LIBS)

# the object files to build
OBJS	= ImportHDF5.o ImportCactusHDF5.o ImportCarpetHDF5.o \
          ImportAHFinderFile.o ImportHDF5Modules.o

# the build rules
ImportHDF5:	$(OBJS)

# The step (below) needs to be performed with a
# different value for DXROOT, otherwise
# the file 'mdf2c.awk' won't be found.
# It is currently (Debian 12) found in
# /usr/share/dx/lib
ImportHDF5Modules.c:	ImportHDF5.mdf
ifeq ($(build_loadable_modules), 'yes')
	/usr/lib/dx/bin/mdf2c -m $^ > $@
else
	/usr/lib/dx/bin/mdf2c $^ > $@
endif

.PHONY: clean
clean:
	rm -f $(OBJS) ImportHDF5Modules.c ImportHDF5
