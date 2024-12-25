# PhDCfd
Code written as part of my PhD work in Cfd on overset grids 

## MPI - OpenMPI
Version 4.1.4 from Debian 12.5

## A++/P++
Using version 0.8.3
from
https://sourceforge.net/p/overtureframework/code/ci/master/tree/packages/

Build: 
	./configure --with-mpi-include=-I/usr/lib/x86_64-linux-gnu/openmpi/include --with-mpi-lib-dirs=-L/usr/lib/x86_64-linux-gnu/openmpi/lib --enable-PXX '--with-mpi-libs=-lmpi -lmpi_cxx' --prefix=/home/stefan/LLNLfiles/A++-P++-0.8.3 --without-PADRE
	make
	make install
	and
	./configure --with-mpi-include=-I/usr/lib/x86_64-linux-gnu/openmpi/include --with-mpi-lib-dirs=-L/usr/lib/x86_64-linux-gnu/openmpi/lib --enable-PXX --with-mpi-libs="-lmpi -lmpi_cxx" --prefix=/home/stefan/LLNLfiles/A++-P++-0.8.3--Debug --enable-CXX_OPT=-g --enable-C_OPT=-g --without-PADRE
	make
	make install
	(for debug build)

## Aztec
Using 2.0. Old download. Should move to Trilinos?

## Petsc
Using 3.18 from Debian 12.5

## Hdf5
Version 1.10.8 from Debian 12.5

## Xcog
Used to generate Overset grids.

Version 2.3 used. Local tar.gz-file.

### Hdf4
Xcog can save to hdf4 format binary files.

Version 4.2.15 from Debian 12.5

### H4toH5tools
Used to convert hdf4 xcog grids to hdf5 format.

Version 2.2.5 
from
https://portal.hdfgroup.org/downloads/h4h5tools/h4h5tools_2_2_5.html
Build:
	CPPFLAGS="-I/usr/include/hdf5/openmpi/ -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/home/stefan/include" LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/" ./configure --prefix=/home/stefan
	make (needed to add hdf5 libs to some makefiles for linking)
	make install


