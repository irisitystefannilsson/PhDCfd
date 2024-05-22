# PhDCfd
Code written as part of my PhD work in Cfd on overset grids 

## Overture
Using version 0.8.3
from
https://sourceforge.net/p/overtureframework/code/ci/master/tree/packages/

Build: 
	./configure --with-mpi-include=-I/usr/lib/x86_64-linux-gnu/openmpi/include --with-mpi-lib-dirs=-L/usr/lib/x86_64-linux-gnu/openmpi/lib --enable-PXX '--with-mpi-libs=-lmpi -lmpi_cxx' --prefix=/home/stefan/LLNLfiles/A++-P++-0.8.3 --without-PADRE
	make
	make install
	and
	./configure --with-mpi-include=-I/usr/lib/x86_64-linux-gnu/openmpi/include --with-mpi-lib-dirs=-L/usr/lib/x86_64-linux-gnu/openmpi/lib --enable-PXX --with-mpi-libs="-lmpi -lmpi_cxx" --prefix=/home/stefan/LLNLfiles/A++-P++-0.8.3--Debug --enable-CXX_OPT=-g --enable-C_OPT=-g --without-PADRE
	

