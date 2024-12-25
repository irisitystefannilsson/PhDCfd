#!/usr/bin/bash

echo "CARTESIAN MESH"
# single proc
echo "NUMPROCS = 1"
for SIZE in 11 21 41 81 161 321
do
    echo ${SIZE}x${SIZE}
    echo "acg"
    mpirun -np 1 ../ins_sharp --fileName ../XcogGrids/square_convergence_test_${SIZE}x${SIZE}.acg --twilight --allNeumann --stoptime 1 2>&1 | grep Error
    echo "h5"
    mpirun -np 1 ../ins_sharp --fileName ../XcogGrids/square_convergence_test_${SIZE}x${SIZE}.h5 --twilight --allNeumann --stoptime 1 2>&1 | grep Error
done

# multi proc
echo "NUMPROCS = 4"
for SIZE in 11 21 41 81 161 321
do
    echo ${SIZE}x${SIZE}
    echo "acg"
    mpirun -np 4 ../ins_sharp --fileName ../XcogGrids/square_convergence_test_${SIZE}x${SIZE}.acg --twilight --allNeumann --stoptime 1 2>&1 | grep Error
    echo "h5"
    mpirun -np 4 ../ins_sharp --fileName ../XcogGrids/square_convergence_test_${SIZE}x${SIZE}.h5 --twilight --allNeumann --stoptime 1 2>&1 | grep Error
done

echo "COMPOSITE MESH"
# single proc
echo "NUMPROCS = 1"
for SIZE in 21 41 81 161 321
do
    echo ${SIZE}x${SIZE}
    echo "acg"
    mpirun -np 1 ../ins_sharp --fileName ../XcogGrids/cGrid_${SIZE}x${SIZE}.acg --twilight --allNeumann --stoptime 1 2>&1 | grep Error
    echo "h5"
    mpirun -np 1 ../ins_sharp --fileName ../XcogGrids/cGrid_${SIZE}x${SIZE}.h5 --twilight --allNeumann --stoptime 1 2>&1 | grep Error 
done

# multi proc
echo "NUMPROCS = 4"
for SIZE in 21 41 81 161 321
do
    echo ${SIZE}x${SIZE}
    echo "acg"
    mpirun -np 4 ../ins_sharp --fileName ../XcogGrids/cGrid_${SIZE}x${SIZE}.acg --twilight --allNeumann --stoptime 1 2>&1 | grep Error
    echo "h5"
    mpirun -np 4 ../ins_sharp --fileName ../XcogGrids/cGrid_${SIZE}x${SIZE}.h5 --twilight --allNeumann --stoptime 1 2>&1 | grep Error 
done
