#!/bin/sh
#PBS -l nodes=2:ppn=16
#PBS -N VF-Jaguar
#PBS -o ../data_out/outdata
#PBS -e ../data_out/outdata-error
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
mpirun -np $NPROCS ./a.out
