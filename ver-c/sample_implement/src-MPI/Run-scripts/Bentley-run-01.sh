#!/bin/sh
#PBS -l nodes=1:ppn=40:bentley-01
#PBS -N VF-Bentley-01
#PBS -o ../data_out/outdata
#PBS -e ../data_out/outdata-error
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
mpirun -machinefile $PBS_NODEFILE -np $NPROCS ./a.out
