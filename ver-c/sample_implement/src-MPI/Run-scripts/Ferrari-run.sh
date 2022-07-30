#!/bin/sh
#PBS -l nodes=2:ppn=16
#PBS -N VF-Ferrari
#PBS -o ../data_out/outdata
#PBS -e ../data_out/outdata-error
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
mpirun -machinefile $PBS_NODEFILE -np $NPROCS ./a.out
