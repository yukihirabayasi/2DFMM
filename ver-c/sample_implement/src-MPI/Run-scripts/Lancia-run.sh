#!/bin/sh
#PBS -l nodes=1:ppn=20
#PBS -N VM-Lancia
#PBS -o ../data_out/outdata
#PBS -e ../data_out/outdata-error
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
mpirun -np $NPROCS ./a.out
