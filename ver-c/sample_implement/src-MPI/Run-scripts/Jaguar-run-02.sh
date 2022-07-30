#!/bin/sh
#PBS -l walltime=1000:00:00
#PBS -l nodes=1:ppn=16:Jaguar-02
#PBS -j oe
#PBS -N VF-Jaguar-02
#PBS -o outdata
#PBS -e outdata-error

cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
mpirun -np $NPROCS ./a.out
