#!/bin/bash
#SBATCH --job-name=channelsApp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=60
#SBATCH --ntasks=60
#SBATCH --output=%x-%j.Out
#SBATCH --partition=hpc

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/shared/home/azureuser/lib/petsc/arch-linux2-c-debug/lib
source /shared/home/azureuser/lib/openfoam/mto/OpenFOAM-6/etc/bashrc
/shared/home/azureuser/lib/openfoam/mto/ThirdParty-6/platforms/linux64Gcc/openmpi-1.10.7/bin/mpirun -np $SLURM_NTASKS MTO_ThermalFluid -parallel > runLOG
