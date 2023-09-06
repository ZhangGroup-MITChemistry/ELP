#!/bin/bash
#SBATCH --job-name=V10_AA              # Job name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=192:00:00

module purge
module load mpi/openmpi-4.0

./initram-v5.sh -f V10_100.gro -o V10_100_CHARMM.gro -to charmm36 -p topol.top -em 1000 -nb 5000 -np 16
