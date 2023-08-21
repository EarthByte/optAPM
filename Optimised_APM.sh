#!/bin/bash

# Name.
#PBS -N optAPM

# Project.
#PBS -P q97

# Queue.
#PBS -q normal

# CPUs.
#
# Ideally the number of seed models (eg, 400) divided by this number
# should be an integer (eg, 400 / 100 = 4).
# This is because each MPI process gets allocated a fixed number of seed models.
# If it's not an integer multiple then some processes get allocated an extra seed
# (which can cause the other processes to idle while waiting).
# For example, if 400 seed models are distributed across 96 cpus then 16 cpus (= 400 % 96)
# are allocated 5 seed models (= (400 // 96) + 1) and 80 cpus are allocated only 4 seed models.
#PBS -l ncpus=96

# Total memory.
#PBS -l mem=190GB

# Time limit.
#PBS -l walltime=6:00:00

# Working directory (set to where job was submitted).
#PBS -l wd

#
# Set up the environment.
#

# Use the system MPI implementation (not conda's MPI).
module load openmpi

# Initialise the shell for conda environments to avoid the error:
#   "CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'."
source ~/.bashrc

# Activate the "optAPM" conda environment in our home Miniconda installation.
conda activate optAPM

#
# Run the job.
#
# Note: It seems "LD_PRELOAD=libmpi.so" is needed to prevent the error:
#       "[LOG_CAT_ML] component basesmuma is not available but requested in hierarchy".
#       See https://www.mail-archive.com/users@lists.open-mpi.org/msg35048.html
mpirun -x LD_PRELOAD=libmpi.so -np $PBS_NCPUS python Optimised_APM.py
