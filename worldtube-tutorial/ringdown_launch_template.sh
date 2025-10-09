#!/bin/bash -
#SBATCH -o spectre.out
#SBATCH -e spectre.out
#SBATCH -J worldtube_ringdown
#SBATCH --no-requeue
#SBATCH --reservation=sxs_standing
#SBATCH --nodes 2
#SBATCH --ntasks-per-node 2
#SBATCH --cpus-per-task 28
#SBATCH --constraint=cascadelake
#SBATCH -p expansion
#SBATCH -t 24:00:00
# Distributed under the MIT License.
# See LICENSE.txt for details.



# Spectre directories
export SPECTRE_BUILD_DIR=/home/nwittek/build-Release
export RUN_DIR=${PWD}
export SPECTRE_INPUT_FILE=${RUN_DIR}/Ringdown.yaml
export SPECTRE_EXECUTABLE=${SPECTRE_BUILD_DIR}/bin/EvolveCurvedScalarWaveKerrSchild3D


cd ${RUN_DIR}

# Set desired permissions for files created with this script
umask 0022

# Set the path to include the build directory's bin directory
export PATH=${SPECTRE_BUILD_DIR}/bin:$PATH

# One thread for communication
CHARM_PPN=$(expr ${SLURM_CPUS_PER_TASK} - 1)


mpirun -n ${SLURM_NTASKS} \
  ${SPECTRE_EXECUTABLE} --input-file ${SPECTRE_INPUT_FILE} \
  ++ppn ${CHARM_PPN} +setcpuaffinity +no_isomalloc_sync 

