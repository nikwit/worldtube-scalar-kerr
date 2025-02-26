#!/bin/bash -
#SBATCH -J {{sim_name}}
#SBATCH --nodes {{num_nodes}}
#SBATCH -t 4:00:00
#SBATCH -p p.debug
#SBATCH -o spectre.out
#SBATCH -e spectre.out
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=36
#SBATCH --no-requeue
# Distributed under the MIT License.
# See LICENSE.txt for details.

# To run a job on Minerva:
#
# 1. Set the -J, --nodes, and -t options above, which correspond to job name,
#    number of nodes, and wall time limit in HH:MM:SS, respectively.
# 2. Set the build directory, run directory, executable name,
#    and input file below.
#    Note: The executable will not be copied from the build directory, so if you
#    update your build directory this file will use the updated executable.
# 3. Optionally, if you need more control over how SpECTRE is launched on
#    Minerva you can edit the launch command at the end of this file directly.
# 4. Submit the job to the queue:
#    ```sh
#    sbatch Minerva.sh
#    ```
#
# Load compiler and MPI modules with explicit version specifications,
# consistently with the versions used to build the executable.
module purge
module load gcc/11
module load impi/2021.7
module load boost/1.79
module load gsl/1.16
module load cmake/3.26
module load hdf5-serial/1.12.2
module load anaconda/3/2021.11
# For plots and visualization
module load paraview/5.10
module load texlive/2021
module load gnuplot/5.4

source /u/guilara/repos/spack/share/spack/setup-env.sh
spack env activate env3_spectre_impi
export CHARM_ROOT=/u/guilara/charm_impi_3/mpi-linux-x86_64-smp
export PATH=$PATH:/u/guilara/charm_impi_3/mpi-linux-x86_64-smp/bin

# Spectre directories
export SPECTRE_HOME=/path/to/my/spectre/home
export SPECTRE_BUILD_DIR=/path/to/my/build/directory
export SPECTRE_RUN_DIR=${PWD}
# Load python environment
source $SPECTRE_HOME/env/bin/activate

echo "Spectre home directory: ${SPECTRE_HOME}"
echo "Build directory: ${SPECTRE_BUILD_DIR}"
echo "Run directory: ${SPECTRE_RUN_DIR}"

# Choose the executable and input file to run
export SPECTRE_EXECUTABLE=${SPECTRE_BUILD_DIR}/bin/EvolveWorldtubeCurvedScalarWaveKerrSchild3D
export SPECTRE_INPUT_FILE=${SPECTRE_RUN_DIR}/input_file.yaml

cd ${SPECTRE_RUN_DIR}

# Set desired permissions for files created with this script
umask 0022

# Set the path to include the build directory's bin directory
export PATH=${SPECTRE_BUILD_DIR}/bin:$PATH

# One thread for communication
CHARM_PPN=$(expr ${SLURM_CPUS_PER_TASK} - 1)
echo "Slurm tasks: ${SLURM_NTASKS}"
echo "Slurm cpus per task: ${SLURM_CPUS_PER_TASK}"
echo "Charm ppn: ${CHARM_PPN}"

module list
spack find

# Run the program:
srun ${SPECTRE_EXECUTABLE} ++ppn ${CHARM_PPN} --input-file ${SPECTRE_INPUT_FILE} +pemap 0-34,36-70 +commap 35,71 
