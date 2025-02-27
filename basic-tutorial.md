# Worldtube tutorial

## Installing the executable

First of all, you will need to compile the worldtube executable in SpECTRE.
The source code can be found [here](https://github.com/sxs-collaboration/spectre). 

Once you have cloned the codebase into your home directory, you'll need to install it. A lot of information on installing SpECTRE can be found [here](https://spectre-code.org/installation.html).

On urania, I recommend you simply use Memo's build. Load the following list of modules:

```
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
module load ninja

source /u/guilara/repos/spack/share/spack/setup-env.sh
spack env activate env3_spectre_impi
export CHARM_ROOT=/u/guilara/charm_impi_3/mpi-linux-x86_64-smp
export PATH=$PATH:/u/guilara/charm_impi_3/mpi-linux-x86_64-smp/bin
export SPECTRE_HOME=/u/guilara/repos/spectre
source $SPECTRE_HOME/env/bin/activate
```

Next, create a build directory called e.g. worldtube-build. From it, you can now run CMake to configure your installation.A sensible example command would be (make sure to replace PATH_TO_YOUR_SPECTRE_DIRECTORY):

```
cmake -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ -D CMAKE_Fortran_COMPILER=gfortran -D CHARM_ROOT=/u/guilara/charm_impi_3/mpi-linux-x86_64-smp -D CMAKE_BUILD_TYPE=Release -D DEBUG_SYMBOLS=ON -D BUILD_SHARED_LIBS=ON -D MEMORY_ALLOCATOR=SYSTEM -D BUILD_PYTHON_BINDINGS=ON -D Python_EXECUTABLE=/u/guilara/repos/spectre/env/bin/python -D Catch2_DIR=/u/guilara/repos/Catch2/install_dir/lib64/cmake/Catch2 -D MPI_C_COMPILER=/mpcdf/soft/SLE_15/packages/skylake/impi/gcc_11-11.2.0/2021.7.1/bin/mpigcc -D MPI_CXX_COMPILER=/mpcdf/soft/SLE_15/packages/skylake/impi/gcc_11-11.2.0/2021.7.1/bin/mpig++ -D MPI_Fortran_COMPILER=/mpcdf/soft/SLE_15/packages/skylake/impi/gcc_11-11.2.0/2021.7.1/bin/mpigfortran PATH_TO_YOUR_SPECTRE_DIRECTORY -G Ninja
```

To build the worldtube project, simply run:

```
ninja EvolveWorldtubeCurvedScalarWaveKerrSchild3D -j 32
```

## The launch script

First of all, you will need a launch script. I have included a sample launch script that should be usable called `launch_template.sh`. 
You will have to adjust these lines for it work:

```
export SPECTRE_HOME=/path/to/my/spectre/home
export SPECTRE_BUILD_DIR=/path/to/my/build/directory
```

The number of nodes and timeout is controlled right at the beginning of the script:

```
#SBATCH --nodes 2
#SBATCH -t 4:00:00
```

Note that this currently uses the debug queue on urania and can therefore only use 2 nodes for a maximum of 2 hours. To use the regular queue, simply delete the line `#SBATCH -p p.debug`.

## The input template

All SpECTRE executables use an input file that determines the runtime configuration of the simulation. I have included a template input file called `input_template.yaml`. You will note that this currently has a lot of values in double angular brackets, e.g. `{{particle_charge}}`. This is because the file is a jinja template which is loaded by a python script that will fill in these unknowns.

I will briefly summarize some of the most important options here:

- `AnalyticData`: This determines the initial conditions of the simulation. The executable uses the `ZerothOrderPuncture` class which provides a decent approximation for initial conditions based on the geodesic orbit of the scalar charge. However, it is also possible to simply use zero initial data by setting `ParticleCharge` to 0.

- `BackgroundSpacetime`: The background used to solved the massless Klein-Gordon equation, i.e. the scalar wave equation. You should just be able to change the spin and mass. However, this will break a few things in the code that assume a central Schwarzschild black hole of mass 1 at rest. See [Notes on Kerr](#notes-on-kerr).

- `InitialRefinement`: Gives the initial h-refinements of the BinaryCompactObject domain.
 
- `IntialGridPoints`: Gives the initial p-refinements of the BinaryCompactObject domain. These are computed based on the `P` set by the python script. In general P of 0 or even -1 is a good place to start. Note that if the very last value (the number of grid points in the radial direction in the outermost shell) is set to a value larger than 12, the outer boundary can become unstable. If this should be necessary, simply use more h refinement there, though it is already fairly high

- `ObjectA/B`: Settings for the excision spheres. The inner radius corresponds to the initial excision radius. The excision spheres are surrounded by spheres with radius `OuterRadius` which are in turn surrounded by cubes with side length equal to the initial orbital radius. The outer radius should therefore be chosen to be roughly in the middle of the inner radius and half of the initial orbital radius.

- `Envelope` and `OuterShell`: These are the spheres surrounding the black hole and scalar charge. The radii can be adjusted according to the required orbit.

- `OrbitRadius`: These trigger when the orbit of the particle crosses a certain orbit radius. Unfortunately, I had no luck experimenting with automatically adjusting the time step so I use this trigger to adjust the time step manually. As the amount the grid is squeezed depends only on the orbital radius, I re-adjust the time step using this trigger. There are then groups of radii which use a certain CFL safety factor as I found that I need to adjust it for different orbital radii. In general, it is a good idea to run one orbit with loose CFL safety factors. If the simulation blows up because the time step is too large, lower the safety factor of the corresponding radii near where it blew up.

- `InterpolationTargets`: Here we set up the observation of the waveform. The `Spheres` target observes the scalar field Psi on different extraction spheres and writes them to file. In postprocessing, we can then compute the different spherical harmonics. Th `LMax` will give you the accuracy of the spherical harmonic projection, corresponding to the amount of collocation points that will be written to file on each extraction sphere. Even if you only want to observe the first 2 or 3 modes you will need to give a higher `LMax` to accurately resolve the lower modes. The `Radii` gives the different radii of the extraction spheres. I have found that giving radii close to the outer boundary can lead to data not being written to file for longer periods during the simulation so this should probably be avoided.

## The launch file
The launch file that will actually run a job is given by `launch_runs.py`. It creates a new directory and copies in the filled out templates of the submit script and the input file. Most important is the `config_dict` in the function `submit_job`. It contains all the parameters in the input file that were not set but will be set by jinja. I added some small comments to what each of them mean.

At the moment it is set up to run the eccentric orbit given in Figure 2 of our recent [letter](https://arxiv.org/pdf/2410.22290).

### Orbital parameters
The orbit is determined by the initial position and velocity of the particle. The configuration script is set up in a way so that the particle always starts out at the apoapsis on the x-axis with a positive velocity in the y direction. The simplest initial orbit would be circular. Simply set some initial orbital radius and set the particle velocity to one over the square root of that particle for a circular orbit. Make sure that the `Domain` option such as `OuterRadius` are sensible for the orbital radius you have chosen. 

For more complicated orbits, I usually set the apoapsis and the semi-latus rectum. The python script then attempts a root find of the Schwarzschild geodesic to get the velocity in the y direction that this corresponds to.

### Excision sphere radii
The radii of the excision spheres are changed dynamically according to equation 3 of our recent [letter](https://arxiv.org/pdf/2410.22290), corresponding to a smoothly broken power law where the excision spheres do not grow indefinitely but approach an asymptotic value. The transition radius is given by the parameter r0 (sometimes called rb) which is not very intuitive. The script therefore determines the value of r0 based on the desired worldtube radius at the ISCO. This is also done using a root search.

## Restarting a simulation
On urania, simulations can only go for at most 24 hours. The executable will therefore write a checkpoint just before that time and cancel the job. To restart a simulation, it is easiest to use the `restart_job.py` script and simply running `python restart_job.py /path/to/simulation`. 


## Analyzing a simulation
The simulations are in general analyzed with python. First, the python bindings have to be compiled using `ninja all-pybindings`. The libraries created then have to be added to your `PYTHONPATH` so python knows where to look for them. The path is given by `your_build_directory/bin/python`. 

I have included a script that contains a bunch of functions that may or may not be useful for the analysis of the simulations. Most important are the functions `extract_sim_data` which returns a large dictionary of the data collected about the particle and the self force as well as `extract_sphere_data` which computes the spherical harmonic modes. These are given in Spherepack format and the expected modes can be computed using the function `ylm`.


## Speed up hack
It is possible to accelerate the worldtube simulations by making sure that the elements neighboring the worldtube get more weight when they are distributed to the cores of the cluster. Unfortunately, this is rather hacky and I cannot include this on the main spectre branch. I have therefore included a patch `speedup.patch` that can be applied before compiling and should speed up the simulation by a factor of ~2.


## Understanding the executable
Getting the scalar charge orbit to work on Kerr orbits involves actually changing (and understanding) the executable.

A good place to understand how the worldtube executable functions is to consider the actions that the worldtube and the elements are doing each time step. This is given by a compile time list called `step_actions`. For the worldtube these are listed in `SingletonChare.hpp` and for the elements they are listed in `EvolveWorldtubeCurvedScalarWave.hpp`. The algorithm they follow is described in this [paper](https://arxiv.org/abs/2403.08864). It is probably a good idea to read at least the comments of all of the `ElementActions` and `SingletonActions` in the `src/Evolution/Systems/CurvedScalarWave/Worldtube/` folder. 

# Notes on Kerr

## The puncture field
Moving from Schwarzschild to Kerr is mostly a matter of adjusting the puncture field. On the C++ side, this is computed in the `PunctureField.hpp` file. The puncture field is a function of the particle's coordinate position and time derivatives of it and is evaluated on the worldtube boundary.

It is subtracted from the numerical field to obtain the regular field which is then send to the worldtube. This happens in `SendToWorldtube.hpp`. It then gets iterated to compute the acceleration in `IteratePunctureField.hpp`. Finally, it gets added on back to the regular field to create boundary conditions for the domain in `ReceiveWorldtubeData.hpp`.

If the self force is turned on, the particle is no longer moving on a geodesic and the puncture field requires the so-called acceleration terms which are added to the puncture field, see `IteratePunctureField.hpp`. These terms need more functions that are computed such as the self force data which is calculated by the worldtube in `IterateAccelerationTerms.hpp` and sent to the neighboring elements which compute the updated puncture field.

### Generating a new puncture field


## Derivatives of Kerr Metric
The additional acceleration terms require computation of the metric derivatives. These are used to compute the acceleration terms in `IterateAccelerationTerms.hpp` The derivatives are computed in `KerrSchildDerivatives.cpp` and need to be adjusted to include the spin of the central black hole. Take a look at `Test_KerrSchildDerivatives.cpp` for a unit test that checks the analytic derivatives against a numerical derivative.


## Generic Orbits
At the moment, the orbit is fixed in the equatorial plane. The worldtube excision sphere is via through a series of time dependent coordinate maps. In `UpdateFunctionsOfTime`, the functions of time controlling these maps are set according to the integrated worldtube position and velocity each time step. In `Test_UpdateFunctionsOfTime`, it is checked that providing the position and velocity of the particle to this action correctly moves the excision sphere. The `UpdateFunctionsOfTime` therefore need to be generalized to work for any position and velocity outside the equatorial plane. I think this can just be done by changing the `Rotation` function of time. The test can be readily adjusted to ensure that this works as intended.

