# FRTPhase

FRTPhase is a small toolset that was developed in connection
with my master thesis about an analytic approach to the **F**irst 
**R**eturn **T**ime **Phase**. Trajectories of
isotropic planar stochastic oscillators, that are defined
by a set of stochastic differential equations (SDEs), are
simulated. In doing so, the times of first first passage through an arbitrary
surface in state space are detected. Thus, this toolset enables 
one to verify such surfaces as isochrons of the 
FRT Phase. The toolset is
divided into two parts: a Python toolbox and a C++ project.

The Python toolbox "PyTools" consists of a package "Toolbox"
and the applications that implement the functionality provided
by the package. 

The C++ project consists of three libraries and applications that
implement the simulations or rather in-silico "Experiments".

In "ExampleFiles" examples of the input and output files needed
or created by the applications are provided.

## Dependencies

**PyTools**: 
* Python 3
* NumPy
* SciPy
* [h5py](https://www.h5py.org/)

**FRTPhase**:
* [cmake](https://cmake.org/): building the repository
* [Boost](https://www.boost.org/): unit testing, 
  json file parsing, logging
    
* [HDF5](http://h5cpp.org/): binary file IO, 
  C++ **and** C API required
  
* [Open MPI](https://www.open-mpi.org/software/ompi/v4.1/):
parallelization on a local and cluster level

## Installation and Build

Clone the repository by typing

`git clone https://github.com/holzhauk/FRTPhase.git`.

Then, create a build directory and build the project 
using cmake

`mkdir FRTPhase-build/ && cd FRTPhase-build`

`cmake ../FRTPhase/ && make`.

## Usage pipeline

You can use the libraries supplied with this repository to
create your own applications or you can use the 
"Experiments" out-of-the-box. The following lines explain the 
intended pipeline using the included applications and example files.

Model and simulation specifications are defined in .json
configuration files (see ExampleFiles/configs). 

Python applications are used to create the binary .h5 
files containing representations of the isosurfaces.

The isosurfaces can be verified as such with 
*IsoSurfaceVeridication*. This programm was implemented 
using the Message Passing Interface for parallelization.
This example shows local execution using 3 threads (cpus).

`FRTPhase/Experiments/IsoSurfaceVerification$ 
mpirun -np 3 --use-hwthread-cpus 
./IsoSurfaceVerification 
../../../ExampleFiles/configs/config_verification.json`

Serial correlation coefficients between the time 
intervals are measured with *IsoSurfaceCorrelation*

`FRTPhase/Experiments/IsoSurfaceCorrelation& 
./SerialCorrelation 
../../../ExampleFiles/configs/config_correlation.json`

The results are written to .h5 files containing first passage time 
values, statistics and information about how this data was 
generated (path to config file and isosurface file, 
simulation parameters).

## Documentation

A comprehensive documentation is in progress and 
is going to be released after the paper associated 
with this project has been published.

## Authors

Konstantin Holzhausen




