---
title: Code Features
has_children: False
parent: Tristan-MP Features and Code Structure
nav_order: 1
---

## Overview
Tristan-mp is implemented in Fortran 95, and is designed in a modular way, so that new relevant features can be added to the code with relative ease. In its present version, and without any alterations to the source code, Tristan-mp can simulate 2 different problems: **UPDATE THIS: a shock simulation (with arbitrary plasma parameters defined in the input file) and, and a generic periodic run with two electron - ion counter-streaming beams (Weibel instability).**

The code numerically solves Maxwell's equations, and the relativistic equations of motion for the particles. It is a massively parallel code, used and tested on high-performance computational clusters, and makes use of the Message Passing Interface ([MPI](https://computing.llnl.gov/tutorials/mpi/)) library, and the Hierarchical Data Format 5 ([HDF5](https://portal.hdfgroup.org/display/support)) libraries for outputs. There is an ongoing effort to vectorize the code using [OpenMP](https://www.openmp.org/). These are third party libraries that can be easily compiled for the specific system where the code is to be used (and are usually available in most of the high-performance clusters by default). On systems where there is not a specific implementation of the MPI library to be used, one possible implementation that can be downloaded and installed is the [Open MPI library](https://www.open-mpi.org/). In general, we find the none-free Intel MPI library to be slightly faster in practice.


Below you can find a general description of the code structure, and how the code is divided among the several modules. Also, you can find on this page general performance information and general rules of thumb to use when setting up a simulation. For specific information on how to setup all the keywords in the input file visit the Input Structure page. Specific instructions on how to compile the code can be found here, and specific instructions on how to run the code are available here.


## General code structure

For Tristan-mp, each logical module corresponds to a Fortran 95 module, which corresponds to a specific implementation file. To run the code, you need to compile it, producing the correct binary files (two different binary files can be produced, for the 2D and 3D versions of the code). A root directory for the run must then be created, and the input file must be copied into it.


**UPDATE THIS** The present version of the code, which is available from github (please contact the authors), consists of 20 F90 files, 1 .h file, 3 .c files, a MakefileD), and 2 sample input files. The actual number of files, however, might vary over time with different versions of the code being released. The distribution of the code among files can be roughly categorized in 7 sections: main files, main algorithms, outputs, initialization, communication, auxiliary modules, and user modules. We proceed to specify what can be found on each section.


### Main files

The main files of the code are the ones containing the calls to initialize all variables, and the main loop. These are `tristan.F90` and `tristanmainloop.F90`. These are the top level files, which control the code flow; if you are actually seeking an in-depth knowledge on how the code works, you might want to start by looking at these two files, as they represent the highest level of abstraction.

The file `tristan.F90` is the main project file and actually just calls three functions whose names are self explanatory: `initialize()`, `mainloop()`, and `finalize()`. The `mainloop()` function is the core of Tristan-mp and is defined in `tristanmainloop.F90`; this function controls the flow of the code during the main part of the simulation. Here, functions are called to update the electromagnetic fields and to push particles around the simulation space, and to deposit charge on the grid. These functions actually implement the numerical solution to the physical equations solved in Tristan-mp and are defined in the next-subsection.


### Main algorithms

The files and modules of this abstract section deal mostly with the numerical implementation of the physical equations solved by the code. For a particle code such as Tristan-mp, this involves solving the electromagnetic field equations (Maxwell's equations), and using the Electric and Magnetic fields calculated on the computational mesh to advance the particle's velocities in time via the Lorentz force equation. Additional numerical steps that must be performed include interpolating the fields from the computational mesh to the particles positions, using a first order interpolating function, and depositing the current derived from the particles' positions and velocities back to the same mesh. Also, a digital filtering function is used to suppress non-physical short-wavelength modes that result from the finite difference techniques employed.

These modules thus contain the core of the code, and are defined in the files `fields.F90` and `particles.F90`. The `file fields.F90` contains functions that are directly related to the electromagnetic field solvers, including the digital filtering functions. The file `particles.F90` includes all the functions that directly iterate over all the particle population, namely the velocity and position advance functions, and the functions used to deposit current on the mesh.


### Outputs

The outputs of the code are all in hdf5 files. All the output of the code is controlled and performed in the `output.F90` file. Here you can find the functions that create the output files and process the data to be dumped to the files.


### Initialization

Initialization procedures are done by initializing the MPI communication library, reading the input for the run, and then requesting that each module reads the variables that are relevant to it. This behavior is implemented and controlled mainly in the `initialize.F90` file, with corresponding implementations in each particular module. After reading the input variables, functions defined in each module are called in sequence for initialization. For example, the `initialize_particles()` function is called from `initialize.F90`; this function is actually defined and implemented in `particles.F90`. The purpose of this module is thus to centralize the initialization logic in one place, leaving the details of the actual module initialization to the modules themselves.


### Communication

Currently, the inter-process communication functions for each module are defined in the modules themselves (e.g. functions used to transfer parts of the computational mesh between processes are defined in `fields.F90`). This behavior however, might change in the near future, depending on how the code expands, and separate fieldcomm and particlecomm modules might be created. The purpose of the communications module file communications.F90 is to define general MPI-related functions that are useful for all the communication routines defined throughout the code. At this point this function just defines functions used to time the execution of the code, and the initialization of the random seed used in the random number generator of the code. The file `mpidummy.F90` simply re-defines the MPI functions used in the code, and will only be relevant if the code is compiled in a system where no MPI library distribution is available. The use of this feature is discouraged, as it is a simple matter of installing a free distribution of the MPI libraries to have the full functionality of the code available.


### Auxiliary modules
**OUT OF DATE, UPDATE**

The auxiliary module files are the ones handling the auxiliary features of the code, some essential, others not. These include the files `aux.F90`, `domain.F90`, `restart.F90`, `inputparser.F90`, `fparser.F90`, `system.F90`, `selectprt.F90`, and `splittestpart.F90`. These two last files define independent programs used for the particle track diagnostics of the code. In the current version of the code, the file `aux.F90` just defines the random number generator. The file `domain.F90` defines functions that allow the computational domain to grow in time, which is used for efficiency in the shock simulations. The file `inputparser.F90` links to `par.c` and `par.h` which are actually responsible for parsing the input file. The fparser.F90 module is generica mathematical function parser, that allows the user to specify arbitrary functions in the input file (which can then be used to setup the simulation). Finally, the `restart.F90` file defines functions that allow the complete state of a simulation to be stored on a disk at a given point and restarted later, and the file system.F90 defines system specific auxiliary functions.


### User module

The user file, `user_*.F90`, such as `user_shock.F90` or `user_weibel.F90`, is used for readily expanding the behavior of the code without having to alter any of the code-flow logic. The user defines routines for loading and injecting particles here, as well special boundary conditions on fields and particles (e.g. reflecting walls). The purpose of the user module is to make the expansion and generalization of the code more accessible, in particular when this does not involve modifying the underlying numerical methods used to solve the physical equations. Expansion of the code by adding other modules and modifying the core structure of the code might also be attempted, but will probably require more time, and a more complete understanding of the underlying structure of the entire code.
