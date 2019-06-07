**THIS SECTION IS OLD AND DEPRECATED. PLEASE SEE [THIS](https://github.com/PrincetonUniversity/tristan-mp-pu/wiki/Downloading-and-Compiling-Tristan) INSTEAD**

## Compilation
There are three different steps involved in the compilation of Tristan-mp. The first two steps are the compilation of the MPI and HDF5 libraries if they are not already available on the system you are going to run. For most cases, there will be a preferred system-specific MPI implementation for most of the high-performance computational clusters where you should be running Tristan-mp. You should always use this MPI library implementation; just be aware that Fortran bindings to the library are necessary (they usually are included by default). Regarding the HDF5 libraries, when available they still might not have the required fortran bindings and/or parallel output bindings, which are a requirement of the code. As such you might have to compile your own version of the HDF5 library. The following two sections show you how to do this, but they are not intended as a substitute for the documentation of the libraries, and assume a basic knowledge of Unix systems. For more information read the documentation available with the MPI and HDF5 library distributions.

## MPI library
The preferred MPI distribution to use (other than any system specific distribution available), is the Open MPI library available for download [here](https://www.open-mpi.org/). After downloading the source files for the library (as a .tar.gz file or similar), you will need to copy this to the remote system where it is to be installed. To do that you will usually use a command such as scp, as in `scp downloaded-file username@remotesystem:` . You might be asked for your password to the remote system, and this should copy the file to your home directory there. Next step is to login into the system (usually via ssh).

You will now have to untar / unzip the .tar.gz file, using a command such as tar -xzf filename or similar, and this will create the installation directory. Change directories to the installation directory.

Now you will have to figure out the appropriate configuration options for the library. You will do this by typing `./configure --help`. At the very least, you will have to change the final installation directory, and make sure that Fortran bindings are included in the distribution (they usually are by default). Also, you will have to configure the correct C and Fortran compilers for the library to use (such as gcc and gfortran, or icc and ifort). These should be the compilers that are recommended for the system you are using.

After figuring out the configuration options, you will just have to run: ./configure (passing in the options), make and make install. This will produce the final library installation diretory with the .dylib (and / or .a files), and the executable files mpicc and mpif90, among others. These two will be necessary for the next step.

## HDF5 library
Installing the HDF5 libraries is similar to the process described above for the MPI libraries. You first download the source files for the library [here](https://portal.hdfgroup.org/display/support). After that, you copy to the remote system, and untar / unzip as described above.

Once logged in the system, the first step is to figure out if the mpicc and mpif90 compiler wrappers are available in the system. When you type `which mpicc` and `which mpif90`, you should see the complete path to these executables; the path should correspond to the appropriate MPI library that you will be using (either a system one, or the one you installed). For this to work you might have to add the path for mpicc and mpif90 to your path variable, if you did not do so already. You do this by typing 
```
setenv PATH ${PATH}:path-to-the-mpi-distribution-bin-folder 
```
(for csh, tcsh, or similar), or 
```
export PATH=${PATH}:path-to-the-mpi-distribution-bin-folder 
``` 
(for sh, bash, ksh, or similar).

Next, you will have to determine the correct options to pass to `./configure`, as in the installation for the MPI library above. At the very least, you will configure a correct destination folder (a sub-directory on your home folder), and turn on Fortran bindings and parallel bindings. Also, the compilers to use will be mpicc and mpif90 for C and Fortran, respectively.

After this, you will run `./configure` (passing in the options), `make` and `make install`. This will produce the final library installation directory with the .dylib (and / or .a files), and the executable file h5pfc, among others. This file will be necessary for the next step.

Example:
```bash
./configure --help
```
this lists the options
```bash
export CC=mpicc 
export CXX=mpicc
export FC=mpif90
configure --prefix=/home/username/bin/hdf5 --enable-parallel --enable-fortran
make install
```

## Tristan-mp compilation
After downloading the code source, and copying the source to your cluster, you can compile the code to produce binaries for 2D or 3D version.

As a first step, you have to make sure that the h5pfc tool is available in the HDF5 library distribution you are using. This file should be in the bin directory of the library distribution. If it is not, you will have to recompile the HDF5 libraries following the instructions above. Once you know where this executable is located, add it to your path, so that it is accessible by the Tristan-mp Makefile (and Makefile3D). You usually do that by executing 
```
setenv PATH ${PATH}:path-to-the-library-bin-folder
``` 
(in c type shells - csh, tcsh and such), or 
```
export PATH=${PATH}:path-to-the-library-bin-folder 
```
(in sh, bash, ksh, or other similar shells). It is a good idea to add these commands to your `.bashrc` or `.cshrc` file.

If when you type which h5pfc you get the path to the binary file, you are set to go. Just change directory to the Tristan-mp distribution's source folder and type `make` to get the 2D version of the code).

**THIS IS IMPORTANT!**

The code is organized in a modular way, so that in most cases the user does not need to edit any of the source files except for the `user_*.F90` file. Two examples, `user_weibel.F90` and `user_shock.F90` are provided. 

Corresponding input files are also in the user directory.
The choice of which user file will be compiled is made by modifying the USER variable in the source/Makefile or Makefile3D.

The standard workflow for a new problem definition is to copy one of the example `user_*` files to a new `user_myrun.F90` file (where `myrun` is something descriptive) and editing this user file. There are routines for initializing particles, injecting particles during the run, and special particle and field boundary conditions which are listed in the user_*.F90 file. When ready, there is no reason to edit the Makefile, instead run
```bash
make clean
make USER_FILE=user_myrun
``` 
and set USER=user_myrun, or whatever is the name of your user file omitting .F90

IT IS VERY EASY TO FORGET THE USER_FILE flag, and the code will compile the old user module. Not very useful.
The name of the user file is printed in the beginning of the out file when the code runs. So, check that if the result does not make sense.

Note that you might also have to alter some of the pre-processor and/or processor options under the keyword CUSTOM in the Makefile (Makefile and Makefile3D), which reads something like `CUSTOM= -DMPI -DHDF5 -DserIO -fpp -DtwoD`. This is just defining keywords to be passed to the pre-processor, and since the way to do this is compiler dependent you might have to change it. The code in the current version of the repository will compile by default using gcc and gfortran.

*Some notes:*
When working with both 2D and 3D version, make sure to run "make clean" or "make -f Makefile3D clean" before switching to compilation in different number of dimensions. Some of the libraries have the same names, but are compiled with different flags.

Not doing this will in general lead to a compilation failure.
