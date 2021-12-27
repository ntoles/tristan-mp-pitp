# tristan-mp-pitp
TRISTAN-MP parallel electromagnetic 3D particle-in-cell code.
Pre-release version for testing. 

Developed by: Anatoly Spitkovsky, Luis Gargate, Jaehong Park, Lorenzo Sironi. 
Based on original TRISTAN code by Oscar Buneman. 

Refer to the [wiki](http://ntoles.github.io/tristan-mp-pitp) for more extensive documentation. 
Warning: some of the info on the wiki is obsolete. 

Organization and suggested workflow: 
The code is modular and is organized in such a way that in most cases 
the user would not need to edit any of the main source files in the main 
tristan-mp directory. 

All of the user-specific configuration should be possible to confine
to user_* routines. There are three example user configuration files 
in main directory, showing a counterstreaming Weibel instability
setup, two-stream instability and a collisionless shock simulation. 
There are also sample input files.

When using code from github, clone it to your local machine, 
switch to "master" branch. You can create your branch off of master.

On local machine:
git clone https://github.com/ntoles/tristan-mp-pitp.git

copy some example files to start with

cp user_weibel.F90 user_mysetup.F90

cp ../Makefile Makefile.mysetup

edit user_mysetup.F90
edit Makefile.mysetup to add USER_FILE=user_mysetup 
(no need for extension F90)

To compile
cd source directory 
make -f Makefile.mysetup clean
make -f Makefile.mysetup

3D version is enabled when -DtwoD flag is omitted from the Makefile. 
 
You need to have parallel HDF5 library installed with intel or GNU compilers, 
which will create h5pfc alias for the compiler. Some instructions for 
installation are on wiki page. For Macs brew seems to work fine with gfortran:

$ brew install gcc
$ brew install openmpi --enable-mpi-thread-multiple
$ brew install hdf5 --with-fortran --with-mpi

This will produce tristan-mp2d executable. 

To run:
Make a run directory somewhere outside the git-controlled source directory. 
Copy the executable tristan-mp2d there. 
Copy example submit and input files from directory 
(see wiki page for example submit
files for clusters; you don't need these on your desktop/laptop).
 
Input file has to be named "input" in the run directory, or the executable takes -i option. 
E.g.: 
./tristan-mp2d -i input.weibel
(for MPI, it can be, e.g.: srun -n 16 ./tristan-mp2d -i input.weibel)
Note that you need to edit the input file to set the number of domain sub-partitions
 sizex * sizey * sizez be equal to the total number of cores to be used. sizez = 1 in 2D. 

Edit submit and input files for your case and according to the queue policy of your cluster. 
In 3D, the domain is split in y and z directions, 
with "sizey" chunks in y direction and "total # of cpus/sizey" chunks in z direction. 

>qsub submit 
or other appropriate submission command. 

When running, the code will create subdirectories output and restart.
The output is single HDF5 files (single per time step). 
Currently we provide routines to interactively view output using python.
https://github.com/pcrumley/Iseult

It requires anaconda to run. We had good experience with anaconda 4.0.0 
on Mac, but not later. The older version is available on anaconda's website. 

To launch the vis tool, run (path to Iseult)Iseult/iseult.py .
This will open interactive windows. Right click on plots to get more options. 

load.py is a script that loads HDF5 files into a python dictionary, 
that can be accessed as d[i]['bz'], where i is the file number. 

There are also older IDL routines, which are available on request. 

Good luck! 

