*This section assumes you have already [logged into perseus](https://github.com/PrincetonUniversity/tristan-mp-pu/wiki/Logging-in-to-perseus) and followed all the instructions to ensure you have loaded all the modules required to compile tristan-MP*

## Downloading Tristan-MP from Github
Inside of perseus, cd to your tigress directory
```bash
cd ~/tig
```
Download tristan-mp from github.com to perseus by cloning the directory.
```bash
git clone https://github.com/PrincetonUniversity/tristan-mp-pu.git
```
Go into the `tristan-mp-pu` folder
```bash
cd tristan-mp-pu
```
## Creating your own Tristan-MP branch and user files
Create your own branch of Tristan-MP where you can edit your user-files
```bash
git checkout -b YOUR_NAME
```
The structure of Tristan is such that you should not have to edit the main code to change the set-up of the files (more details [here](https://github.com/PrincetonUniversity/tristan-mp-pu/wiki/Code-Features)). Instead, the simulation set-up is changed using user_files that define important functions in the code. You should use your own user files. Make a new directory where you will save your user files    
```bash    
mkdir user_NAME
```
Copy the code you want to work with. We'll set-up Tristan-mp to run a shock problem. Each problem has its own user file and input file.
```bash
cp user_pu/user_shock_pu.F90 user_NAME/user_myshock.F90
cp user_pu/input.shock_pu user_NAME/input.myshock
```
There is a line in user_myshock.f90 that writes the userfile name in the output while tristan-mp is running. This is a very convenient feature that lets you know what you compiled your run with, but must be updated for every new user file. Find the line. The easiest way to do this with nano is `nano user_NAME/user_myshock.F90`. Inside of nano, you can search for a string with `ctrl+W` then type 'user_shock_pu'. Change the line from
```fortran
if(rank.eq.0) print *, "Using user file user_shock_pu.F90"              
```
to
```fortran
if(rank.eq.0) print *, "Using user file user_NAME/user_myshock.F90"              
```
Save and close the file in nano by pressing `ctrl+x` and then `y` and `enter` to the prompts.

## Compiling
Compilation in Tristan-MP is handled with [make](https://www.gnu.org/software/make/) and a Makefile. *Make* is a very full-featured program with many options. Luckily, you shouldn't have to mess too much with the Makefile at first.

`cd` into  the `~/tig/tristan-mp-pu` folder. We will want to compile tristan-mp, but before doing so it is good practice to first clean the directory in case you changed some files: type 
```bash
make clean
```
While not always necessary, I think it is good habit to always run `make clean` before running `make`. Now compile tristan-mp pointing it towards your user file.
```bash
make USER_FILE=user_NAME/user_myshock
```
The executable "tristan-mp2d" should appear in the newly created `~tig/tristan-mp/exec` directory, check that it is there.
```bash
ls exec
```

If instead you get a warning stating something along the lines of

```bash
h5pcc -O3  -c code/system.c -o obj/system.o
make: h5pcc: Command not found
make: *** [obj/system.o] Error 127
```

That means you have not loaded the proper modules to compile the code. Please ensure you are logged into perseus by typing ```bash hostname```. If you are logged into perseus, check that you have loaded all the neccessary modules by typing ```bash module list``` you should find the following modules ```intel, intel-mpi, hdf5/intel-16.0/intel-mpi/1.8.16``` See more [here](https://github.com/PrincetonUniversity/tristan-mp-pu/wiki/Logging-in-to-perseus). Remember to run ```make clean``` before trying to recompile.