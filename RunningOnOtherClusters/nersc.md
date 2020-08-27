---
title: Nersc Cori
has_children: False
parent: Running On Other HPC Clusters
nav_order: 1
---

#  Logging into Cori Nodes
First you need to get a NIM account and password with MFA. Once you have that
set up, ssh into cori.
```bash
ssh user@cori.nersc.gov
```
All of your simulations should be run from your `$SCRATCH` directory, so cd
there. *IMPORTANT* NERSC will delete files left in the $SCRATCH directory every
~2 weeks or so (check official policy if you need to)
```bash
cd $SCRATCH
```

For the most part, running `tristan-mp` on Cori is simple, except Cori no longer
supports the older version of `HDF5` that `tristan-mp` requires, so you'll need
to compile your own. Here are the instructions to compile `tristan-mp` with
`hdf5-1.8` on Cori.

Download the `hdf5-1.8.21.tar` to your `$SCRATCH` folder: located [here](https://portal.hdfgroup.org/display/ support/HDF5+1.8.21) and unzip it.
```bash
tar -xvf hdf5-1.8.21.tar
```

You'll need to compile and install `HDF5`. The easiest way to do that is create a
file called `install.sh` in the unzipped `HDF5` directory. The file should contain
the following
```bash
#!/bin/bash
H5_HOME=/PATH/TO/hdf5-1.8.21/ # modify this to your path to hdf5

./configure --prefix=$H5_HOME --enable-fortran --enable-parallel \
--enable-java --enable-shared CFLAGS="-fPIC -Ofast" \
FCFLAGS="-fPIC -Ofast" LDFLAGS="-dynamic" FC=ftn CC=cc \
CXX=CC --enable-build-mode=production --enable-unsupported

make -j 4 # parallel compile
make install # install
```
Finally you'll compile tristan following the instructions
[here](/tristan-mp-pitp/GettingStarted/Downloading-and-Compiling-Tristan),
except you will need to change your tristan-mp `Makefile` to point to version
of `HDF5` you just compiled. You can keep all parts the same except change the
following flags and variables where they occur in the original `Makefile`
```bash
cc = icc
FC = /PATH/TO/hdf5-1.8.21/bin/h5pfc # replace with the actual path
LD = /PATH/TO/hdf5-1.8.21/bin/h5pfc
PERFORMANCE = -O3 -I/PATH/TO/hdf5-1.8.21/include -L/PATH/TO/hdf5-1.8.21/lib
```


## Submit Scripts for NERSC

You should look at the official guide to running jobs on Cori, but here is a
good template of a slurm script to run on the Haswell cores.
```bash
#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=175
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell
#SBATCH -t 48:00:00                                                             
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=USER.EMAIL@gmail.com
srun -n 5600 ./tristan-mp2d -i input.a > out
```
