**This page is deprecated from previous instructions on tristan-MP. Please see [this](https://github.com/PrincetonUniversity/tristan-mp-pu/wiki/running-your-first-tristan-mp-simulation) instead.**

## Running Tristan-mp

To run the code you need first to have the binary files prepared for the system in question. For instructions on how to compile, please refer to the Compilation section.

What follows are just general procedures that you will have to do in order to run the code; for more details on how to prepare an input file you should refer to the Input Structure section, or / and to the two examples provided (the Weibel instability, and the Collisionless shock run).

## Preparing a run on your local machine

After having compiled the code you should have two binary files, by default named tristan-mp2d and tristan-mp3d (if you are not interested in 3D runs, you do not have to compile the 3D version). Prepare the input file to run the problem that you want, and then create a directory and copy both the binary and the input file there. The input file must be named "input". After changing directories to the run directory, just run 
```
mpirun -n <ncpus> <name-of-executable> 
``` 
(e.g. `mpirun -n 2 ./trist-della2d` for a 2D run using two cpus).

## Preparing a run for a regular high-performance cluster

As a first step, please refer to the specific instructions for the cluster where you want to run. These instructions are not to supersede any specific site instructions, but usually running this kind of code in a high-performance cluster goes like this.

The first step to run the code is to create a root folder for the run; this folder should be accessible from all the CPU's (cores) in the run. Usually, there is a specific file system on which to run codes. You should then create a folder there using `mkdir name-of-the-folder`. The next step is to copy the pre-prepared input file for the specific run into the root folder. This input file must be named input.

Next, copy the correct version of the binary file also to the root folder you created `cp path-to-binary/tristan-mp2d path-to-run-folder/`.

Finally, create a submit file, which should be used to submit the job to the batch-system (usually PBS or MOAB), which will look something like this:

```
#!/bin/bash
#PBS -q regular
#PBS -m abe
#PBS -j oe
#PBS -l nodes=8:ppn=8
#PBS -l walltime=12:00:00


cd $PBS_O_WORKDIR

export PBS_JOBNAME="my1024"
echo $PBS_JOBNAME > jobname
echo $PBS_JOBID > jobid
cat $PBS_NODEFILE > hosts.$PBS_JOBID
mpiexec --mca btl ^tcp ./tristan-mp2d > out
```

## System dependent tasks
For completeness it should now be stated that the final step to run the code is probably to submit the run script (as the one above) to the correct queue. This is usually done in most clusters by issuing the command `qsub <submit-filename>` (for torque-type queuing systems). The specific command might vary, and the contents required to be in the run script file will most certainly vary as well. To determine what to do for each specific system you must consult the cluster user documentation / webpage.

## Code output
The code saves log messages in the file "out", sequentially numbered output files in a subdirectory "output", and restart files in a subdirectory "restart". The code overwrites present files, and creates subdirectories if they don't exist. 