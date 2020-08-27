---
title: Running Tristan-MP
has_children: False
parent: Getting Started
nav_order: 3
---

*This section assumes you have already [logged into perseus](Logging-in-to-perseus.md) and followed the instructions in [Downloading and Compiling Tristan-mp](Downloading-and-Compiling-Tristan.md)*
## Getting ready
Now, let's run your first shock simulation!

`cd` back to your `tigress` folder
```bash
cd ~/tig
```
Generally, you'll want to run each simulation in its own directory. Let's make a directory for your first run, and `cd` into it
```bash
mkdir firstrun
cd firstrun
```

Now, we need to add our input file and executable to the `firstrun` folder.
```bash
cp ~/tig/tristan-mp-pu/exec/tristan-mp2d ./
cp ~/tig/tristan-mp-pu/user_NAME/input.myshock ./
```

Check the contents of `firstrun` folder by typing `ls`. You should be able to confirm files were copied.

## Editing the Input file
The parameters of the shock run e.g. size of the box, number of particles, etc. is set using a input file. The input file has different sections demarcated by a `<SECTION NAME>` Let's go over each section of the `input.myshock` file. Follow along by typing `nano input.shock`.
### Node configuration
```
#
#	Tristan-mp input file
#
#

<node_configuration>

sizex	= 16    # number of cpus in x direction
sizey 	= 1	# number of cpus in the y direction
```

This section defines the number of MPI jobs AKA cores/CPUs used in each direction. In general, it's best practice to completely fill nodes on your system. The total number of cores required is `sizex*sizey`. On `perseus` each node has 28 cores. Change `sizex` to 4  and `sizey` to 7. Later, when we submit the job we'll need to tell the scheduler we'll take 28 cores.

### Time
```
<time>
last 	= 600000      # last timestep
c	= .45	      # velocity of light in comp. units
                      # this defines the timestep
timespan = 60000000   # time, in seconds, available to run the problem,
                      # to not save restart before the end
```
In general you shouldn't need to change this section. `last` is the timestep where the program will terminate. `c` is the speed of light in the grid. In PIC simulations `c` is traditionally set to 0.45 which means a photon will propagate 0.45 cells in one time step. This value guarantees that the Courant condition is satisfied. **I'm not sure that `timespan` does**.  Since we're just running a short test run, let's set last to 15000.

### Grid
```
<grid>

mx0 	= 10000         # number of actual grid points in the x direction
my0 	= 128           # number of actual grid points in the y direction
mz0 	= 1 	        # ... (ignored for 2D simulations)
```

The part of the input file that sets the size of the simulation domain in cells. Since we're running the simulation for 15000 time steps, lets set mx0 to `10000`. In production runs, you should take advantage of the expanding box, but for now, let's keep things simple. Also `sizey` should be a factor of `my0`. Since we chose `sizey=7` set `my0=133`

### Load balancing
```
<dynamic_grid>

dynmx           = 1      # 1 for dynamic allocation of mx
dynmxstart      = 1      # starting timestep
dynmxinterval   = 200    # update interval
mxmin           = 7      # minimum mx/cpu

dynmy           = 0
dynmystart      = 1
dynmyinterval   = 200
mymin           = 7
```

This part of the input file are related to load balancing the code. In general, you want a roughly equal number of particles on each code, so each core is working as hard as all the others. The simulation has to wait until all cores are done with the current lap before it can start the next lap. Initially each core takes a equal sized part of the grid, but if there are strong density fluctuations, some parts of the grid become more expensive than others. If you choose, tristan-MP can try to dynamically re-allocate the grid. The parameters `dynmx` and `dynmy` are set to 1, tristan will reallocate in the `x` and `y` direction respectively. `dynm*start` Is the starting lap to reallocate the grid, this should almost always be set to 1. `dynm*interval` chooses how many laps between load balancing. `m*min` is the minimum size the dynamic decomposition will allow. If you make it less than 7 you'll spend too much time calculating ghost zones. Don't change anything here.

### Field solver
```
<algorithm>

highorder	= 0	# 0 -- 2nd order FDTD field integration; 1 -- 4th order;
                        # don't use 1 for non-relativistic flows

Corr		= 1.0	# correction for the speed of light
ntimes		= 32    # number of passes for smoothing filter (current)
```

These parts of the input change parts of the field solver. Setting `highorder` to 1 if particles are moving relativistic flows is useful because it helps with numerical Cherenkov. `Corr` is a special parameter, it is basically a way to makes E&M fields to travel slightly faster than the speed of light. For relativistic flows, you can set `Corr=1.025`. `ntimes` is the number of times to apply a 2D Gaussian smoothing filter to the current before depositing the current on the grid. It helps with high frequency fluctuations which can be caused by noise, numerical Cherenkov, or electrostatic fluctions from too low `ppc0` and/or resolution. Set `ntimes` to 4.

### Restart
```
<restart>

irestart = 0                # 1 to restart the simulation from saved
                            # restart/*.d files.
intrestart = 10000          # How often to save restart files.
                            # They overwrite previous *d files.
laprestart = 0	            # if different from 0, restart from a
                            # named restart file,
                            # saved at timestep laprestart
namedrestartint = 80000000  # interval for saving named restart files,
                            # like restart.lap01234.d
```

This part of the file is for restarting simulations after they have ended. If you want to restart a simulation, you need to resubmit the job from the same directory as before but with `irestart` equal to 1. For now, let's not do anything.

### Output
```
<output>

intervaal  = 500    # plot interval
pltstart   = 0      # starting iteration for first plot

istep      = 2	    # downsampling factor for grid output
stride	   = 20	    # particle stride for particle output

###############################################################
writetestlec    = 0         # write test electrons for tracking
dlaplec         = 90        # interval
teststartlec    = 1000000   # starting lap
testendlec      = 6000000

writetestion    = 0         # write test ions for tracking
dlapion         = 600       # interval
teststartion    = 1000000   # starting lap
testendion      = 6000000
###############################################################
```

Here we can change how often we output files. `interval` is how often the data will be saved. `ptlstart` is the first  time to start saving. When saving fields data, it is downsampled by a factor `istep` in all directions. Similarly the particle arrays are downsampled by a factor `stride`. The stuff inside of the `#` is related to particle tracking, don't worry about it for now.

### Boundaries
```
<boundaries>

periodicx	= 0	# periodic boundaries in the x direction?
periodicy	= 1	# periodic boundaries in the y direction?
periodicz	= 1     # periodic boundaries in the z direction?
```

Pretty self-explanatory. Change nothing.

### Domain
```
<domain>

enlarge		= 1   # if 1, enlarge box in the x direction if
                      # injector is close to right wall
```

Setting `enlarge` to 1 makes the box expand to the right with time. For now, let's set it to 0.

### B-field setup
```
<fields>

btheta	= 70	# bfield angle bphi=0 -> bz,
                # bphi=90 in x-y plane, btheta=0-> parallel
bphi	= 90 	#
```

We're going to run an unmagnetized shock so `btheta` and `bphi` are meaningless.

### Particles
```
<particles>

sigma	= 0.1    # magnetization number (omega_ce/omega_pe)^2,
                 # including gamma0 for inertia
maxptl0 = 1e9    # max number of particles in the simulation
ppc0 	= 16     # number of particles per cell


delgam 	= 1.e-4	 # delta gamma for ions
                 # delgam = k T_i / m_i c^2,
                 # T_i is ion temperature
me	= 1.	 # electron mass
mi	= 100. 	 # ion mass
                 # (actually mass to charge ratio)

gamma0	= 0.1	 # flow drift gamma. If < 1,
                 # interpreted as v/c.
		 # the drift of background plasma
                 # in the negative x direction
c_omp	= 10	 # electron skin depth in cells
```

Here we edit some particle & field quantities. `sigma` is the magnetization and it's definition shows up in `binit` in the user file. Let's run an unmagnetized shock; set `sigma` to zero. `maxptl0` is the total number of particles the simulation can handle, for very long runs, you may need to increase it. Keep it at 1 billion for now. `ppc0` is the total number of particles in each initialized cell. `ppc0=16` means there are 8 electrons and 8 ions. `delgam` is the temperature of the ions, and the mass of the ions is controlled via `me` and `mi`. In theory one could change both of these quantities, but you should generally just change `mi`. Let's run a pair shock, so set `mi` to 1. `gamma0` is the flow speed in the leftward direction. Set `gamma0` to 0.5, i.e., half the speed of light.

### Problem
```
<problem>

distr_dim = 3           # Dimensionality of particle distirbution;
                        # 2 means Tx=Ty=T, Tz=0; 3 means Tx=Ty=Tz=T
	                # if distr_dim = 2 and sigma is non-zero,
                        # code will reset it to 3.


temperature_ratio = 1   # T_e/T_i

external_fields = 0     # if nonzero, add external nonevolving fields
                        # to mover; defined in user file
sigma_ext = 0.          # strength of external magnetization,
                        # (omega_ce/omega_pe)^2,
                        # including gamma0

user_part_bcs=1         # call particle_bc_user routine from user file,
                        # specify particle bcs like walls
wall = 1                # left wall on?
wallgam = 0             # gamma of the moving reflecting wall on the left.
                        # If < 1, read as v/c. Set to 0 or 1 to
                        # not move the wall.
left_wall_speedup = 1.  # remove downstream a little faster;
                        # used for moving left wall, 1 is no speedup.

betainj = 0.2           # how fast injector is receding
betainj_interval = 2000 # how often to check file adjust_params where betainj
                        # can be modified as the code runs
rightclean = 1          # jumping right injector, on average moving at betainj
rightwall = 1           # reflecting particle wall at injector
```
These are the parts of the input file that are defined in the user file. I don't want to explain every parameter, their definitions can be deduced by reading the source code of your user file. Let's just focus on the ones you will change. `left_wall_speedup` periodically removes parts of the downstream of the shock. Set to zero. `beta_inj` is the average speed of the expanding box. I am not sure if it matters when `enlarge` is set to zero, but set to zero anyway. `righclean` turns on the jumping injector. The injectors moves faster than `beta_inj` on every timestep to the right but jumps backwards periodically to ensure that the average speed is `beta_inj`. Set to zero. `rightwall` puts a wall on the right side of the box. For cold flows it doesn't matter. Set to zero.

You are now done editing the parameter file. Close and save your `nano` session using `ctrl+x` and typing `y` and `enter`.

## Submitting jobs

To run a big job on a cluster, you cannot simply set it to run. Instead you send your job to a scheduler using a submit script. The scheduler then adds your job to a queue. The scheduler, not you, starts running your code. The specifics of the scheduler depends on the cluster, and you should always read the recommended submit file for the machine you are using. Perseus uses [SLURM](https://researchcomputing.princeton.edu/education/online-tutorials/getting-started/introducing-slurm).

Slurm needs a `submit` shell script. The script tells slurm how many resources you need and what job you want run. We'll create a `submit` file now. By opening a new empty text file called `submit` in nano

```bash
nano submit
```

The first line of your submit file just tells the schedule you want to interpret it using bash. Add this to your file
```bash
#!/bin/bash
```

Next, we need to tell slurm how many resources we will need. In the input file we set `sizex=4`, `sizey=7`, meaning we will start `4x7=28` MPI jobs or as slurm calls them 'tasks'. Each cpu will have its own task, the default for slurm. You request tasks with the `--ntasks=` flag (or `-n`). Add a line to your submit script so slurm knows you need 28 tasks.
```bash
#SBATCH --ntasks=28
```

It is nice to also specify the number of nodes, although technically not always 100% necessary. Each node has 28 cores, and you most likely will want to fully fill nodes. You request nodes with the `--nodes=` flag (or `-N`). Add a line to make slurm ask for 1 node
``` bash
#SBATCH --nodes=1
```

You should add some lines that makes slurm email you when your job starts and stops. Add the following lines, and change it to your personal email. Slurm may start are stop your job at unexpected times. If it is important to know when your job starts, you may want to choose an email you that will send you notifications.
```bash
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=YOUR_EMAIL@THEINTER.NET
```

Finally we need to tell slurm how long you want your job to run for using the `--time` flag (or `-t`). If your job ends before the time is up or crashes, slurm deallocates the resources. Let's run for 1 hour. Acceptable time formats include 'minutes', 'minutes:seconds', 'hours:minutes:seconds', 'days-hours', 'days-hours:minutes' and 'days-hours:minutes:seconds'. The maximum time you can run a job on perseus is 2 days.

```bash
#SBATCH --time=01:00:00
```
Now we have given slurm all of our options, and we are ready for it to run the code. In general, you'll need to make sure that your code will run after you have loaded the modules that you compiled tristan-mp with to run (i.e., intel, intel-mpi and hdf5). Add the following lines to your submit script.

```bash
module load intel
module load intel-mpi
module load hdf5/intel-16.0/intel-mpi/1.8.16
```

Now we are ready to run tristan in our `submit` script. Add the line
```bash
srun ./tristan-mp2d -i input.myshock > out
```

And submit the job by running
```bash
sbatch submit
```
Slurm will respond with `Submitted batch job <XXXXXXX>`, where
**IGNORE REMAINING PART**

## Checking in on running jobs

nano -w submit (enables you to modify the submit file)

Make sure that `-i input.shock` is included in the last line of the `submit` file, right before the ` > out`.


After completing the parameter settings, type:
```bash
sbatch submit
```
Type `showq |grep netID` to view the progress.


If you want to cancel a job, find the `slurm - xxxxx` file in the run directory and type:
```bash
scancel xxxxx
```
