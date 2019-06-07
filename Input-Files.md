**NEEDS UPDATING: VERY OUT OF DATE**

When `tristan-mp` is run, it looks for a file named `input` inside of the current directory. The `input` file contains all of the 

## Input file structure

This page describes in detail what each keyword in the input file does, the possible values and the default values assumed by the code.

Below there is a sample input file for reference purposes. This input specifies a 2D simulation with 800 cells in the x direction, and 256 cells in the y direction. The simulation domain is decomposed in 16 blocks (CPU's / cores), along the y direction. The maximum number of particles supported is 6 x 10^7, and the physical problem being run is a Weibel instability simulation in external magnetic field. Strength of external field is set by
parameter sigma_ext. The velocity of light is 0.45 (simulation units), and the magnetization of the flow (sigma) is 0. The electron collisionless skin depth ( c / omega_pe) is 8, and the bulk flow velocity of the particles has a gamma of 100. Also, there are 10 simulation particles per cell in the box initially, the electron mass is set to 1, and the ion mass is set to 1. For a detailed description of each keyword in the input file check the Detailed Keyword Description section.

The input file is structured around named sections, which should stand on its own line, have a unique name, and be surrounded by <> (as in <node_configuration>). Any line starting with a # is a comment and is ignored by the input parser (any information after the # in a given line is ignored). Furthermore, any variables which are specified but not explicitly read in the code are just ignored, and if some keyword is omitted, a default value is usually assumed (which might lead to unexpected results). There are currently 10 sections that are currently being parsed: <node_configuration>, <time>, <grid>, <algorithm>, <restart>, <output>, <boundaries>, <domain>, <fields>, <particles>, and <problem>. Each section is concerned with a specific part of the algorithm, or with the control of some specific behavior of the code.

## Sample Input File

```
#
#

<node_configuration>

sizey = 16# number of cpus in the y direction

<time>

last = 2000# last timestep

c= .45# velocity of light in comp. units (this defines the timestep)

timespan= 86400# time, in seconds, available to run the problem

<grid>

mx0 = 800 # number of actual grid points in the x direction
my0 = 256# number of actual grid points in the y direction
mz0 = 256 # ... (ignored for 2D simulations)

<algorithm>

conserv = 1 # charge-conservative current deposition -- the only available option
highorder= 0 # 0 -- 2nd order FDTD field integrateion; 1 -- 4th order;
# don't use 1 for non-relativistic flows

Corr= 1.025 # correction for the speed of light

ntimes= 4# number of passes for smoothing filter (current)
cleanfld= 0# number of passes for smoothing filter (fields). don't use.
cleanint= 10# interval for field cleaning; don't use.

cooling= 0# cool particles? ; not implemented
acool= 10.# cooling parameter for particles
splitparts = 0# split particles to improve statistics?

<restart>

irestart= 0# 1 to restart the simulation from saved restart/*.d files.
intrestart= 1000# how often to save restart files. They overwrite previous *d files.
laprestart= 0# if different from 0, restart from a named restart file, saved at timestep laprestart
namedrestartint = 1000000 # interval for saving named restart files, like restart.lap01234.d

<output>

interval = 20# plot interval
torqint= 2000000# interval for outputs at different resolution (currently broken)
pltstart= 0# starting iteration for first plot

istep= 1# downsampling factor for grid output
istep1= 4# downsampling factor for grid output every torqint steps
stride= 10# particle stride for particle output

writetestpart= 0# write test particles?
selectprt= 0# re-trace the same selected particles?

<boundaries>

periodicx= 1# periodic boundaries in the x direction? Choose 0 for radiative boundaries.
periodicy= 1# periodic boundaries in the y direction?
periodicz= 1# periodic boundaries in the z direction?

<domain>

enlarge= 0# if 1, enlarge box in the x direction if injector is close to right wall?
movwin = 0# if 1, use moving window
shiftinterval = 20# how often to apply moving window (in steps)
shiftstart = 1000# at what step to start shifting moving window
movwingam = 5. # gamma factor of moving window. If > 10000, it moves at c.
# if < 1, it is interpreted as v/c.
<fields>

btheta= 85# bfield angle bphi=0 -> bz, bph=90 in x-y plane, bth=0-> parallel
bphi= 90#

<particles>

sigma= 0.# magnetization number (omega_c/omega_p)^2, including gamma0
maxptl0 = 6e7# max number of particles in the simulation
ppc0 = 10# number of particles per cell

delgam = 1.e-4# delta gamma (temperature control)
me= 1.# electron mass
mi= 1.# ion mass (actually mass to charge ratio)

gamma0= 100.# flow drift gamma. If < 1, interpreted as v/c.

c_omp= 8# electron skin depth in cells

<problem>
caseinit = 1 #can be used to select subcases of the problem. Not used.

#density_profile=0. #x,y,z dependent function distributing initial particle weight
temperature_ratio = 1 # T_e/T_i

external_fields = 1 # if nonzero, add external nonevolving fields to mover
sigma_ext = 0.1 # strength of external magnetization,(omega_c/omega_p)^2,
# including gamma0


user_part_bcs=1 # call particle_bc_user routine from user file, specify particle bcs like walls

wall = 0 # use reflecting wall? Position is set in user file.
wallgam = 0. # gamma of the moving reflecting wall. If < 1, read as v/c. Set to 0 or 1 to not move the wall.
```

Keyword Specification

    mx0, my0, and mz0

Actual number of grid points in x,y and z directions

    MAXPTL0

<detailed description goes here>

.... 