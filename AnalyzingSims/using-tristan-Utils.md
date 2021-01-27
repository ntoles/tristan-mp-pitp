There is an automatic object relational mapping that exists for Pic Simulations here: https://pcrumley.github.io/tristanUtils/tristanSim.html 

Here is the documentation repeated here: 

TristanSim is an object relational mapping (ORM) that exposes an API to access 
*tristan-MP* simulations. All of the data attributes are lazily evaluated. 
The first time you access it the attributes are loaded from the disk, but the next time it is cached.

## Code structure & implementation

### Basic Tristan Outputs
The TristanSim object takes a path to a directory upon initialization. Then it looks in that directory for any
file matching names `flds.tot.*`, `prtl.tot.*`, `param.*` and `spect.*`. If all 4 files are present it creates
and output point for this timestep. Here's a quick example:
```python
from tristanSim import TristanSim
myRun = TristanSim('/path/to/tristan/output/')
```
The simulation object has been built. All output points are in the simulation object. For accessing it,
you can think of it as a simple list


Since you can treat the simulation object as if it just a list, we can access any of the output points using simple
list operators. e.g., `len(myRun)` gives the number of output files, access the first
one by `myRun[0]` or you can iterate over them.
```python
for out in myRun:
    # Do something for each output point here
```

To access the data on the disk you have to look at one of the objects in an output point. For instance
```python
myRun[4]
myRun[4]._flds
```
Fields is a pointer to the 5th flds.tot file in your output dir. You can see the attributes saved to the disk here by
```python
print(myRun[4]._flds.keys())
# you can get any of those attributes by typing e.g.,
myRun[4]._flds.ex
# or more simply
myRun[4].ex
```
`myRun[4].ex` is lazily evaluated. The first time it is read from the hdf5 file,
but afterwards it is in memory. If you want to delete it and reload `myRun[4].reload()`
Then you can access it from the disk by typing `myRun[4].ex`.
If you want to force all the field output to load and cache type `myRun.loadAllFields()`
or `myRun.loadAllPrtls()` to load all the prtl outputs. In practice you shouldn't have to worry too much
about the memory access if you are ok with it being lazily evaluated.

`myRun[4]._prtl`, `myRun[4]._spect` and `myRun[4]._param` are similarly defined, but if there are no collisions, any attribute should be able to be accessed as `myRun[4].c_omp` or whatever. There are a few collisions, and then you have to which one you want to default to. For instance `dens` is in `spect.*` and `flds.tot.` we default to the `flds` value. This can be set in the `__init__()` funcion of TristanSim, in the `self._colisionFixer` dictionary. If you want to add additional hdf5 output files to look for, you can do so in the `__init__()` function of the TristanSim class.

### Fancy Examples
If you have a suite of runs you had run with [automater.py](automater.md),
this class comes in handy.

Building an iterable list of all your runs in a directory:
```python
import matplotlib.pyplot as plt
import numpy as np
from tristanSim import TristanSim

# Point to the directory where the suite of runs
# In this example I have 9 different simulations
# saved in ../batchTristan
outdir = '../batchTristan'

# Let's create a list that will hold all our simulation
# instances.
runs = []

# We'll also name each run based on the directory
# it resides.
runNames = []

for elm in os.listdir(outdir):
    elm = os.path.join(outdir, elm)
    if os.path.isdir(elm):
        elm = os.path.join(elm,'output')
        if os.path.exists(elm):
            # get the directory before 'output/'
            dirName = os.path.split(elm)[0]
            dirName = os.path.split(dirName)[-1]
            runs.append(TristanSim(elm))
            runNames.append(dirName)

# Now, any run can be accessed
# e.g. runs[4], with a name runNames[4]
# to access the 3rd output time of run4: e.g. runs[4][3]

# NOTE: because runs is a 1D list of objects, treating
# it as a 2D list, e.g. runs[4,3] WILL NOT WORK.
```

Plotting ex of the 5th output timestep of each simulation in a 3 x 3 grid. Because of how I set up the automater, all output[i] are at the same physical time.

```python
fig = plt.figure()
axes = fig.subplots(3,3).flatten()
#print(axes)
j = 0
for run, name in zip(runs, runNames):
    ax = axes[j]
    istep = run[5].istep
    comp = run[5].c_omp
    ex = run[5].ex[0,:,:]
    imSize = [0, ex.shape[1], 0, ex.shape[0]]
    ax.imshow(ex,
        extent=(x*istep/comp for x in imSize),
        origin = 'lower')
    ax.set_title(name)
    #plt.colorbar()
    j += 1
    if j == len(axes):
        break
#plt.savefig('test.png')
plt.show()
```

Let's plot all the total electron energy as function of time for each run, where the colors, line-styles,
and marker styles depend on c_omp, ppc and ntimes.

```python
# First get all of the unique values of c_omp, ppc and
# ntimes from our suite of runs.

c_omp_val = list(set([r[0].c_omp for r in runs]))
ppc_val = list(set([r[0].ppc0 for r in runs]))
ntimes_val = list(set([r[0].ntimes for r in runs]))

# Lists that store what the linestyles will be.
ms = ['.', 'x', '4', '8']
ls = ['-', '--', ':', '-.']
color = ['b', 'r', 'g', 'y']

fig = plt.figure()
for run in runs:
    # In this example, we have fast moving test
    # particles that have negative indices we don't
    # want to count towards this energy.
    plt.plot([o.time for o in run],
     [np.average(o.gammae[o.inde>0]-1) for o in run],
     c = color[ppc_val.index(run[0].ppc0)],
     linestyle = ls[ntimes_val.index(run[0].ntimes)],
     marker = ms[c_omp_val.index(run[0].comp)],
     markersize = 10
    )
plt.show()
```

### Tracking Prtls.

The tristanSim class can also track particles. It requires a particular
`output.F90` of tristan, but if the files are saved properly,
you'll find all of the tracked particles in a trackedLecs and
trackedIon object. The first call builds the database which may take
awhile. It is saved afterwards.

### Examples
Let's plot a single tracked electron for these runs we didn't track the ions.

Focus just on one run for simplicity
```python
myRun = runs[0]

# plot t vs gamma for a random prtl

choice = np.random.randint(len(myRun.trackedLecs))
randPrtl = myRun.trackedLecs[choice]

# Each prtl has the following attributes: 'x', 'y', 'u',
# 'v', 'w', 'gamma', 'bx', 'by', 'bz', 'ex', 'ey', 'ez'

plt.plot(randPrtl.t, randPrtl.gamma)
plt.show()
```

This is nice, but let's say you want to find the 10 highest energy
prtls you saved. trackedLecs allows you to sort the database.

As an example, let's  sort by energy. You can pass any function here
```python
myRun.trackedLecs.sort(lambda x: np.max(x.gamma))
# now plot the botton N
for prtl in myRun.trackedLecs[:-10]:
    plt.plot(prtl.t, prtl.gamma, 'lightgray')

for prtl in myRun.trackedLecs[-10:]:
    plt.plot(prtl.t, prtl.gamma, 'k')

plt.show()
```
You can also apply a mask to your particle to choose ones
matching a certain criteria. You can pass any function
that returns a truthy value
```python
myRun.trackedLecs.mask(lambda x: np.max(x.gamma)>10.1)
plt.subplot(211)
for prtl in myRun.trackedLecs:
    plt.plot(prtl.t, prtl.gamma)
# Masks are applied successively. However you can unmask.
# Let's plot all the other prtls
myRun.trackedLecs.unmask()
myRun.trackedLecs.mask(lambda x: np.max(x.gamma)<10.1)
plt.subplot(212)
for prtl in myRun.trackedLecs:
    plt.plot(prtl.t, prtl.gamma)

plt.show()
```

All of this code can be found in `examples.py` in the main source directory.

## Support for other pic sims.

It's easy enough to tweak tristanSim to create an ORM for other simulations. I have done so for 
[PICTOR](https://rahuliitk.wixsite.com/pictor/getting-started), although unfortunately it 
uses the same keys that pictor uses. For example `sim[n].ex` in the `TristanSim` class becomes `sim[n].Ex` in the `PictorSim` class, and `sim[n].xi` becomes `sim[n].x[sim[n].flv==1]`.
