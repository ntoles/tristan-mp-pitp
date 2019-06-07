## Below is an example script to do some rudimentary analysis of tristan-MP outputs. 

The first part of the script we will load [numpy](https://www.numpy.org), the default python library for numerical work, [matplotlib](https://matplotlib.com), the default python 2D plotting library, and [h5py](https://www.h5py.org/), the library that allows us to access our data output files. We'll also load the [os](https://docs.python.org/3/library/os.html) library which will let us interact with our directories etc.

```python
import numpy as np
import matplotlib.pyplot as plt
import h5py, os
```

First let's have a variable for our output directory
```python
outdir = '/PATH/TO/YOUR/output/'
```

Tristan saves 4 files every output, `param.***`, `flds.tot.***`, `prtl.tot.***` and `spect.***`, for each output timestep; e.g. param.001, param.002, ... It's useful to look at our output directory and make sure all the timesteps that has all four of these files in our outdir

First we will make a list of all of the files in the output directory that start with 'param', 'flds', 'prtl' or 'spect'
```python
fnum= filter(lambda x: x.split('.')[0] in ['param', 'flds', 'prtl', 'spect'], os.listdir(outdir))
```
fnum is an iterable that has the values= ['param.001', 'flds.001',... ,'spect.NNN'] We only need the last bit of the file that is the output number
```python
fnum = list(map(lambda x: x.split('.')[-1], fnum))
```
Fnum is a list that looks like ['001', '001', ... , 'NNN', 'NNN']. We will want the outputs that have all 4 output files. Convert fnum into a set so that only unique values are saved. Then, only keep values that have 4 elements in fnum

```python
fnum = sorted(filter(lambda x: fnum.count(x) == 4, set(fnum)))
```
fnum is a list (technically an iterable) that contains all of the '***' endings to the outfiles, but only 1 of each (['001', '002', ...]). Soon we'll iterate over all the files to plot Bsq vs time, but beforehand let's see what values each outfiles have saved.

We can see the keys available to us by loading the file with h5py

```python
for elm in ['param.', 'flds.tot.', 'prtl.tot.', 'spect.']:
    with h5py.File(outdir+elm+fnum[0], 'r') as f:
        print('The quantities available in the ' + elm + 'files are...')
        print([key for key in f.keys()])
```

See tristan-mp output.F90 file, or ask a friend to figure out what each of these keys are.

Let's say we want to plot the total B-field as a function  of time. First make a list to hold the Btot values and time values.
```python
Btot = []
t = []
```
Iterate over the output numbers.
```python
for num in fnum:
    # For each output, we have to load the param file to get the time
    with h5py.File(outdir+'param.'+num, 'r') as f:
        # Append it to our list
		t.append(f['time'][0])

    # Now, load the fields file and get the B value
    with h5py.File(outdir+'flds.tot.'+num, 'r') as f:
        # Let's assume the simulation is 2D
		# The simulation saves Bx, By, Bz separately
		# so we must square them, take square root, then
		# sum over the box.
        Btot.append(np.sum(np.sqrt(f['bx'][0,:,:]**2+f['by'][0,:,:]**2+f['bz'][0,:,:]**2)))
```
We can plot our list using matplotlib

```python
plt.plot(t, Btot)
# Add labels to the plot
plt.xlabel(r'$\omega_{pe}t$')
plt.ylabel(r'$B^2 \ [{\rm arb. unit}]$')
plt.title('Total Magnetic energy vs time')
plt.show()
```

Here's the script is in its entirety:

```python
import numpy as np
import matplotlib.pyplot as plt
import h5py, os

#First let's have a variable for our output directory

outdir = '/PATH/TO/YOUR/output/'


#Tristan saves 4 files every output, `param.***`, `flds.tot.***`, `prtl.tot.***` and `spect.***`, for each output timestep; e.g. param.001 param.002, ... It's useful to look at our output directory and make sure all the timesteps that has all four of these files in our outdir

#First we will make a list of all of the files in the output directory that start with 'param', 'flds', 'prtl' or 'spect'

fnum= filter(lambda x: x.split('.')[0] in ['param', 'flds', 'prtl', 'spect'], os.listdir(outdir))

#We only need the last bit of the file that is the output number

fnum = list(map(lambda x: x.split('.')[-1], fnum))

#We will want the outputs that have all 4 output files. Convert fnum into a set so that only unique values are saved. Then, only keep values that have 4 elements in fnum

fnum = sorted(filter(lambda x: fnum.count(x) == 4, set(fnum)))

#Now fnum is a list (technically an iterable) that contains all of the '***' endings to the outfiles. Soon we'll iterate over all the files to plot Bsq vs time, but beforehand let's see what values each outfiles have saved.

#We can see the keys available to us by loading the file with h5py


for elm in ['param.', 'flds.tot.', 'prtl.tot.', 'spect.']:
    with h5py.File(outdir+elm+fnum[0], 'r') as f:
        print('The quantities available in the ' + elm + 'files are...')
        print([key for key in f.keys()])

#See tristan-mp output.F90 file, or ask a friend to figure out what each of these keys are.

#Let's say we want to plot the total B-field as a function  of time. First make a list to hold the Btot values and time values.

Btot = []
t = []
 
# Iterate over the output numbers.

for num in fnum:
    # For each output, we have to load the param file to get the time
    with h5py.File(outdir+'param.'+num, 'r') as f:
        # Append it to our list
        t.append(f['time'][0])

    # Now, load the fields file and get the B value
    with h5py.File(outdir+'flds.tot.'+num, 'r') as f:
        # Let's assume the simulation is 2D
        # The simulation saves Bx, By, Bz separately
        # so we must square them, take square root, then
        # sum over the box.
        Btot.append(np.sum(np.sqrt(f['bx'][0,:,:]**2+f['by'][0,:,:]**2+f['bz'][0,:,:]**2)))

#Now we can plot our list using matplotlib


plt.plot(t, Btot)
# Now we should add labels to the plot
plt.xlabel(r'$\omega_{pe}t$')
plt.ylabel(r'$B^2 \ [{\rm arb. unit}]$')
plt.title('Total Magnetic energy vs time')
plt.show()
```