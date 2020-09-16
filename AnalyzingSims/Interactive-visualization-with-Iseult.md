---
title: Interactive Visualization with Iseult
has_children: False
parent: Analyzing Simulations
nav_order: 2
---

After running a simulation, one wants to interactively examine its outputs. We
have provided a GUI written by Patrick Crumley (patrick.crumley@gmail.com),
named [iseult](https://github.com/pcrumley/Iseult)


It is capable of plotting many field and particle quantities in interactive way.
![Iseult Set-up](https://ntoles.github.io/tristan-mp-pitp/assets/IseultPanels.png)


An example visualization of a Tristan-MP simulation. It's basically a skinned
version  of the interactive matplotlib figure, where you can go through your
timesteps by  pressing the arrows left and right on the bottom. All graphs (and
graph types) are configureble by right clicking on the subplot
![Iseult Setup](https://ntoles.github.io/tristan-mp-pitp/assets/panelSettings.png)

Things like the number of rows and columns can be changed by clicking the
settings button on the bottom:
![Iseult Setup](https://ntoles.github.io/tristan-mp-pitp/assets/generalSettings.png)

If you get a set-up you like and want to save for the future, go to the file
menu and click 'Save Current State.' This saves a config file in the
.iseult_configs file. The config files are directly editable (you have access to
change things that you can't in the GUI), and they can be loaded in future using
'Preset Views' menu.

![Iseult Setup](https://ntoles.github.io/tristan-mp-pitp/assets/fileOpts.png)

Also in the file menu is the option to save movies using ffmpeg. If ffmpeg is
not  installed, you can make a movie by using the 'record' checkbox on the
bottom and click play. This creates a 'Movie' folder and saves pngs everytime
Iseult reloads a plot.  You can then put these together into a
movie/gif/whatever.

Written by:

Patrick Crumley, patrick.crumley@gmail.com, based on Jaehong's Tristan analysis
IDL script.

Python packages required to run Iseult include: Anaconda 2019/3, matplotlib 3.0
& its required dependencies, python 3.7.3, h5py. Will not work with older
versions of anaconda3 and matplotlib 2.0 or older.

To use the movie saving feature you'll also need to install ffmpeg.

Iseult should work on Windows, MacOS & Linux, although it may run better on
MacOS, if you replace the first line of `iseult.py` with `#! /usr/bin/env pythonw`

To run Iseult on tigressdata type the following:
```bash
$ module load anaconda3/2019.3
$ cd /path/to/Iseult/
$ chmod +x ./iseult.py
$ ./iseult.py
```
