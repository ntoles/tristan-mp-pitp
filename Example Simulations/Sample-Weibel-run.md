---
title: Weibel Simulation
has_children: false
parent: Example Simulations
nav_order: 1
---
Below is an example of a periodic box with counter-moving cold beam, which is unstable to the Weibel instability.


![Bz](https://ntoles.github.io/tristan-mp-pitp/assets/sample_weibelBz.jpg)

![Density](https://ntoles.github.io/tristan-mp-pitp/assets/sample_weibel_dens.jpg)



You can reproduce these results similar to these by following the instructions [here](/GettingStarted/Downloading-and-Compiling-Tristan) and [here](/GettingStarted/Running-your-first-Tristan-MP-simulation) but compiling with the command

`make USER_FILE=user/user_weibel`

When submitting the job, use the input file located at `~/user/input.weibel` and run it using:

`srun ./tristan-mp2d -i input.weibel > out`
