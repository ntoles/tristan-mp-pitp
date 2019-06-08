Below is an example of a periodic box with counter-moving cold beam, which are unstable to the Weibel instability. 

<img src="_img/sample_weibelBz.jpg" alt="Bz" class="inline"/>
<img src="_img/sample_weibel_dens.jpg" alt="Density" class="inline"/>


You can reproduce these results similar to these by following the instructions [here](Downloading-and-Compiling-Tristan) and [here](Running-your-first-Tristan-MP-simulation) but compiling with the command

`make USER_FILE=user/user_weibel`

When submitting the job, use the input file located at `~/user/input.weibel` and run it using:

`srun ./tristan-mp2d -i input.weibel > out`
