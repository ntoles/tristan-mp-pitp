Below is an example of a periodic box with counter-moving cold beam, which are unstable to the Weibel instability. 

[[_img/sample_weibelBz.jpg|alt="Bz"]]

[[_img/sample_weibel_dens|alt="Density"]]

You can reproduce these results similar to these by following the instructions [here](https://github.com/PrincetonUniversity/tristan-mp-pu/wiki/Downloading-and-Compiling-Tristan) and [here](https://github.com/PrincetonUniversity/tristan-mp-pu/wiki/Running-your-first-Tristan-MP-simulation) but compiling with the command

`make USER_FILE=user/user_weibel`

When submitting the job, use the input file located at `~/user/input.weibel` and run it using:

`srun ./tristan-mp2d -i input.weibel > out`
