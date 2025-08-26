
# Origami Simulation with Collision and Friction

This is a physics based Origami simulation aimed at simulating the behavior of actuated creases including approximating contact and friction.
This project is written as a Bachelor's thesis with the [Computational Design Laboratory](https://cdl.ethz.ch/) at ETH.  
It takes [.fold](https://github.com/edemaine/fold) files as [Crease Pattern](https://en.wikipedia.org/wiki/Crease_pattern) input, as well as .json parameter files and actuation profiles.
Much of the mathematics behind this simulation is based on [Origami Simulator](https://github.com/amandaghassaei/OrigamiSimulator) and [Incremental Potential Contact](https://ipc-sim.github.io/).

## How to run locally

Using CMakeLists.txt you can execute the following commands to run the simulation locally:

```
mkdir release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

## Simulation Arguments

The simulation takes the following arguments:
`./origami-simulator <Crease Pattern>(.fold format) <Parameters>(.json) <Fold Timeline>(.json)`
The argument priority goes from left to right. So if only two arguments are provided, it's assumed to be Crease Pattern and Parameters. Any missing argument will use a default instead.

## Video Export

To start the simulation, press "s". To export the current simulation frames press "f".

Use the following command inside `results/tmp` to export the current simulation frames to video:
`ffmpeg -framerate 30 -i %05d.png -c:v libx264 -crf 18 -preset slow -pix_fmt yuv420p ../myout.mp4`

## Python Helper scripts

These scripts can be found in src/Python-Helperscripts.

### 3d_print_maker

This script takes as input a hard-coded path to a .fold file and creates a 3D-printable .obj file made up of three layers, a bottom and top PLA layer and a middle TPU layer.

### fold-labeller

This script takes as input a hard-coded path to a .fold file and outputs a .png with all the edges labelled. Helpful for writing fold-timelines.

### function-timeline-maker

This script can generate Fold Timelines (or parts of a fold timeline) where certain folds are moving as a sine-wave.

### strain-plotter

Takes as input a path to a strain file (per default these are saved at data/simulation_strains) and generates a strain plot from it
