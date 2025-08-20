
# Origami Simulation with Collision and Friction

This is a physics based Origami simulation aimed at simulating the behaviour of actuated creases including approximating contact and friction.
This project is written as a Bachelor's thesis with the [Computational Design Laboratory](https://cdl.ethz.ch/) @ ETH.  
It takes [.fold](https://github.com/edemaine/fold) files as [Crease Pattern](https://en.wikipedia.org/wiki/Crease_pattern) input, as well as .json parameter files and actuation profiles.
Much of the mathematics behind this simulation are based on [Origami Simulator](https://github.com/amandaghassaei/OrigamiSimulator).

## How to run locally

Use the following commands to set up this simulation locally:

```
mkdir release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```
Good Exports:
ffmpeg -framerate 30 -i %05d.png -c:v libx264 -crf 18 -preset slow -pix_fmt yuv420p ../myout.mp4

and then call the executable with the arguments `./origami-simulator <Crease Pattern> <Parameters> <Actuation Profile>`
