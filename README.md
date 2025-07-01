
# Origami Simulation with Collision and Friction

This is a physics based Origami simulation aimed at simulating the behaviour of actuated creases including approximating contact and friction.
This project is written as a Bachelor's thesis with the [Computational Design Laboratory](https://cdl.ethz.ch/) @ ETH.  
It takes [.fold](https://github.com/edemaine/fold) files as [Crease Pattern](https://en.wikipedia.org/wiki/Crease_pattern) input, as well as .json parameter files and actuation profiles.
Much of the mathematics behind this simulation are based on [Origami Simulator](https://github.com/amandaghassaei/OrigamiSimulator).

## How to run locally

Use the following commands to set up this simulation locally:

```
mkdir build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

and then call the executable with the arguments `./origami-simulator <Crease Pattern> <Parameters> <Actuation Profile>`
