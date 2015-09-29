# Five-Link Biped - Example problem for TrajOpt

The model for this problem is for the robot RABBIT, taken from the 2003 paper by Westervelt, Grizzle, and Koditschek: "Hybrid Zero Dynamics of Planar Biped Walkers"

The dynamics are derived here in the file `Derive_Equations.m` using the matlab symbolic toolbox. This script also generates the gradients of all functions used by the trajectory optimization.

The entry-point script for this example is `MAIN.m` which can be used to run the trajectory optimization using a variety of techniques.

The default objective function is a simple sum of the integral of torque-squared at each joint. This produces nice smooth trajectories. I've also included two other objective functions, which can be used instead by simply commenting out the correct line in MAIN.m. They produce trajectories that are not as nice, tending to drag the swing foot along the ground. These can be fixed by adding more constraints to the solution.
