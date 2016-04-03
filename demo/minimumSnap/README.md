# README.md  --  Pendulum Kinematics

This directory contains a few example problems, each showing how to compute optimal trajectories that contain high derivatives of the angle of the pendulum, hence the name kinematics. 

Since trajectory optimization requires the problem to be in first-order form, we need to use a chain integrator to express these high-derivatives in terms of state variables. This can be a bit confusing, so I've written up this simple example to demonstrate the concept.

## Chain Integrator

This demo shows how to construct a simple chain integrator, for a system with trivial dynamics.

## Minimum Acceleration

This shows the minimum-acceleration trajectory that satisfies the boundary value problem. This is a bit easier to understand than the following demos.


## Minimum Jerk (derivative of acceleration)

This shows the minimum-jerk trajectory that satisfies the boundary value problem. 
This objective function is widely used in trajectory generation for robot arms, althought the formulation is typically different.

## Minimum Snap (second derivative of acceleration)

This shows the minimum-snap trajectory that satisfies the boundary value problem. This objective function is often seen in quad-rotor trajectory generation.