# OptimTraj - Trajectory Optimization for Matlab
OptimTraj is a matlab library designed for solving continuous-time single-phase trajectory optimization problems. I developed it while working on my PhD at Cornell, studying non-linear controller design for walking robots.

## What sort of problems does OptimTraj solve?

#### Examples:
- [Cart-pole swing-up](https://youtu.be/kAlhKJlu7O8): Find the force profile to apply to the cart to swing-up the pendulum that freely hanges from it.
- Compute the gait (joint angles, rates, and torques) for a walking robot that minimizes the energy used while walking.
- Find a minimum-thrust orbit transfer trajectory for a satellite.

#### Details:

OptimTraj finds the optimal trajectory for a dynamical system. This trajectory is a sequence of controls (expressed as a function) that moves the dynamical system between two points in state space. The trajectory will minimize some cost function, which is typically an integral along the trajectory. The trajectory will also satisfy a set user-defined constraints.

OptimTraj solves problems with
- continuous dynamics
- boundary constraints
- path constraints
- integral cost function
- boundary cost function

All functions in the problem description can be non-linear, but they must be smooth (C2 continuous).


## Features:

- __Easy to install -__ no dependencies outside of Matlab (for base functionality)
- __Lots of examples -__ look at the `demo/` directory to see for yourself!
- __Readable source code -__ easy to debug your code and figure out how the software works
- __Analytic gradients -__ most methods support analytic gradients
- __Rapidly switch methods -__ choose from a variety of methods:
    - direct collocation
        - trapezoid
        - Hermite-Simpson (seperated)
    - direct multiple shooting
        - 4th-order Runge-Kutta
    - global (pseudospectral) collocation
        - Chebyshev (Lobatto)  --  (requires [chebfun](http://www.chebfun.org/))

## Installation:
1. Clone or download the repository
2. Add the top level folder to your Matlab path
3. (Optional) Clone or download [chebfun](http://www.chebfun.org/) (needed for global collocation)
4. Done!


## Usage:
- Call the function `optimTraj` from inside matlab.
- `optimTraj` takes a single argument: a struct that describes your trajectory optimization problem.
- `optimTraj` returns a struct that describes the solution. It contains a full description of the problem, the transcription method that was used, and the solution (both as a vector of points and a function handle for interpolation).
- For more details, type `help optimTraj` at the command line, or check out some of the examples in the `demo/` directory.

## Contribute:
This code is still under development, and will be from now until at least May 2016. Please contact me if you have any comments or suggestions, or create a pull request if you would like to add content.

If you are interested in contributing, here are a few possible things to do:
- Create additional demo problems
- Identify holes in the documentation
- Report bugs
- Implement new methods or features

## Contributions:

- [__Will Wehner__](https://github.com/wwehner) wrote the code that enables analytic gradients in the multiple shooting method (4th-order Runge-Kutta).
