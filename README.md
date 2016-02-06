# TrajOpt - Trajectory Optimization for Matlab
TrajOpt is a matlab library that I developed while working on my PhD at Cornell, studying controller design for walking robots. Some key features of this library:

- __Easy to install -__ no dependencies outside of Matlab (for multiple shooting and direct collocation)
- __Easy to use -__ check out the examples to see for yourself
- __Readable source code -__ easy to debug your code and figure out how the software works
- __Analytic gradients -__ direct collocation methods support analytic gradients
- __Rapidly switch methods -__ choose from a variety of methods:
    - direct collocation
        - trapezoid
        - Hermite-Simpson (seperated)
    - multiple shooting
        - 4th-order Runge-Kutta
    - global (pseudospectral) collocation
        - Chebyshev (Lobatto)

## Installation:
1. Clone or download the repository
2. Add the top level folder to your Matlab path
3. Clone or download [chebfun](http://www.chebfun.org/) (not needed for direct collocation or multiple shooting)
4. Done!


## Usage:
- Call the function `trajOpt` from inside matlab.
- `trajOpt` takes a single argument: a struct that describes your trajectory optimization problem.
- `trajOpt` returns a struct that describes the solution. It contains a full description of the problem, the transcription method that was used, and the solution (both as a vector of points and a function handle for interpolation).
- For more details, type `help trajOpt` at the command line, or check out some of the examples in the `demo/` directory.

## TrajOpt can solve problems with:
- continuous dynamics
- boundary constraints
- path constraints
- integral cost function
- boundary cost function

## Contribute:
This code is still under development, and will be from now until at least May 2016. Please contact me if you have any comments or suggestions, or create a pull request if you would like to add content.

If you are interested in contributing, here are a few possible things to do:
- Create additional demo problems
- Identify holes in the documentation
- Report bugs
- Implement new methods or features

## Contributions:

- __Will Wehner__ wrote the code that enables analytic gradients in the multiple shooting method (4th-order Runge-Kutta). 

