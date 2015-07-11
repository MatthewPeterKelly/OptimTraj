# TrajOpt - Trajectory Optimization for Matlab
TrajOpt is a matlab library that I developed while working on my PhD at Cornell, studying controller design for walking robots. Some key features of this library:
- __Easy to install -__ no dependancies besides Matlab and the optimization toolbox
- __Easy to use -__ check out the examples to see for yourself
- __Readable source code -__ easy to debug your code and figure out how the software works
- __Rapidly switch methods -__ same function call for a variety of methods, including:
    - trapazoid rule (1st-order direct transcription)
    - hermite-simpson seperated (4th-order direct transcription)
    - multi-segment chebyshev (high-order orthogonal collocation)

# Installation:
1. Clone or download the repository
2. Add the top level folder to your Matlab path
3. Done!

# Usage:
- Call the function `trajOpt` from inside matlab. 
- `trajOpt` takes a single argument: a struct that describes your trajectory optimization problem. 
- `trajOpt` returns a struct that describes the solution. It contains a full description of the problem, the transcription method that was used, and the solution (both as a vector of points and a function handle for interpolation).
- For more details, type `help trajOpt` at the command line, or check out some of the examples in the `demo/` directory.
 
# TrajOpt can solve problems with:
- continuous dynamics (required to at least second-order)
- boundary constraints 
- path constraints
- integral cost function
- boundary cost function