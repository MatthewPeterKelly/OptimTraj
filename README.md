# OptimTraj - Trajectory Optimization for Matlab

[![DOI](https://zenodo.org/badge/40544279.svg)](https://zenodo.org/badge/latestdoi/40544279)

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

## Documentation:

The best way to learn OptimTraj is by working through a few of the examples in the `demo/` directory. There is also a `.pdf` user guide that can be found in the `/docs` directory.

For more background on trajectory optimization in general, I wrote a [tutorial paper](https://epubs.siam.org/doi/10.1137/16M1062569) that includes derivations and references for most methods implemented here, along with a variety of practical suggestions and debugging tips. Finally, I have a [tutorial webpage](http://www.matthewpeterkelly.com/tutorials/trajectoryOptimization/index.html) for trajectory optimization.

## Citing OptimTraj

The best way to cite OptimTraj is using the DOI assigned the the library source code:

[![DOI](https://zenodo.org/badge/40544279.svg)](https://zenodo.org/badge/latestdoi/40544279)

```bibtex
@software{Kelly_OptimTraj_Trajectory_Optimization_2022,
  author = {Kelly, Matthew Peter},
  doi = {10.5281/zenodo.7430524},
  month = {12},
  title = {{OptimTraj: Trajectory Optimization for Matlab}},
  url = {https://github.com/MatthewPeterKelly/OptimTraj},
  version = {1.7},
  year = {2022}
}
```

Alternatively, you can cite the [tutorial paper](https://epubs.siam.org/doi/10.1137/16M1062569) that accompanies OptimTraj:
```bibtex
@article{kelly2017introduction,
  title={An Introduction to Trajectory Optimization: How to do your own Direct Collocation},
  author={Kelly, Matthew},
  journal={SIAM Review},
  volume={59},
  number={4},
  pages={849--904},
  year={2017},
  publisher={SIAM}
}
```

## Contribute:
Contributions are welcome! Feel free to reach out if you're planning a bigger submission, or just go ahead and make a pull request for smaller contributions.

If you are interested in contributing, here are a few possible things to do:
- Create additional demo problems
- Identify holes in the documentation
- Report bugs
- Implement new methods or features
- Look through the open issues

## Contributions:

- [__Will Wehner__](https://github.com/wwehner) implemented analytic gradients in the 4th-order Runge-Kutta method.

- [__Conrad McGreal__](https://github.com/cjmcgreal) contributed the minimum-time 3D quadrotor example.
