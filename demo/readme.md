# README.txt -- TrajOpt/demo

This directory contains a collection of example problems for trajectory optimization, all solved using TrajOpt. Each example contains a MAIN.m file, which is the entry-point file for the example. 

If the dynamics or constraints are complicated, then a script Derive_*.m is provided to use the symbolic toolbox to derive these equations. Any files with the autoGen_fileName.m are created by the Matlab symbolic toolbox, and should not be edited. In each case, they will be called by a regular function named fileName.m.

If you find any errors, have comments, or would like to suggest an example problem, just send me an email. Contact Info can be found at either page:
-->  https://github.com/MatthewPeterKelly
-->  www.matthewpeterkelly.com


List of Examples:    (Easy --> Hard)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pendulum:
  --> a simple demo showing basic functionality. Shows how to use analytic gradients.

cart-pole:
  --> A standard benchmark for trajectory optimization. Easy to solve for most choices of parameters. 

simpleWalker:
  --> A double-pendulum model for a walking robot. Pretty easy to set up and solve. 

toyCar:
  --> A fun toy problem finding an optimal trajectory for a car driving over a hilly landscape. Turns out to be a bit tricky due to discontinuity in the solution to the problem. 

acrobot:
  --> A standard trajectory optimization problem. Somewhat tricky to solve, but easy to set up.

goddard roacket:
  --> Hard to solve, easy to set up. Optimal rocket thrust trajectory.

five-link biped:
  --> Easy to solve, very difficult to set up. Requires advanced knowledge of dynamics and hybrid systems.

pointMass:
  --> A pathological example. Shows various ways to handle discontinuities.



