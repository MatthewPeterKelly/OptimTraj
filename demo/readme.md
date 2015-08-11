# README.md -- TrajOpt/demo

This document gives an overview of the demos, as well as ideas for planned ones in the future.

Each example has a script called MAIN.m, which is the entry-point function for that example.

## simpleWalker
This example uses TrajOpt to find a periodic gait for a simple (double pendulum) model of walking. The model has a continuous torque actuator at the hip, and an impulsive heel-strike (no push-off impulse)

## goddardRocket
Cannonical optimal control problem from rocket engineering. Find the minimal fuel usage to attain some desired height, given atmospheric drag and a simple model of a rocket.

## acrobot
An acrobot is a double pendulum with a motor between the two links. This example finds the minimum torque trajectory to transition the system from hanging down to balancing inverted.

## whyWalkRun - planned
Reproduce the results of Manoj Srinivasavan's and Andy Ruina's paper on modeling cost of transport for a simple model of walking. They discover that walking is more efficient at low speeds, while running is better at high speeds. Neat!

## cartPole - planned
Find the minimal-torque path to swing up a cart-pole.
