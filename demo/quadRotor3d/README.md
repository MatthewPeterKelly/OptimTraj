# README.md -- 3D Quadcopter

Plant Model:
    Aircraft with arbitrary number of 'thrusters' (one thruster = one motor + one propeller)
    Thrust and torque is calculated based on propeller RPM. 
    See function headers dynQuadRotor3D.m & dynBodyFrame.m for full details.

Dynamics:
    Model of a quadcopter in 3D (technically 6 DOF: X, Y, Z, pitch, roll, yaw).  
    Control is RPM normalized to throttle.

Objective:
    1. Minimum time solution for boundary value problem with control and kinematic limits.
    TODO: add other objective functions

Entry-point script: MAIN.m