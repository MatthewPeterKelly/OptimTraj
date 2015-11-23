% TEST_openLoopDynamics.m
%
% This script performs some basic checks on the equations of motion.
%
% For example, the total energy should be constant to the tolerance of the
% integrator if the applied torque (u) is zero.
%
% If m1 >> m2, then q should behave like a simple pendulum
%

clc; clear;

%%%% Set up the simulation
z0 = [
    0.0;   %horizontal position
    (pi/180)*80;  %pendulum angle (wrt gravity)
    0.3;   %horizontal velocity
    0.5];  %pendulum angular rate

tSpan = [0,1.5];

p.m1 = 1.0;  % (kg) Cart mass
p.m2 = 0.3;  % (kg) pole mass
p.g = 9.81;  % (m/s^2) gravity 
p.l = 0.5;   % (m) pendulum (pole) length 

%%%% Function Handles
ctrlFun = @(z)( zeros(size(z(1,:))) );  %Passive controller for now
dynFun = @(t,z)( cartPoleDynamics(z, ctrlFun(z), p) );

%%%% Simulate the system!
options = odeset(...
    'RelTol',1e-8, ...
    'AbsTol',1e-8);
sol = ode45(dynFun, tSpan, z0, options);

%%%% Unpack the simulation
t = linspace(tSpan(1), tSpan(2), 200);
z = deval(sol,t);
u = ctrlFun(z);

%%%% Plots:
figure(1); clf;
plotPendulumCart(t,z,u,p);


%%%% Draw Trajectory:
[p1,p2] = cartPoleKinematics(z,p);

figure(2); clf; 
nFrame = 5;  %Number of frames to draw
drawCartPoleTraj(t,p1,p2,nFrame);




