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
    (pi/180)*120;  %pendulum angle (wrt gravity)
    0.0;   %horizontal velocity
    0.0];  %pendulum angular rate

tSpan = [0,2];

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

x = z(1,:);
q = z(2,:);
dx = z(3,:);
dq = z(4,:);

u = ctrlFun(z);

[energy, potential, kinetic] = cartPoleEnergy(z,p);

%%%% Plots:
figure(2); clf;

subplot(3,2,1);
plot(t,x)
ylabel('x')

subplot(3,2,2);
plot(t,q)
ylabel('q')

subplot(3,2,3);
plot(t,dx)
ylabel('dx')

subplot(3,2,4);
plot(t,dq)
ylabel('dq')

subplot(3,2,5)
plot(t,u)
ylabel('u')
xlabel('t')

subplot(3,2,6); hold on;
plot(t,potential,'b-','LineWidth',1)
plot(t,kinetic,'r-','LineWidth',1)
plot(t,energy,'k-','LineWidth',2)
ylabel('energy')
xlabel('t')
legend('potential','kinetic','total');

