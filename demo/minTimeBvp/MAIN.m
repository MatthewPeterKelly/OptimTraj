% MAIN - Minimum Time Boundary Value Problem
%
% Solve a minimum-time boundary value problem with simple dynamics (chain
% integrator) and limits on the state and control. Scalar trajectory.
%
% Here we will solve a scalar trajectory, where the position, velocity, 
% and acceleration are states. The jerk (derivative of acceleration) will
% be the only control.
% 

clc; clear;
addpath ../../

% Kinematic Limits:
xLim = [0, 4]; % position
vLim = [-2, 2]; % velocity
aLim = [-4, 4]; % acceleration
jLim = 5*[-8, 8]; % jerk 

% Boundary value problem:
xBegin = xLim(1);  % initial state
vBegin = 0;
aBegin = 0;
xFinal = xLim(2);  % final state
vFinal = 0;
aFinal = 0;

% User-defined dynamics and objective functions
problem.func.dynamics = @(t,x,u)( scalarChainIntegrator(x,u) );
problem.func.bndObj = @(t0,x0,tF,xF)( tF - t0 ); % minimum time  -- primary objective
problem.func.pathObj = @(t,x,u)( 0.001*u.^2 ); %minimum jerk  -- regularization

% Problem boundsTime
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 0.1;
problem.bounds.finalTime.upp = 10;

problem.bounds.state.low = [xLim(1); vLim(1); aLim(1)];
problem.bounds.state.upp = [xLim(2); vLim(2); aLim(2)];
problem.bounds.initialState.low = [xBegin; vBegin; aBegin];
problem.bounds.initialState.upp = [xBegin; vBegin; aBegin];
problem.bounds.finalState.low = [xFinal; vFinal; aFinal];
problem.bounds.finalState.upp = [xFinal; vFinal; aFinal];

problem.bounds.control.low = jLim(1);
problem.bounds.control.upp = jLim(2); 

% Guess at the initial trajectory
problem.guess.time = [0,2];
problem.guess.state = [[xBegin; vBegin; aBegin], [xFinal; vFinal; aFinal]];
problem.guess.control = [0, 0];

% Select a solver:
problem.options(1).method = 'trapezoid';
problem.options(1).trapezoid.nGrid = 8;
problem.options(2).method = 'trapezoid';
problem.options(2).trapezoid.nGrid = 16;
problem.options(3).method = 'hermiteSimpson';
problem.options(3).hermiteSimpson.nSegment = 15;

% Solve the problem
soln = optimTraj(problem);
t = soln(end).grid.time;
q = soln(end).grid.state(1,:);
dq = soln(end).grid.state(2,:);
ddq = soln(end).grid.state(3,:);
u = soln(end).grid.control;

% Plot the solution:
figure(1); clf;

subplot(4,1,1)
plot(t,q)
ylabel('q')
title('Minimum-time boundary value problem');

subplot(4,1,2)
plot(t,dq)
ylabel('dq')

subplot(4,1,3)
plot(t,ddq)
ylabel('ddq')

subplot(4,1,4)
plot(t,u)
ylabel('dddq')


