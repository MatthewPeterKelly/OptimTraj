% MAIN - Point Mass
%
% Finds the optimal trajectory to slide a point-mass across a 1d
% frictionless plane, using a variety of cost functions.
%
% This script optimizes the trajectory using a "smoothed" version of the
% abs() in the objective function. 
%
% The Integral{abs(power)} cost function is difficult to optimize for two
% reasons:
%   1) The objective function is non-smooth. This is addressed here by
%   directly smoothing the objective. The alternative method is to
%   introduce slack variables, as illustrated in MAIN_cstWork.
%
%   2) The second problem is that the solution itself is non-smooth. This
%   means that the piece-wise polynomial representation will fail to
%   accurately represent the solution, making the optimization difficult.
%   One solution to this problem is to add additional smoothing terms to
%   the cost function. Integral of of the input squared is good, as is the
%   integral of the input rate squared.
%
%

clc; clear;
addpath ../../

alpha = 1e0;  %abs() smoothing parameter   1e5 = heavy smoothing,  ~no smoothing 1e-5
beta = 0;   %torque-squared smoothing.

% User-defined dynamics and objective functions
problem.func.dynamics = @(t,x,u)( dynamics(x,u) );
problem.func.pathObj = @(t,x,u)( obj_smoothWork(x,u,alpha, beta) );


% Problem bounds
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 1.0;
problem.bounds.finalTime.upp = 1.0;

problem.bounds.state.low = [0; -inf];
problem.bounds.state.upp = [1; inf];
problem.bounds.initialState.low = [0;0];
problem.bounds.initialState.upp = [0;0];
problem.bounds.finalState.low = [1;0]; 
problem.bounds.finalState.upp = [1;0];

problem.bounds.control.low = -50; %-inf;
problem.bounds.control.upp = 50; %inf;

% Guess at the initial trajectory
problem.guess.time = [0,1];
problem.guess.state = [0, 0; 1, 0];
problem.guess.control = [1, -1];

% Options for fmincon
problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'GradObj','on',...
    'GradConstr','on',...
    'DerivativeCheck','off');   %Fmincon automatic gradient check

problem.options.method = 'trapezoid';
problem.options.trapezoid.nGrid = 40;
problem.options.defaultAccuracy = 'medium';

% Solve the problem
soln = trajOpt(problem);
t = soln.grid.time;
q = soln.grid.state(1,:);
dq = soln.grid.state(2,:);
u = soln.grid.control;

% Plot the solution:
figure(2); clf;

subplot(3,1,1)
plot(t,q)
ylabel('pos')
title('Move Point Mass');

subplot(3,1,2)
plot(t,dq)
ylabel('vel')

subplot(3,1,3)
plot(t,u)
ylabel('force')


