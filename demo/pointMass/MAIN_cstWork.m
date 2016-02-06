% MAIN - Point Mass
%
% Demonstrates how to use slack variables for an objective function that
% includes an abs()
%
% The Integral{abs(power)} cost function is difficult to optimize for two
% reasons:
%   1) The objective function is non-smooth. This is addressed here by
%   introducing a pair of slack variables and a path constraint. An
%   alternative method is shown in MAIN_smoothWork.m, that directly smooths
%   
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

% User-defined dynamics and objective functions
problem.func.dynamics = @(t,x,u)( cstDyn(x,u) );
problem.func.pathObj = @(t,x,u)( obj_cstWork(u) );
problem.func.pathCst = @(t,x,u)( cstSlackPower(x,u) );

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

uMax = 20;
problem.bounds.control.low = [-uMax;zeros(2,1)];  %Two slack variables
problem.bounds.control.upp = [uMax;inf(2,1)];

% Guess at the initial trajectory
problem.guess.time = [0,1];
problem.guess.state = [0, 0; 1, 0];
problem.guess.control = [0, 0;zeros(2,2)]; %Two slack variables

% Options for fmincon
problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'GradObj','on',...
    'GradConstr','on',...
    'DerivativeCheck','off');   %Fmincon automatic gradient check

problem.options.method = 'trapezoid';
problem.options.trapezoid.nGrid = 100;
problem.options.defaultAccuracy = 'medium';

% Solve the problem
soln = trajOpt(problem);
t = soln.grid.time;
q = soln.grid.state(1,:);
dq = soln.grid.state(2,:);
u = soln.grid.control;

% Plot the solution:
figure(3); clf;

subplot(4,1,1)
plot(t,q)
ylabel('pos')
title('Move Point Mass');

subplot(4,1,2)
plot(t,dq)
ylabel('vel')

subplot(4,1,3)
plot(t,u(1,:))
ylabel('force')

subplot(4,1,4);
plot(t,u(2:3,:))
ylabel('slack')



