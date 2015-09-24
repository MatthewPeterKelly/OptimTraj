% MAIN - Pendulum
%
% Demonstrates simple swing-up for a single pendulum with a torque motor.
% This is an easy problem, used for demonstrating how to use analytic
% gradients with trajOpt.
%

% Physical parameters of the pendulum
p.m = 1;
p.g = 9.81;
p.l = 0.4;
p.c = 0.1;

% User-defined dynamics and objective functions
problem.func.dynamics = @(t,x,u)( dynamics(t,x,u,p) );
problem.func.pathObj = @(t,x,u)( objective(t,x,u,p) );
problem.func.pathCst = @(t,x,u)( pathConstraint(t,x,u,p) );

% Problem bounds
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 0.5;
problem.bounds.finalTime.upp = 1.5;

problem.bounds.state.low = [-2*pi; -inf];
problem.bounds.state.upp = [2*pi; inf];
problem.bounds.initialState.low = [0;0];
problem.bounds.initialState.upp = [0;0];
problem.bounds.finalState.low = [pi;0];
problem.bounds.finalState.upp = [pi;0];

problem.bounds.control.low = -inf;
problem.bounds.control.upp = inf;

% Guess at the initial trajectory
problem.guess.time = [0,1];
problem.guess.state = [0, pi; pi, pi];
problem.guess.control = [0, 0];

% Options for fmincon
problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e5,...
    'GradObj','on',...
    'GradConstr','off',...
    'DerivativeCheck','on');

problem.options.method = 'trapazoid';
problem.options.defaultAccuracy = 'medium';

% Solve the problem
soln = trajOpt(problem);
t = soln.grid.time;
q = soln.grid.state(1,:);
dq = soln.grid.state(2,:);
u = soln.grid.control;

% Plot the solution:
figure(1); clf;

subplot(3,1,1)
plot(t,q)
ylabel('q')
title('Single Pendulum Swing-Up');

subplot(3,1,2)
plot(t,dq)
ylabel('dq')

subplot(3,1,3)
plot(t,u)
ylabel('u')


