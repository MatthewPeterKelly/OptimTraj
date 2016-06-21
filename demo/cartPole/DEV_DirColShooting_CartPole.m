% MAIN.m
%
% Solve the cart-pole swing-up problem

clc; clear;
addpath ../../

p.m1 = 2.0;  % (kg) Cart mass
p.m2 = 0.5;  % (kg) pole mass
p.g = 9.81;  % (m/s^2) gravity
p.l = 0.5;   % (m) pendulum (pole) length

dist = 0.8;  %How far must the cart translate during its swing-up
maxForce = 100;  %Maximum actuator forces
duration = 2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( cartPoleDynamics(x,u,p) );
problem.func.pathObj = @(t,x,u)( u.^2 );  %Force-squared cost function

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up problem bounds                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = duration;
problem.bounds.finalTime.upp = duration;

problem.bounds.initialState.low = zeros(4,1);
problem.bounds.initialState.upp = zeros(4,1);
problem.bounds.finalState.low = [dist;pi;0;0];
problem.bounds.finalState.upp = [dist;pi;0;0];

problem.bounds.state.low = [-2*dist;-2*pi;-inf;-inf];
problem.bounds.state.upp = [2*dist;2*pi;inf;inf];

problem.bounds.control.low = -maxForce;
problem.bounds.control.upp = maxForce;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Initial guess at trajectory                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.guess.time = [0,duration];
problem.guess.state = [problem.bounds.initialState.low, problem.bounds.finalState.low];
problem.guess.control = [0,0];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e5);

% problem.options.method = 'trapezoid';
% problem.options.method = 'hermiteSimpson';
% problem.options.method = 'rungeKutta';
% problem.options.method = 'chebyshev';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = optimTraj(problem);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Display Solution                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Unpack the simulation
t = linspace(soln.grid.time(1), soln.grid.time(end), 150);
z = soln.interp.state(t);
u = soln.interp.control(t);

%%%% Plots:

%%%% Draw Trajectory:
[p1,p2] = cartPoleKinematics(z,p);

figure(2); clf;
nFrame = 9;  %Number of frames to draw
drawCartPoleTraj(t,p1,p2,nFrame);


%%%% Show the error in the collocation constraint between grid points:
%
if strcmp(soln.problem.options.method,'trapezoid') || strcmp(soln.problem.options.method,'hermiteSimpson')
    % Then we can plot an estimate of the error along the trajectory
    figure(5); clf;
    
    % NOTE: the following commands have only been implemented for the direct
    % collocation(trapezoid, hermiteSimpson) methods, and will not work for
    % chebyshev or rungeKutta methods.
    cc = soln.interp.collCst(t);
    
    subplot(2,2,1);
    plot(t,cc(1,:))
    title('Collocation Error:   dx/dt - f(t,x,u)')
    ylabel('d/dt cart position')
    
    subplot(2,2,3);
    plot(t,cc(2,:))
    xlabel('time')
    ylabel('d/dt pole angle')
    
    idx = 1:length(soln.info.error);
    subplot(2,2,2); hold on;
    plot(idx,soln.info.error(1,:),'ko');
    title('State Error')
    ylabel('cart position')
    
    subplot(2,2,4); hold on;
    plot(idx,soln.info.error(2,:),'ko');
    xlabel('segment index')
    ylabel('pole angle');
end

%%%% Plot the state and control against time
figure(1); clf;
plotPendulumCart(t,z,u,p);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% MAIN - Pendulum
%
% Demonstrates simple swing-up for a single pendulum with a torque motor.
% This is an easy problem, used for demonstrating how to use analytic
% gradients with trajOpt.
%

clc; clear;
addpath ../../

% Physical parameters of the pendulum
p.k = 1;  % Normalized gravity constant
p.c = 0.1;  % Normalized damping constant

% User-defined dynamics and objective functions
problem.func.dynamics = @(t,x,u)( dynamics(x,u,p) );
problem.func.pathObj = @(t,x,u)( pathObjectiveTest(t,x,u) );

% bound objective for testing
xF_target = [pi;0];
problem.func.bndObj = @(t0,x0,tF,xF)( boundObjective(xF,xF_target) );

% Problem bounds
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 0.5;
problem.bounds.finalTime.upp = 2.5;

problem.bounds.state.low = [-2*pi; -inf];
problem.bounds.state.upp = [2*pi; inf];
problem.bounds.initialState.low = [0;0];
problem.bounds.initialState.upp = [0;0];
problem.bounds.finalState.low = [pi;0];
problem.bounds.finalState.upp = [pi;0];

problem.bounds.control.low = -5; %-inf;
problem.bounds.control.upp = 5; %inf;

% Guess at the initial trajectory
problem.guess.time = [0,1];
problem.guess.state = [0, pi; pi, 0];
problem.guess.control = [1, 0];

% Options for fmincon
problem.options(1).nlpOpt = optimset(...
    'Display','iter',...
    'GradObj','on',...
    'GradConstr','on',...
    'DerivativeCheck','on',...
    'MaxFunEvals',6000);   %Fmincon automatically checks derivatives

method =  'trapezoid';
problem.options(1).method = method;
problem.options(1).shooting = 'on';
problem.options(1).defaultAccuracy = 'medium';
problem.options(1).crtldefect = 'on';
% problem.options(1).(method).shooting = 'on';
% problem.options(1).(method).defaultAccuracy = 'medium';
% problem.options(1).(method).crtldefect = 'off';

problem.options(2) = problem.options(1);


% Solve the problem
soln = optimTraj(problem);
soln = soln(end);
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


