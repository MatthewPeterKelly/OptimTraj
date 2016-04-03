% MAIN - chain integrator
%
% Problem statement:
%
% Find the minimum-snap trajectory that moves a system between two boundary
% points. Note that snap is the 4th derivative of position. Since the
% dynamics are in first-order form, we need to include position, velocity,
% acceleration, jerk in our state vector. We then set the control to be the
% snap of the trajectory.
%

clc; clear;
addpath ../../..

%%%% Boundary-value problem:

t0 = 0;        %initial time
x0 = [1;0];    %initial position
dx0 = [0;0];   %initial velocity
ddx0 = [0;0];  %initial acceleration
dddx0 = [0;0]; %initial jerk (derivative of acceleration)
z0 = [x0;dx0;ddx0;dddx0];   %Full initial state

tF = 1;        %final time
xF = [0;1];    %final position
dxF = [0;0];   %final velocity
ddxF = [0;0];  %final acceleration
dddxF = [0;0]; %final jerk (derivative of acceleration)
zF = [xF;dxF;ddxF;dddxF];  %full final state


%%%% Construct bounds struct, given problem specifications

problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0;
problem.bounds.finalTime.low = tF;
problem.bounds.finalTime.upp = tF;

problem.bounds.initialState.low = z0;
problem.bounds.initialState.upp = z0;
problem.bounds.finalState.low = zF;
problem.bounds.finalState.upp = zF;


%%%% Construct a simple initial guess (linear between boundary)
problem.guess.time = [t0, tF];
problem.guess.state = [z0, zF];
problem.guess.control = [zeros(size(x0)), zeros(size(xF))];


%%%% Define dynamics and objective functions:

% Enforce the chain integrator dynamics:
problem.func.dynamics = @(t,z,u)(  dynamics(z,u)  );

% Minimize the integral of the snap-squared along the trajectory.
% Sum along each dimension of the state space. 
problem.func.pathObj = @(t,z,u)( sum(u.^2,1) );  


%%%% Select the method of choice:

% problem.options.method = 'trapezoid';
% problem.options.method = 'hermiteSimpson';
problem.options.method = 'chebyshev';
% problem.options.method = 'rungeKutta';
% problem.options.method = 'gpops';    % requires license for GPOPS-II


%%%% Solve!
soln = trajOpt(problem);


%%%% Unpack the solution

tGrid = soln.grid.time;
xGrid = soln.grid.state(1:2, :);
dxGrid = soln.grid.state(3:4, :);
ddxGrid = soln.grid.state(5:6, :);
dddxGrid = soln.grid.state(7:8, :);
ddddxGrid = soln.grid.control;

t = linspace(tGrid(1), tGrid(end), 100);
z = soln.interp.state(t);
x = z(1:2,:);
dx = z(3:4,:);
ddx = z(5:6,:);
dddx = z(7:8,:);
ddddx = soln.interp.control(t);

%%%% Plot the trajectory against time
figure(1); clf;

subplot(5,1,1); hold on;
plot(t,x)
plot(tGrid,xGrid,'ko','MarkerSize',8,'LineWidth',2);

subplot(5,1,2); hold on;
plot(t,dx)
plot(tGrid,dxGrid,'ko','MarkerSize',8,'LineWidth',2);

subplot(5,1,3); hold on;
plot(t,ddx)
plot(tGrid,ddxGrid,'ko','MarkerSize',8,'LineWidth',2);

subplot(5,1,4); hold on;
plot(t,dddx)
plot(tGrid,dddxGrid,'ko','MarkerSize',8,'LineWidth',2);

subplot(5,1,5); hold on;
plot(t,ddddx)
plot(tGrid,ddddxGrid,'ko','MarkerSize',8,'LineWidth',2);


