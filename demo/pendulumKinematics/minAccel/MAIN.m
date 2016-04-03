% MAIN  --  minimum acceleration trajectory
%
% For a simple pendulum:
%
% x = position
% v = velocity
% u = torque
%
% ddx = f(x,dx,u);     <-- dynamics
%
% cost = integral(  ddx^2  );     <-- cost function
%
% subject to:
%   x(0) = 0;
%   x(1) = pi;
%   dx(0) = 0;
%   dx(1) = pi;
%
% How to pose as a standard trajectory optimization problem?
%
% dx = v1;
% dv1 = f(x,v1,u1)
%
% v2 == v1;   % <-- Key line. 
% dv2 = u2;
% cost = integral(  u2^2  );
%
%
% NOTES:
%   
%   z = [x;v1;v2];
%   u = [u1;u2];
%
clc; clear;
addpath ../../..

%%%% Specify boundary conditions
t0 = 0;
tF = 5;    

maxTorque = 1.2;
z0 = [0;0;0];
zF = [pi;0;0];


%%%% Pack up boundary conditions
problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0;

problem.bounds.finalTime.low = tF;
problem.bounds.finalTime.upp = tF;

problem.bounds.initialState.low = z0;
problem.bounds.initialState.upp = z0;

problem.bounds.finalState.low = zF;
problem.bounds.finalState.upp = zF;

problem.bounds.control.low = [-maxTorque; -inf];
problem.bounds.control.upp = [maxTorque; inf];

%%%% Initialize trajectory with a straight line
problem.guess.time = [t0,tF];
problem.guess.state = [z0, zF];
problem.guess.control = [zeros(2,1), zeros(2,1)];

%%%% Pack up function handles
problem.func.dynamics = @(t,z,u)(  dynamics(z,u)  );
problem.func.pathObj = @(t,z,u)(  pathObjective(u)  );
problem.func.pathCst = @(t,z,u)(  pathConstraint(z)  );

%%%% Choice of solver:
method = 'trapezoid';

switch method
    case 'chebyshev'
        problem.options.method = method;
        problem.options.chebyshev.nColPts = 25;
    case 'trapezoid'
        problem.options.method = method;
    case 'hermiteSimpson'
        problem.options.method = method;
        problem.options.hermiteSimpson.nSegment = 15;
        problem.options.nlpOpt.MaxFunEvals = 5e4;
    case 'gpops'
        problem.options.method = 'gpops';
    otherwise
        error('invalid method')
end

%%%% Solve
soln = trajOpt(problem);


%%%% Unpack the solution

tGrid = soln.grid.time;
xGrid = soln.grid.state(1, :);
v1Grid = soln.grid.state(2, :);
v2Grid = soln.grid.state(3, :);
u1Grid = soln.grid.control(1, :);
dv2Grid = soln.grid.control(2, :);

t = linspace(tGrid(1), tGrid(end), 100);
z = soln.interp.state(t);
u = soln.interp.control(t);
x = z(1,:);
v1 = z(2,:);
v2 = z(3,:);
u1 = u(1,:);
dv2 = u(2,:);

%%%% Plot the trajectory against time
figure(1); clf;

subplot(2,2,1); hold on;
plot(t,x)
plot(tGrid,xGrid,'ko','MarkerSize',8,'LineWidth',2);
title('angle')

subplot(2,2,2); hold on;
plot(t,v1)
plot(t,v2)
plot(tGrid,v1Grid,'ko','MarkerSize',8,'LineWidth',2);
plot(tGrid,v2Grid,'ko','MarkerSize',8,'LineWidth',2);
title('angular rate')
legend('v1','v2')

subplot(2,2,3); hold on;
plot(t([1,end]),[1,1]*maxTorque,'k--','LineWidth',1);
plot(t([1,end]),-[1,1]*maxTorque,'k--','LineWidth',1);
plot(t,u1)
plot(tGrid,u1Grid,'ko','MarkerSize',8,'LineWidth',2);
title('torque')

subplot(2,2,4); hold on;
plot(t,dv2)
plot(tGrid,dv2Grid,'ko','MarkerSize',8,'LineWidth',2);
title('angular acceleration')





