% MAIN  --  minimum jerk* trajectory
%
% *jerk = derivative of acceleration
%
% For a simple pendulum:
%
% x = position
% v = velocity
% u = torque
%
% ddx = f(x,dx,u);     <-- dynamics
%
% cost = integral(  dddx^2  );     <-- cost function
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
% dv2 = a2;
% da2 = u2;   % jerk = derivative of acceleration
% cost = integral(  u2^2  );
%
%
% NOTES:
%   
%   z = [x;v1;v2;a2];
%   u = [u1;u2];
%
clc; clear;

%%%% Specify boundary conditions
t0 = 0;
tF = 1;

z0 = [0;0;0;0];
zF = [pi;0;0;0];

maxTorque = 15;

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

%%%% Options
problem.options.method = 'chebyshev';

%%%% Solve
soln = trajOpt(problem);


%%%% Unpack the solution

tGrid = soln.grid.time;
xGrid = soln.grid.state(1, :);
v1Grid = soln.grid.state(2, :);
v2Grid = soln.grid.state(3, :);
a2Grid = soln.grid.state(4, :);
u1Grid = soln.grid.control(1, :);
da2Grid = soln.grid.control(2, :);

t = linspace(tGrid(1), tGrid(end), 100);
z = soln.interp.state(t);
u = soln.interp.control(t);
x = z(1,:);
v1 = z(2,:);
v2 = z(3,:);
a2 = z(4,:);
u1 = u(1,:);
da2 = u(2,:);

%%%% Plot the trajectory against time
figure(1); clf;

subplot(5,1,1); hold on;
plot(t,x)
plot(tGrid,xGrid,'ko','MarkerSize',8,'LineWidth',2);
title('position (angle)')

subplot(5,1,2); hold on;
plot(t,v1)
plot(t,v2)
plot(tGrid,v1Grid,'ko','MarkerSize',8,'LineWidth',2);
plot(tGrid,v2Grid,'ko','MarkerSize',8,'LineWidth',2);
title('velocity (angular rate)')
legend('v1','v2')

subplot(5,1,3); hold on;
plot(t,a2)
plot(tGrid,a2Grid,'ko','MarkerSize',8,'LineWidth',2);
title('acceleration')

subplot(5,1,4); hold on;
plot(t([1,end]),[1,1]*maxTorque,'k--','LineWidth',1);
plot(t([1,end]),-[1,1]*maxTorque,'k--','LineWidth',1);
plot(t,u1)
plot(tGrid,u1Grid,'ko','MarkerSize',8,'LineWidth',2);
title('torque')

subplot(5,1,5); hold on;
plot(t,da2)
plot(tGrid,da2Grid,'ko','MarkerSize',8,'LineWidth',2);
title('jerk')





