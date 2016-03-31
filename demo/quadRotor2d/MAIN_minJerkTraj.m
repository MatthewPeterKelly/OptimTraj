% MAIN  --  Quad-Rotor  --  Minimal-Jerk trajectory
%
% Fin the minimal jerk-squared trajectory to move the quad-rotor from an
% arbitrary state to the origin. Note that jerk is the derivative of
% acceleration.
%
% NOTES:
%   X = [x;y;q] = [x pos, y pos, angle] = configuration
%  dX = [dx;dy;dq] = [x vel, y vel, angle rate] = rate
% ddX = [ddx;ddy;ddq] = acceleration
%
% PROBLEM:
% 
% ddX = f(X,dX,u);     <-- dynamics
%
% cost = integral(  dddX^2  );     <-- cost function
%
% subject to:
%   X(0) = X0;
%   X(1) = XF;
%   dX(0) = dX0;
%   dX(1) = dXF;
%
% How to pose as a standard trajectory optimization problem?
%
% dX = V1;
% dV1 = f(X,V1,U1)
%
% V2 == V1;   % <-- Key line. 
% dV2 = A2
% dA2 = U2;
% cost = integral(  U2^2  );
%
% z = [X;V1;V2;A2]
% u = [U1;U2]
%


clc; clear;

addpath ../../

% Dynamics paramters
p.g = 9.81; % (m/s^2) gravity
p.d = 0.3;  % (m) half-width of quad rotor
p.m = 0.2;  % (m) half-mass of the quad rotor

% Trajectory Parameters:
duration = 2;
uMax = 5*p.g*p.m;

% Initial State:
X0 = [1;0;0];   %  initial configuration
dX0 = zeros(3,1);  % initial rates
ddX0 = zeros(3,1);  % initial acceleration
z0 = [X0; dX0; dX0; ddX0];  % initial state

XF = [0;0;0];   % final configuration
dXF = zeros(3,1);  % final rates
ddXF = zeros(3,1);  % final acceleration
zF = [XF; dXF; dXF; ddXF];  % final state

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,z,u)( dynJerk(z,u,p) );
problem.func.pathObj = @(t,z,u)( pathObj(u) );  %accel-squared cost function
problem.func.pathCst = @(t,z,u)( pathCst(z) );

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up problem bounds                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = duration;
problem.bounds.finalTime.upp = duration;

problem.bounds.initialState.low = z0;
problem.bounds.initialState.upp = z0;
problem.bounds.finalState.low = zF;
problem.bounds.finalState.upp = zF;

problem.bounds.control.low = [-uMax*[1;1];  -inf(3,1)];   %[torque, accel]
problem.bounds.control.upp = [uMax*[1;1];  inf(3,1)];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Initial guess at trajectory                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.guess.time = [0,duration];
problem.guess.state = [z0, zF];
problem.guess.control = [p.g*p.m*ones(2,2); zeros(3,2)];



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e5);


% problem.options.method = 'trapezoid'; 
% problem.options.trapezoid.nGrid = 40;
% 
% problem.options.method = 'hermiteSimpson';  
% problem.options.hermiteSimpson.nSegment = 10;

problem.options.method = 'chebyshev';
problem.options.chebyshev.nColPts = 15;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = trajOpt(problem);



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Display Solution                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Unpack the simulation
t = linspace(soln.grid.time(1), soln.grid.time(end), 150);

z = soln.interp.state(t);
x = z(1,:);
y = z(2,:);
q = z(3,:);
X1 = z(1:3,:);
V1 = z(4:6,:);
V2 = z(7:9,:);
A2 = z(10:12,:);

u = soln.interp.control(t);
u1 = u(1,:);
u2 = u(2,:);
J2 = u(3:5,:);

[dObj,uStar] = pathObj(u);


%%%% Plots:

figure(2); clf;

subplot(2,2,1)
plot(t,X1);
legend('x','y','q')
title('configuration')

subplot(2,2,3)
plot(t,V2);
legend('x','y','q')
title('rates')

subplot(2,2,2)
plot(t,A2);
legend('x','y','q')
title('acceleration')

subplot(2,2,4)
plot(t,J2);
legend('x','y','q')
title('jerk')



% Configuration trajectories
figure(1); clf;

subplot(2,2,1); hold on;
plot(t,x);
xlabel('t')
ylabel('x')
title('Minimum acceleration-squared trajectory')

subplot(2,2,2); hold on;
plot(t,y);
xlabel('t')
ylabel('y')

subplot(2,2,3); hold on;
plot(t,q);
xlabel('t')
ylabel('q')

subplot(2,2,4); hold on;
plot(t,u1);  plot(t,u2);
xlabel('t')
ylabel('u')
legend('u1','u2');


