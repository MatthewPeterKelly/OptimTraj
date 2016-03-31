% MAIN  --  Quad-Rotor  --  Minimal-Force trajectory
%
% Fin the minimal torque-squared trajectory to move the quad-rotor from an
% arbitrary state to the origin.
%
% NOTES:
%   X = [x;y;q] = [x pos, y pos, angle] = configuration
%  dX = [dx;dy;dq] = [x vel, y vel, angle rate] = rate
% ddX = [ddx;ddy;ddq] = acceleration
%

clc; clear;

addpath ../../

% Dynamics paramters
p.g = 9.81; % (m/s^2) gravity
p.d = 0.3;  % (m) half-width of quad rotor
p.m = 0.2;  % (m) half-mass of the quad rotor

% Trajectory Parameters:
duration = 1;
uMax = 5*p.g*p.m;

% Initial State:
X0 = [1;0;0];   %  initial configuration
dX0 = zeros(3,1);  % initial rates
z0 = [X0; dX0];  % initial state

XF = [0;0;0];   % final configuration
dXF = zeros(3,1);  % final rates
zF = [XF; dXF];  % final state

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,z,u)( dynamics(z,u,p) );
problem.func.pathObj = @(t,z,u)( sum(u.^2,1) );  %Force-squared cost function


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

problem.bounds.control.low = -uMax*[1;1];
problem.bounds.control.upp = uMax*[1;1];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Initial guess at trajectory                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.guess.time = [0,duration];
problem.guess.state = [z0, zeros(6,1)];
problem.guess.control = p.g*p.m*ones(2,2);



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e5);


% problem.options.method = 'trapezoid'; 
% problem.options.trapezoid.nGrid = 40;

problem.options.method = 'hermiteSimpson';  
problem.options.hermiteSimpson.nSegment = 30;



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
dx = z(4,:);
dy = z(5,:);
dq = z(6,:);

u = soln.interp.control(t);
u1 = u(1,:);
u2 = u(2,:);


%%%% Plots:
figure(1); clf;

subplot(2,2,1); hold on;
plot(t,x);
xlabel('t')
ylabel('x')
title('Minimum force-squared trajectory')

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


