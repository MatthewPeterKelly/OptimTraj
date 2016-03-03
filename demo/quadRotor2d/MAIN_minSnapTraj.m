% MAIN  --  Quad-Rotor  --  Minimal-Snap trajectory
%
% Fin the minimal torque-squared trajectory to move the quad-rotor from an
% arbitrary state to the origin.
%

clc; clear;

warning(['This problem appears to be numerically bad, as posed here.',...
    'If you find a bug that is causing the problem, let me know']);

% Dynamics paramters
p.g = 9.81; % (m/s^2) gravity
p.d = 0.3;  % (m) half-width of quad rotor
p.m = 0.2;  % (m) half-mass of the quad rotor

% Trajectory Parameters:
duration = 1;
uMax = 5*p.g*p.m;

% Initial State:
z0 = zeros(3*4,1);  %[pos, vel, acc]
z0(1) = 1.0;  %x0 

% Final State;
zF = zeros(3*4,1);  %[pos, vel, acc]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( integratorChainDynamics(x,u(3:5,:)) );
problem.func.pathObj = @(t,x,u)( sum(u(3:5,:).^2,1) );  %Jerk-Squared cost function
problem.func.pathCst = @(t,x,u)( dynPathCst(x,u,p) );  % Dynamics


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

problem.bounds.control.low = -[uMax*[1;1]; inf(3,1)];
problem.bounds.control.upp = [uMax*[1;1]; inf(3,1)];


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
    'MaxFunEvals',1e5,...
    'MaxIter',50);


% problem.options.method = 'trapezoid'; 
% problem.options.trapezoid.nGrid = 40;

problem.options.method = 'hermiteSimpson';  
problem.options.hermiteSimpson.nSegment = 4;



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = trajOpt(problem);


warning(['This problem appears to be numerically bad, as posed here.',...
    'If you find a bug that is causing the problem, let me know']);


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
figure(3); clf;

subplot(2,2,1); hold on;
plot(t,x);
xlabel('t')
ylabel('x')
title('Minimum snap-squared trajectory')

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


