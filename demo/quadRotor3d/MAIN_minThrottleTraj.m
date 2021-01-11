% MAIN  --  Quad-Rotor  --  Minimal-Throttle trajectory
%
% Find the minimal throttle-squared trajectory to move the quad-rotor from an
% arbitrary state to the origin.
%
% NOTES:
%   X = [x;y;z;p;r;w] = [x pos, y pos, z pos, pitch att, roll att, yaw att] = configuration
%  dX = [dx;dy;dz;dp;dr;dw] = [x vel, y vel, z vel, pitch rate, roll rate, yaw rate] = rate
% ddX = [ddx;ddy;ddz;ddp;ddr;ddw] = acceleration
%

clc; clear;

addpath ../../ ./utilities ./test

% Define environmental and plant model params
[p] = loadPlant_QuadRotor3d(); 

% Trajectory Parameters:
duration = 3;

% Initial State:
X0 = [1;0;0;0;0;0];   % initial configuration
dX0 = zeros(6,1);     % initial rates
z0 = [X0; dX0];  % initial state

XF = [0;0;0;0;0;0];   % final configuration
dXF = zeros(6,1);     % final rates
zF = [XF; dXF];  % final state

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,z,u)( dynQuadRotor3d(z,u,p) );
problem.func.pathObj = @(t,z,u)( sum(u.^2,1) );  % Throttle-squared cost function


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

problem.bounds.control.low = -p.uMax*[1;1;1;1];
problem.bounds.control.upp = p.uMax*[1;1;1;1];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Initial guess at trajectory                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.guess.time = [0,duration];
problem.guess.state = [z0, zeros(12,1)];
problem.guess.control = ones(4,2);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e5);

problem.options.method = 'trapezoid'; 
problem.options.trapezoid.nGrid = 16;

% problem.options.method = 'hermiteSimpson';  
% problem.options.hermiteSimpson.nSegment = 30;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = optimTraj(problem);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Display Solution                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Plot the solution:
plotQuadRotor3d(soln)

