%MAIN.m  --  simple walker trajectory optimization
%
% This script sets up a trajectory optimization problem, and then solves it
% using trajOpt.
%
%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Parameters for the dynamics function                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
param.dyn.m1 = 10;  %hip mass
param.dyn.m2 = 1;  %foot mass
param.dyn.g = 9.81;  %gravity
param.dyn.l = 1;  %leg length


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( dynamics(x,u,param.dyn) );

problem.func.pathObj = @(t,x,u)( costFun(u) );

problem.func.endObj = [];

problem.func.pathCst = [];

problem.func.endCst = [];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
t0 = 0;  tF = 1;
problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0;
problem.bounds.finalTime.low = tF;
problem.bounds.finalTime.upp = tF;

% State: [q1;q2;dq1;dq2];

problem.bounds.state.low = [-pi/2; -pi/2; -4; -4];
problem.bounds.state.upp = [ pi/2;  pi/2;  4;  4];

stepAngle = 0.3;
problem.bounds.initialState.low = [ stepAngle; -stepAngle; -4; -4];
problem.bounds.initialState.upp = [ stepAngle; -stepAngle;  4;  4];
problem.bounds.finalState.low = [-stepAngle;  stepAngle; -4; -4];
problem.bounds.finalState.upp = [-stepAngle;  stepAngle;  4;  4];

problem.bounds.control.low = -inf;
problem.bounds.control.upp =  inf;



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% For now, just assume a linear trajectory between boundary values

problem.guess.time = [t0, tF];

x0 = 0.5*(problem.bounds.initialState.low + problem.bounds.initialState.upp);
xF = 0.5*(problem.bounds.finalState.low + problem.bounds.finalState.upp);
problem.guess.state = [x0, xF];

problem.guess.control = [0, 0];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.nlpOpt = [];   %Use default options for fmincon

problem.options.method = 'trapazoid';
problem.options.nGrid = 20;   


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = trajOpt(problem);

