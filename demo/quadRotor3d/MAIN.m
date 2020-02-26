% MAIN - Minimum Time Boundary Value Problem
%
% Solve a minimum-time boundary value problem for a 3D (6 DOF) quadcopter with limits on the state and control. 
%
% The control is the throttle, u, which acts as normalized RPM, where 0 < u < 1 and 0 < RPM < maxRPM for each motor.
% 

clc; clear;
addpath ../../ ./utilities ./test

% Define environmental and plant model params
% Enviromental params
p.g = -9.81 ; % World Coords is XYZ = [East, North, Up], i.e. gravity is a negative number
p.rho = 1.225 ; % air density during flight (kg/m^3) 

% Inertial params
p.m = 5 ; 
p.I = [0.625 0 0; 0 0.625 0; 0 0 1.25] ; % inertia tensor coords: 
p.cg = [0 0 0] ; % (m) location of center of gravity

% Propulsion system params
% propulsion system parameters shared for all motors:
qRP.d_prop = 0.305*ones(4,1) ; % propeller diameter (m)
qRP.maxThrust = 25*ones(4,1) ; % thrust at 100% throttle (N)
qRP.maxRPM = 10000*ones(4,1) ; % RPM at 100% throttle (RPM)
qRP.maxTorque = ones(4,1) ;  % torque at 100% throttle (Nm)
qRP.thrustLocations = [0.5 0 0; 0 0.5 0; -0.5 0 0; 0 -0.5 0]; % motor locations (each row one motor in coords: [port, nose, top] 
qRP.thrustAxes = repmat([0 0 1],4,1) ; % thrust axes of each motor in coords port, nose, top.
qRP.isSpinDirectionCCW = [1; 0; 1; 0] ; % bool to reverse motor spin direction around 'thrustAxes'.
plotflag = 1 ; 
[p.propulsion] = definePropulsionModel(qRP,plotflag); 

% Trajectory Parameters:
uMax = 1 ; % maximum control input; when u = 1; RPM = maxRPM.

% Boundary value problem:
initialState = zeros(12,1) ; % initialize 
finalState = zeros(12,1) ;   % initialize
finalState(1) = 10 ; % assign non-zero state values.

% User-defined dynamics and objective functions
problem.func.dynamics = @(t,x,u)( dynQuadRotor3d(x,u,p) );
problem.func.bndObj = @(t0,x0,tF,xF)( tF - t0 ); % minimum time  -- primary objective
problem.func.pathObj = @(t,x,u)( sum(0.001*u.^2) ); %minimum jerk  -- regularization

% Problem boundsTime
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 1;
problem.bounds.finalTime.upp = 10;

problem.bounds.state.low = -100*ones(size(initialState)) ;
problem.bounds.state.upp = 100*ones(size(initialState)) ; 
problem.bounds.initialState.low = initialState;
problem.bounds.initialState.upp = initialState;
problem.bounds.finalState.low = finalState;
problem.bounds.finalState.upp = finalState;

problem.bounds.control.low = [0;0;0;0] ;
problem.bounds.control.upp = [uMax;uMax;uMax;uMax] ;   

% Guess at the initial trajectory
problem.guess.time = [0,5];
problem.guess.state = [initialState, finalState];
problem.guess.control = ones(4,2);

% Select a solver:
problem.options(1).method = 'trapezoid';
problem.options(1).trapezoid.nGrid = 8;
problem.options(2).method = 'trapezoid';
problem.options(2).trapezoid.nGrid = 16;

% Example syntax to run 'hermiteSimpson' solver.  Can take a while to run:  
% problem.options(3).method = 'hermiteSimpson';
% problem.options(3).hermiteSimpson.nSegment = 15;

% Solve the problem
soln = optimTraj(problem);
t = soln(end).grid.time;
q = soln(end).grid.state(1,:);
dq = soln(end).grid.state(2,:);
ddq = soln(end).grid.state(3,:);
u = soln(end).grid.control ;

% Plot the solution:
plotQuadRotor3d(soln)



