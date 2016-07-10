%MAIN.m  --  solve swing-up problem for acrobot
%
% This script finds the minimum torque-squared trajectory to swing up the
% acrobot robot: a double pendulum with a motor between the links
%
%

clc; clear;
addpath ../../

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Parameters for the dynamics function                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dyn.m1 = 1;  % elbow mass
dyn.m2 = 1; % wrist mass
dyn.g = 9.81;  % gravity
dyn.l1 = 0.5;   % length of first link
dyn.l2 = 0.5;   % length of second link

t0 = 0;
tF = 2.0;  %For now, force it to take exactly this much time.
x0 = [0;0];   %[q1;q2];  %initial angles   %Stable equilibrium
xF = [pi;pi];  %[q1;q2];  %final angles    %Inverted balance
dx0 = [0;0];   %[dq1;dq2];  %initial angle rates
dxF = [0;0];  %[dq1;dq2];  %final angle rates
maxTorque = 20;  % Max torque at the elbow  (GPOPS goes crazy without this)

%  * The optimal trajectory is not actually constrained by the maximum
%  torque. That being said, GPOPS goes numerically unstable if the torque
%  is not bounded. This does not seem to be a problem with the other
%  methods.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( acrobotDynamics(x,u,dyn) );

problem.func.pathObj = @(t,x,u)( u.^2 );  %Simple torque-squared

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0;
problem.bounds.finalTime.low = tF;
problem.bounds.finalTime.upp = tF;

% State: [q1;q2;dq1;dq2];

problem.bounds.state.low = [-2*pi; -2*pi; -inf(2,1)];
problem.bounds.state.upp = [ 2*pi;  2*pi;  inf(2,1)];

problem.bounds.initialState.low = [x0; dx0];
problem.bounds.initialState.upp = [x0; dx0];
problem.bounds.finalState.low = [xF; dxF];
problem.bounds.finalState.upp = [xF; dxF];

problem.bounds.control.low = -maxTorque;
problem.bounds.control.upp = maxTorque;



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%%%% Run the optimization twice: once on a rough grid with a low tolerance,
%%%% and then again on a fine grid with a tight tolerance.

method = 'trapezoid'; %  <-- this is robust, but less accurate
% method = 'direct'; %  <-- this is robust, but some numerical artifacts
% method = 'rungeKutta';  % <-- slow, gets a reasonable, but sub-optimal soln
% method = 'orthogonal';    %  <-- this usually finds bad local minimum
% method = 'gpops';      %  <-- fast, but numerical problem is maxTorque is large

switch method
    case 'direct'
        problem.options(1).method = 'trapezoid';
        problem.options(1).trapezoid.nGrid = 20;
        
        problem.options(2).method = 'trapezoid';
        problem.options(2).trapezoid.nGrid = 40;
        
        problem.options(3).method = 'hermiteSimpson';
        problem.options(3).hermiteSimpson.nSegment = 20;
        
    case 'trapezoid'
        problem.options(1).method = 'trapezoid';
        problem.options(1).trapezoid.nGrid = 20;
        problem.options(2).method = 'trapezoid';
        problem.options(2).trapezoid.nGrid = 40;
        problem.options(3).method = 'trapezoid';
        problem.options(3).trapezoid.nGrid = 60;
        
    case 'rungeKutta'
        problem.options(1).method = 'rungeKutta';
        problem.options(1).defaultAccuracy = 'low';
        
        problem.options(2).method = 'rungeKutta';
        problem.options(2).defaultAccuracy = 'medium';
        
    case 'orthogonal'
        problem.options(1).method = 'chebyshev';
        problem.options(1).chebyshev.nColPts = 9;
        
        problem.options(2).method = 'chebyshev';
        problem.options(2).chebyshev.nColPts = 18;
    case 'gpops'
        problem.options(1).method = 'gpops';
        
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Start with a linear trajectory between four key frames:
% 0  --  initial configuration
% A  --  back swing
% B  --  front swing
% F  --  final configuration
%

tA = t0 + 0.25*(tF-t0);
xA = [-pi/2; 0];
dxA = [0;0];

tB = t0 + 0.75*(tF-t0);
xB = [pi/2; pi];
dxB = [0;0];

problem.guess.time = [t0, tA, tB, tF];
problem.guess.state = [[x0;dx0], [xA; dxA],[xB; dxB], [xF;dxF]];
problem.guess.control = [0, 0, 0, 0];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = optimTraj(problem);

% Interpolate the solution on a uniform grid for plotting and animation:
tGrid = soln(end).grid.time;
t = linspace(tGrid(1),tGrid(end),100);
z = soln(end).interp.state(t);
u = soln(end).interp.control(t);


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Plot the solution                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%HINT:  type help animate to figure out how to use the keyboard to interact
%with the animation (slow motion, pause, jump forward / backward...)

% Animate the results:
A.plotFunc = @(t,z)( drawAcrobot(t,z,dyn) );
A.speed = 0.25;
A.figNum = 101;
animate(t,z,A)

% Plot the results:
figure(1337); clf; plotAcrobot(t,z,u,dyn);

% Draw a stop-action animation:
figure(1338); clf; drawStopActionAcrobot(soln(end),dyn);


