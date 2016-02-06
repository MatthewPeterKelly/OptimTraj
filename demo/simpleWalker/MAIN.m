%MAIN.m  --  simple walker trajectory optimization
%
% This script sets up a trajectory optimization problem for a simple model
% of walking, and solves it using TrajOpt. The walking model is a double
% pendulum, with point feet, no ankle torques, impulsive heel-strike (but
% not push-off), and continuous hip torque. Both legs have inertia. Cost
% function is minimize integral of torque-squared.
%
%

clc; clear;
addpath ../../

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Parameters for the dynamics function                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
param.dyn.m = 10;  %leg mass
param.dyn.I = 1;  %leg inertia about CoM
param.dyn.g = 9.81;  %gravity
param.dyn.l = 1;  %leg length
param.dyn.d = 0.2;  %Leg CoM distance from hip

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( dynamics(x,u,param.dyn) );

problem.func.pathObj = @(t,x,u)( costFun(u) );

problem.func.bndCst = @(t0,x0,tF,xF)( periodicGait(xF,x0,param.dyn) );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
t0 = 0;  tF = 1;
problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0;
problem.bounds.finalTime.low = tF;
problem.bounds.finalTime.upp = tF;

% State: [q1;q2;dq1;dq2];

problem.bounds.state.low = [-pi/3; -pi/3; -inf(2,1)];
problem.bounds.state.upp = [ pi/3;  pi/3;  inf(2,1)];

stepAngle = 0.2;
problem.bounds.initialState.low = [stepAngle; -stepAngle; -inf(2,1)];
problem.bounds.initialState.upp = [stepAngle; -stepAngle;  inf(2,1)];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% For now, just assume a linear trajectory between boundary values

problem.guess.time = [t0, tF];

stepRate = (2*stepAngle)/(tF-t0);
x0 = [stepAngle; -stepAngle; -stepRate; stepRate];
xF = [-stepAngle; stepAngle; -stepRate; stepRate];
problem.guess.state = [x0, xF];

problem.guess.control = [0, 0];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%NOTE:  Here I choose to run the optimization twice, mostly to demonstrate
%   functionality, although this can be important on harder problems. I've
%   explicitly written out many options below, but the solver will fill in
%   almost all defaults for you if they are ommitted.

% method = 'trapezoid';
method = 'hermiteSimpson';
% method = 'chebyshev';
% method = 'rungeKutta';
% method = 'gpops';

switch method
    case 'trapezoid'
        
        % First iteration: get a more reasonable guess
        problem.options(1).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-3,...
            'MaxFunEvals',1e4);   %options for fmincon
        problem.options(1).verbose = 3; % How much to print out?
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options
        
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-6,...
            'MaxFunEvals',5e4);   %options for fmincon
        problem.options(2).verbose = 3; % How much to print out?
        problem.options(2).method = 'trapezoid'; % Select the transcription method
        problem.options(2).trapezoid.nGrid = 25;  %method-specific options
        
    case 'hermiteSimpson'
        
        % First iteration: get a more reasonable guess
        problem.options(1).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-3,...
            'MaxFunEvals',1e4);   %options for fmincon
        problem.options(1).verbose = 3; % How much to print out?
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = 6;  %method-specific options
        
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-6,...
            'MaxFunEvals',5e4);   %options for fmincon
        problem.options(2).verbose = 3; % How much to print out?
        problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(2).hermiteSimpson.nSegment = 15;  %method-specific options
        
        
    case 'chebyshev'
        
        % First iteration: get a more reasonable guess
        problem.options(1).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-3,...
            'MaxFunEvals',1e4);   %options for fmincon
        problem.options(1).verbose = 3; % How much to print out?
        problem.options(1).method = 'chebyshev'; % Select the transcription method
        problem.options(1).chebyshev.nColPts = 9;  %method-specific options
        
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-8,...
            'MaxFunEvals',5e4);   %options for fmincon
        problem.options(2).verbose = 3; % How much to print out?
        problem.options(2).method = 'chebyshev'; % Select the transcription method
        problem.options(2).chebyshev.nColPts = 15;  %method-specific options
     
    case 'rungeKutta'
        problem.options(1).method = 'rungeKutta'; % Select the transcription method
        problem.options(1).defaultAccuracy = 'low';
        problem.options(2).method = 'rungeKutta'; % Select the transcription method
        problem.options(2).defaultAccuracy = 'medium';
        
    case 'gpops'
        problem.options.method = 'gpops';
        problem.options.defaultAccuracy = 'medium';
        
    otherwise
        error('Invalid method!');
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%%% THE KEY LINE:
soln = trajOpt(problem);

% Transcription Grid points:
t = soln(end).grid.time;
q1 = soln(end).grid.state(1,:);
q2 = soln(end).grid.state(2,:);
dq1 = soln(end).grid.state(3,:);
dq2 = soln(end).grid.state(4,:);
u = soln(end).grid.control;

% Interpolated solution:
tInt = linspace(t(1),t(end),10*length(t)+1);
xInt = soln(end).interp.state(tInt);
q1Int = xInt(1,:);
q2Int = xInt(2,:);
dq1Int = xInt(3,:);
dq2Int = xInt(4,:);
uInt = soln(end).interp.control(tInt);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Plot the solution                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

figure(100); clf;

subplot(3,1,1); hold on;
plot(tInt,q1Int,'r-'); plot(tInt,q2Int,'b-');
plot([t(1),t(end)],[0,0],'k--','LineWidth',1);
plot(t,q1,'ro'); plot(t,q2,'bo');
legend('leg one','leg two')
xlabel('time (sec)')
ylabel('angle (rad)')
title('Leg Angles')

subplot(3,1,2); hold on;
plot(tInt,dq1Int,'r-'); plot(tInt,dq2Int,'b-');
plot(t,dq1,'ro'); plot(t,dq2,'bo');
legend('leg one','leg two')
xlabel('time (sec)')
ylabel('rate (rad/sec)')
title('Leg Angle Rates')

subplot(3,1,3); hold on;
plot(t,u,'mo'); plot(tInt,uInt,'m-');
xlabel('time (sec)')
ylabel('torque (Nm)')
title('Hip Torque')



