% MAIN.m  --  Five Link Biped trajectory optimization --
%
% This script sets up and then solves the optimal trajectory for the five
% link biped, assuming that the walking gait is compused of single-stance
% phases of motion connected by impulsive heel-strike (no double-stance or
% flight phases).
%
% Optimize for minimum cost of transport. This code is far more complicated
% to understand than the torque-squared problem, and some aspects of the
% indexing are not as well documented. For example, to get
% torque-rate-squared regularization, the torque is actually included
% inside of the state vector. Additionally, the abs(power) cost function is
% computed using slack variables to prevent discontinuous a discontinuity
% in the objective function.
%
% The equations of motion and gradients are all derived by:
%   --> Derive_Equations.m
%

%%%% NOTE %%%%
%
% This example - at least for the cost of transport optimization - should
% be considered experimental. This code does not pass strict convergence
% tests - The optimization completes successfully with loose tolerances,
% but fails to converge to a unique solution with more tight tolerances.
%
% 

clc; clear;
addpath ../../../

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up parameters and options                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
param = getPhysicalParameters();

param.stepLength = 0.5;
param.stepTime = 0.7;
param.stepHeight = 0.001;  %Foot must clear this height at mid-stance

param.gammaNeg = 1;   %Cost for negative work
param.gammaPos = 1;  %Cost for positive work
param.alpha = 0;   %Torque-squared smoothing parameter;
param.beta = 1e-3;   %TorqueRate-squared smoothing parameter;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics =  @(t,x,u)( dynamics(t,x,u,param) );

problem.func.pathObj = @(t,x,u)( obj_costOfTransport(x,u,param) );

problem.func.bndCst = @(t0,x0,tF,xF)( stepConstraint(x0,xF,param) );

problem.func.pathCst = @(t,x,u)( pathConstraint(x,u,param) );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = param.stepTime;
problem.bounds.finalTime.upp = param.stepTime;

% State: (absolute reference frames)
%   1 = stance leg tibia angle
%   2 = stance leg femur angle
%   3 = torso angle
%   4 = swing leg femur angle
%   5 = swing leg tibia angle

qLow = (-pi/2)*ones(5,1);
qUpp = (pi/2)*ones(5,1);
dqLow = -10*ones(5,1);
dqUpp = 10*ones(5,1);
uMax = 100;  %Nm
uLow = -uMax*ones(5,1);
uUpp = uMax*ones(5,1);
problem.bounds.state.low = [qLow; dqLow; uLow];
problem.bounds.state.upp = [qUpp; dqUpp; uUpp];
problem.bounds.initialState.low = [qLow; dqLow; uLow];
problem.bounds.initialState.upp = [qUpp; dqUpp; uUpp];
problem.bounds.finalState.low = [qLow; dqLow; uLow];
problem.bounds.finalState.upp = [qUpp; dqUpp; uUpp];

problem.bounds.control.low = [-inf(5,1); zeros(10,1)];   % [torque rate; slack]
problem.bounds.control.upp = [inf(5,1); inf(10,1)];

% Disable the stance ankle motor:
problem.bounds.state.low(5+5+1) = 0;
problem.bounds.state.upp(5+5+1) = 0;
problem.bounds.initialState.low(5+5+1) = 0;
problem.bounds.initialState.upp(5+5+1) = 0;
problem.bounds.finalState.low(5+5+1) = 0;
problem.bounds.finalState.upp(5+5+1) = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% For now, just assume a linear trajectory between boundary values

problem.guess.time = [0, param.stepTime];

q0 = [...
    -0.3; % stance leg tibia angle
    0.7; % stance leg femur angle
    0.0; % torso angle
    -0.5; % swing leg femur angle
    -0.6]; % swing leg tibia angle
qF = q0([5;4;3;2;1]);   %Flip left-right

dq0 = (qF-q0)/param.stepTime;
dqF = dq0;

u0 = zeros(5,1); uF = zeros(5,1); %Start with passive trajectory

problem.guess.state = [q0, qF; dq0, dqF; u0, uF];
problem.guess.control = zeros(5+10,2);  


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%NOTE:  Here I choose to run the optimization twice, mostly to demonstrate
%   functionality, although this can be important on harder problems. I've
%   explicitly written out many options below, but the solver will fill in
%   almost all defaults for you if they are ommitted.

% method = 'test1';
% method = 'test4';
method = 'trapezoid';
% method = 'hermiteSimpson';

switch method
    
    case 'test1'
        
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 20;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        problem.options(1).nlpOpt.MaxIter = 500;
        problem.options(1).nlpOpt.TolFun = 1e-3;
%         problem.options(1).nlpOpt.TolX = 1e-6;
        
    case 'test4'
        
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        problem.options(1).nlpOpt.MaxIter = 500;
        problem.options(1).nlpOpt.TolFun = 1e-3;
%         problem.options(1).nlpOpt.TolX = 1e-6;

    case 'trapezoid'
        
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 15;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        problem.options(1).nlpOpt.MaxIter = 1e3;
        problem.options(1).nlpOpt.TolFun = 1e-4;
        
        problem.options(2).method = 'trapezoid'; % Select the transcription method
        problem.options(2).trapezoid.nGrid = 30;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        problem.options(2).nlpOpt.MaxIter = 1e4;
        problem.options(1).nlpOpt.TolFun = 1e-4;
        
    case 'hermiteSimpson'
        
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = 6;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(2).hermiteSimpson.nSegment = 15;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
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
q = soln(end).grid.state(1:5,:);
dq = soln(end).grid.state(6:10,:);
u = soln(end).grid.state(11:15,:);
du = soln(end).grid.control(1:5,:);
sn = soln(end).grid.control(6:10,:);   %Slack variable for negative power
sp = soln(end).grid.control(11:15,:);   % Slack variable for positive power

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Plot the solution                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

Anim.figNum = 1;
Anim.speed = 0.25;
Anim.plotFunc = @(t,q)( drawRobot(q,param) );
Anim.verbose = true;
animate(t,q,Anim);

figure(2); clf;
subplot(1,2,1);
plot(t,q);
legend('q1','q2','q3','q4','q5');
xlabel('time')
ylabel('link angles')
subplot(1,2,2);
plot(t,u);
legend('u1','u2','u3','u4','u5');
xlabel('time')
ylabel('joint torques')





