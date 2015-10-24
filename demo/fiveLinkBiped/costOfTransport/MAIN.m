% MAIN.m  --  Five Link Biped trajectory optimization --
%
% This script sets up and then solves the optimal trajectory for the five
% link biped, assuming that the walking gait is compused of single-stance
% phases of motion connected by impulsive heel-strike (no double-stance or
% flight phases).
%
% Optimize for minimum cost of transport
%
% The equations of motion and gradients are all derived by:
%   --> Derive_Equations.m
%

clc; clear;
addpath ../../../

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up parameters and options                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
param = getPhysicalParameters();

param.stepLength = 0.5;
param.stepTime = 0.7;

param.alpha = 1e-1;   %Torque-squared smoothing parameter;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics =  @(t,x,u)( dynamics(t,x,u,param) );

problem.func.pathObj = @(t,x,u)( obj_costOfTransport(u,param) );

problem.func.bndCst = @(t0,x0,tF,xF)( stepConstraint(x0,xF,param) );

problem.func.pathCst = @(t,x,u)( pathConstraint(x,u) );


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
problem.bounds.state.low = [qLow; dqLow];
problem.bounds.state.upp = [qUpp; dqUpp];
problem.bounds.initialstate.low = [qLow; dqLow];
problem.bounds.initialstate.upp = [qUpp; dqUpp];
problem.bounds.finalstate.low = [qLow; dqLow];
problem.bounds.finalstate.upp = [qUpp; dqUpp];

uMax = 100;  %Nm
problem.bounds.control.low = [-uMax*ones(5,1); zeros(10,1)];   % [control; slack]
problem.bounds.control.upp = [uMax*ones(5,1); inf(10,1)];

% Disable the stance ankle motor:
problem.bounds.control.low(1) = 0;
problem.bounds.control.upp(1) = 0;


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

problem.guess.state = [q0, qF; dq0, dqF];

problem.guess.control = zeros(5+10,2);  %Start with passive trajectory


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%NOTE:  Here I choose to run the optimization twice, mostly to demonstrate
%   functionality, although this can be important on harder problems. I've
%   explicitly written out many options below, but the solver will fill in
%   almost all defaults for you if they are ommitted.

% method = 'test';
method = 'trapazoid';
% method = 'hermiteSimpson';

switch method
    
    case 'test'
        
        problem.options(1).method = 'trapazoid'; % Select the transcription method
        problem.options(1).trapazoid.nGrid = 20;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        problem.options(1).nlpOpt.MaxIter = 5000;
        
    case 'trapazoid'
        
        problem.options(1).method = 'trapazoid'; % Select the transcription method
        problem.options(1).trapazoid.nGrid = 15;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'trapazoid'; % Select the transcription method
        problem.options(2).trapazoid.nGrid = 30;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
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
u = soln(end).grid.control(1:5,:);
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





