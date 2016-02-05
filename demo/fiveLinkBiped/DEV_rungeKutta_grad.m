% DEV_rungeKutta_grad.m
%
% Script for testing and development of the runge-kutta option for TrajOpt
% with analytic gradients.
%
% WARNING:  This script uses adaptive numerical gradients to check the
% analytic gradients. These are SLOW, causing this script to take a long
% time to run. 
% 

clc; clear; 
addpath ../../

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up parameters and options                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
param = getPhysicalParameters();

param.stepLength = 0.5;
param.stepTime = 0.7;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics =  @(t,x,u)( dynamics(t,x,u,param) );

problem.func.pathObj = @(t,x,u)( obj_torqueSquared(u) );

problem.func.bndCst = @(t0,x0,tF,xF)( stepConstraint(x0,xF,param) );

problem.func.pathCst = @(t,x,u)( DEV_pathConstraint(t,x) );


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
problem.bounds.control.low = -uMax*ones(5,1);
problem.bounds.control.upp = uMax*ones(5,1);

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

problem.guess.control = zeros(5,2);  %Start with passive trajectory


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                   Development / Debugging Scripts                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%NOTE:  Here I choose to run the optimization twice, mostly to demonstrate
%   functionality, although this can be important on harder problems. I've
%   explicitly written out many options below, but the solver will fill in
%   almost all defaults for you if they are ommitted.

%NOTE:  There are many methods written out below, just to show different
%   ways to solve the problem.

% method = 'Baseline Solution';
method = 'rungeKutta_test';

switch method
    
    case 'Baseline Solution'  %hermite simpson with analytic gradients
        
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = 6;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'on';
        
        problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(2).hermiteSimpson.nSegment = 15;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';

    
    case 'rungeKutta_test'
        
        problem.options(1).method = 'rungeKutta'; 
        problem.options(1).defaultAccuracy = 'low';
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'on';
        
        problem.options(2).method = 'rungeKutta'; 
        problem.options(2).defaultAccuracy = 'low';
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        problem.options(2).nlpOpt.DerivativeCheck = 'on';
        
        % The following option takes a long time, but verifies accuracy of
        % gradients more robustly than fmincon's internal checks.
        problem.options(2).rungeKutta.AdaptiveDerivativeCheck = 'on';
        
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
u = soln(end).grid.control;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Plot the solution                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Anim.figNum = 1; clf(Anim.figNum);
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





