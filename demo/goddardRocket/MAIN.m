% MAIN.m -- Goddard Rocket 
%
% This script runs a trajectory optimization to find the optimal thrust
% trajectory for the rocket to reach the maximum altitude. Physical
% parameters are roughly based on the SpaceX Falcon 9 rocket. 
%
% Dynamics include variable mass, inverse-square gravity, speed-dependent
% drag coefficient, height dependent air density.
%

clc; clear;

%%%% Assumptions:
% SpaceX Falcon 9 rocket: 
% http://www.spacex.com/falcon9
%
mTotal = 505846;   %(kg)  %Total lift-off mass
mFuel = 0.8*mTotal;  %(kg)  %mass of the fuel
mEmpty = mTotal-mFuel;  %(kg)  %mass of the rocket (without fuel)
Tmax = 5885000;    %(N)   %Maximum thrust

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

h0 = 0;  %Rocket starts on the ground
v0 = 0;  %Rocket starts stationary
m0 = mTotal;  %Rocket starts full of fuel

vF = 0;  %Trying to reach maximum height
mF = mEmpty;  %Assume that we use all of the fuel

hLow = 0;   %Cannot go through the earth
hUpp = inf;  %To the moon!

vLow = 0; %Just look at the trajectory as it goes up
vUpp = inf;  % Go as fast as you can

mLow = mEmpty;
mUpp = mTotal;

uLow = 0;
uUpp = Tmax; %Maximum thrust output

P.bounds.initialTime.low = 0;
P.bounds.initialTime.upp = 0;

P.bounds.finalTime.low = 0;
P.bounds.finalTime.upp = 60*60;

P.bounds.state.low = [hLow;vLow;mLow];
P.bounds.state.upp = [hUpp;vUpp;mUpp];

P.bounds.initialState.low = [h0;v0;m0];
P.bounds.initialState.upp = [h0;v0;m0];

P.bounds.finalState.low = [hLow;vF;mF];
P.bounds.finalState.upp = [hUpp;vF;mF];

P.bounds.control.low = uLow;
P.bounds.control.upp = uUpp;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Initial Guess                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
hGuess = 2e4;   %(m) guess at the maximum height reached
P.guess.time = [0, 180];  %(s)
P.guess.state = [ [h0;v0;m0],  [hGuess;vF;mF] ];
P.guess.control = [uUpp, uLow];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Objective and Dynamic functions                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Dynamics function:
P.func.dynamics = @(t,x,u)( rocketDynamics(x,u) );

% Objective function:
P.func.bndObj = @(t0,x0,tF,xF)( -xF(1)/10000 );  %Maximize final height


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Options and Method selection                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

P.options(1).method = 'trapazoid';
P.options(1).trapazoid.nGrid = 10;
P.options(1).nlpOpt = optimset(...
    'TolFun',1e-3,...
    'Display','iter',...
    'MaxFunEvals',1e4);

P.options(2).method = 'trapazoid';
P.options(2).trapazoid.nGrid = 25;
P.options(2).nlpOpt = optimset(...
    'TolFun',1e-4,...
    'Display','iter',...
    'MaxFunEvals',1e4);

P.options(3).method = 'trapazoid';
P.options(3).trapazoid.nGrid = 35;
P.options(3).nlpOpt = optimset(...
    'TolFun',1e-6,...
    'Display','iter',...
    'MaxFunEvals',1e4);


%%%% NOTES:
% 
% 1) Chebyshev is not a good method for this problem, beause their is a
% discontinuity in solution of the thrust curve. This will cause ringing in
% the trajectory: an artifact, rather than a property of the solution. For
% this reason, a non-global method, such as trapazoid is desirable.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Solve!                                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soln = trajOpt(P);

t = linspace(soln(end).grid.time(1),soln(end).grid.time(end),250);
x = soln(end).interp.state(t);
u = soln(end).interp.control(t);

figure(120);
subplot(2,2,1);
plot(t,x(1,:)/1000)
xlabel('time (s)')
ylabel('height (km)')
title('Maximal Height Trajectory')
subplot(2,2,2);
plot(t,x(3,:))
xlabel('time (s)')
ylabel('mass (kg)')
title('Goddard Rocket')
subplot(2,2,3);
plot(t,x(2,:))
xlabel('time (s)')
ylabel('velocity (m/s)')
subplot(2,2,4);
plot(t,u/1000)
xlabel('time (s)')
ylabel('thrust (kN)')
