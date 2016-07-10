% MAIN.m -- Goddard Rocket
%
% This script runs a trajectory optimization to find the optimal thrust
% trajectory for the rocket to reach the maximum altitude. Physical
% parameters are roughly based on the SpaceX Falcon 9 rocket.
%
% Dynamics include variable mass, inverse-square gravity, speed-dependent
% drag coefficient, height dependent air density.
%
% NOTES:
%   This problem sort of converges, but not very well. I think that there
%   is a singular arc in it that is not being handled correctly. It is
%   still interesting to see as an example of ways in which problems might
%   misbehave.
%

clc; clear;
addpath ../../

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

method = 'trapezoid';
% method = 'rungeKutta';
% method = 'chebyshev';

switch method
    
    case 'trapezoid'
        
        P.options(1).method = 'trapezoid';
        P.options(1).defaultAccuracy = 'low';
        
        P.options(2).method = 'trapezoid';
        P.options(2).defaultAccuracy = 'medium';
        P.options(2).nlpOpt.MaxFunEvals = 2e4;
        P.options(2).nlpOpt.MaxIter = 1e5;
        
    case 'rungeKutta'
        P.options(1).method = 'rungeKutta';
        P.options(1).defaultAccuracy = 'low';
        
        P.options(2).method = 'rungeKutta';
        P.options(2).defaultAccuracy = 'medium';
        
    case 'chebyshev'
        
        P.options(1).method = 'chebyshev';
        P.options(1).defaultAccuracy = 'low';
        
        P.options(2).method = 'chebyshev';
        P.options(2).defaultAccuracy = 'low';
        P.options(2).chebyshev.nColPts = 15;
        
end


%%%% NOTES:
%
% 1) Orthogonal collocation (chebyshev) is not a good method for this problem, beause there is a
% discontinuity in solution of the thrust curve. It still sort of works,
% but will find a sub-optimal answer, or produce ringing.
%
% 2) Why does the 'trapezoid' low resolution version finish so quickly and the medium
% quality one take forever? Hint: Look at the feasibility printout: it is
% cyclical. If you were to plot the solution as a function of iteration,
% you would find that occasionally the discontinuity moves, which causes a
% consistency error in the NLP. Eventually it gets to the "right" answer,
% although it is pretty boring. I suspect that you could get more
% interesting behavior with different constants.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Solve!                                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soln = optimTraj(P);

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
