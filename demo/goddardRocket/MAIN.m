% MAIN.m -- Goddard Rocket 
%
% This script runs a trajectory optimization to find the optimal thrust
% trajectory for the rocket to reach the maximum altitude. Physical
% parameters are roughly based on the SpaceX Falcon 9 rocket. 
%
% Dynamics include variable mass, inverse-square gravity, speed-dependent
% drag coefficient, height dependent air density.
%


%%%% Assumptions:
% SpaceX Falcon 9 rocket: 
% http://www.spacelaunchreport.com/falcon9v1-1.html
%
mFuel = 290000;  %(kg)  %mass of the fuel
mEmpty = 24000;  %(kg)  %mass of the rocket (without fuel)
Tmax = 390000;    %(N)   %Maximum thrust

error('Something is wrong with the dynamics model - it cannot accelerate off the launch pad')
%My guess is an error in the max thrust, with a divide by gravity error 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

h0 = 0;  %Rocket starts on the ground
v0 = 0;  %Rocket starts stationary
m0 = mFuel+mEmpty;  %Rocket starts full of fuel

vF = 0;  %Trying to reach maximum height
mF = mEmpty;  %Assume that we use all of the fuel

hLow = 0;   %Cannot go through the earth
hUpp = inf;  %To the moon!

vLow = 0; %Just look at the trajectory as it goes up
vUpp = inf;  % Go as fast as you can

mLow = mEmpty;
mUpp = mFuel+mEmpty;

uLow = 0;
uUpp = Tmax; %Maximum thrust output

P.bounds.initialTime.low = 0;
P.bounds.initialTime.upp = 0;

P.bounds.finalTime.low = 0;
P.bounds.finalTime.upp = inf;

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
P.func.bndObj = @(t0,x0,tF,xF)( xF(1) );  %Maximize final height


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Options and Method selection                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

P.options(1).method = 'trapazoid';
P.options(1).trapazoid.nGrid = 15;
P.options(1).nlpOpt = optimset(...
    'TolFun',1e-3,...
    'Display','iter',...
    'MaxFunEvals',1e4);

P.options(2).method = 'trapazoid';
P.options(2).trapazoid.nGrid = 25;
P.options(2).nlpOpt = optimset(...
    'TolFun',1e-6,...
    'Display','iter',...
    'MaxFunEvals',1e4);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Solve!                                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soln = trajOpt(P);


