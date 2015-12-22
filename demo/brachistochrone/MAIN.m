% MAIN.m
%
% Brachistochrone Problem - Optimal Control Formulation
%
% Derivation in Wolfram Mathworld:
% -->  Weisstein, Eric W. "Brachistochrone Problem." 
% -->  From MathWorld--A Wolfram Web Resource.
% -->  http://mathworld.wolfram.com/BrachistochroneProblem.html
%
% Problem Statement:
%   Find the shape of a wire between two points such that a bead sliding on
%   the wire traverses the distance in minimum time. 
%
%
% NOTES:
%   t = horizontal position = independent variable
%   x = vertical position = dependent variable
%   u = dx/dt = slope of the wire
%

clc; clear;

xDatum = 0.5;  % 1;
xFinal = -0.1;  % -0.5;

%%%% Objective Function
%
% Derived by Time = Integral { (differential path length)/(speed) }
problem.func.pathObj = @(t,x,u)( sqrt( (1+u.^2)./(xDatum-x) ) );


%%%% System Dynamics:
%
% Derived by definition: u = dx/dt = slope of wire
problem.func.dynamics = @(t,x,u)( u );


%%%% Boundary conditions:
%
%

% Initial horizontal position
problem.bounds.initialTime.low = 0;   
problem.bounds.initialTime.upp = 0;

% Final horizontal position
problem.bounds.finalTime.low = 1;
problem.bounds.finalTime.upp = 1;

% Initial vertical position
problem.bounds.initialState.low = 0;
problem.bounds.initialState.upp = 0;

% Final vertical position
problem.bounds.finalState.low = xFinal;
problem.bounds.finalState.upp = xFinal;

%
%
%%%% ~~~~~~~~~~~~~~~~~~~~


%%%% Initial Guess:
%
%  Straight line between boundary points
problem.guess.time = [0, 1];
problem.guess.state = [0, xFinal];
problem.guess.control = xFinal*[1,1];



%%%% Parameters:

% problem.options.method = 'chebyshev';
% problem.options.chebyshev.nColPts = 6;

problem.options.method = 'hermiteSimpson';
problem.options.hermiteSimpson.nSegment = 6;



%%%% Solve the problem!
%
soln = trajOpt(problem);


%%%% Plot the solution:
t = linspace(soln.grid.time(1), soln.grid.time(end), 100);
x = soln.interp.state(t);

plot(t,x);




