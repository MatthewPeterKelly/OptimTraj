% MAIN.m
%
% Solve the cart-pole swing-up problem  --  minimum time
%
% Note:  This problem is much more difficult to solve than the
% minimum-force version. This is because most of the control trajectory is
% sitting on a constraint: the maximum or minimum control force. This is
% generally true of minimum-time trajectories: they have bang-bang
% solutions. To get the exact solution, you would need to do many steps of
% mesh refinement. Here I only do two iterations, to keep total time
% reasonable. Another problem with minimum-time objective functions is that
% they sometimes have singular arcs: solutions where there is no single
% best control trajectory. This will manifest itself as "chattering" in the
% control trajectory and slow convergence. One solution is to include a
% regularization term, such as force squared with a very small coefficient,
% which forces a unique solution along the singular arc.
%

clc; clear;
addpath ../../

p.m1 = 2.0;  % (kg) Cart mass
p.m2 = 0.5;  % (kg) pole mass
p.g = 9.81;  % (m/s^2) gravity
p.l = 0.5;   % (m) pendulum (pole) length

dist = 1.0;  %How far must the cart translate during its swing-up
maxForce = 50;  %Maximum actuator forces



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( cartPoleDynamics(x,u,p) );
problem.func.pathObj = @(t,x,u)( ones(size(t)) ); 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up problem bounds                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 0.01;
problem.bounds.finalTime.upp = inf;

problem.bounds.initialState.low = zeros(4,1);
problem.bounds.initialState.upp = zeros(4,1);
problem.bounds.finalState.low = [dist;pi;0;0];
problem.bounds.finalState.upp = [dist;pi;0;0];

problem.bounds.state.low = [-2*dist;-2*pi;-inf;-inf];
problem.bounds.state.upp = [2*dist;2*pi;inf;inf];

problem.bounds.control.low = -maxForce;
problem.bounds.control.upp = maxForce;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Initial guess at trajectory                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.guess.time = [0,2];
problem.guess.state = [problem.bounds.initialState.low, problem.bounds.finalState.low];
problem.guess.control = [0,0];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options(1).nlpOpt = optimset(...
    'Display','iter',...
    'TolFun',1e-3,...
    'MaxFunEvals',1e5);
problem.options(1).method = 'trapezoid';
problem.options(1).trapezoid.nGrid = 10;

problem.options(2).nlpOpt = optimset(...
    'Display','iter',...
    'TolFun',1e-6,...
    'MaxFunEvals',1e5);
problem.options(2).method = 'trapezoid';
problem.options(2).trapezoid.nGrid = 30;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = optimTraj(problem);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Display Solution                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Unpack the simulation
t = linspace(soln(end).grid.time(1), soln(end).grid.time(end), 150);
z = soln(end).interp.state(t);
u = soln(end).interp.control(t);

%%%% Plots:

%%%% Draw Trajectory:
[p1,p2] = cartPoleKinematics(z,p);

figure(2); clf;
nFrame = 9;  %Number of frames to draw
drawCartPoleTraj(t,p1,p2,nFrame);


%%%% Show the error in the collocation constraint between grid points:
%
if strcmp(soln(end).problem.options.method,'trapezoid') || strcmp(soln(end).problem.options.method,'hermiteSimpson')
    % Then we can plot an estimate of the error along the trajectory
    figure(5); clf;
    
    % NOTE: the following commands have only been implemented for the direct
    % collocation(trapezoid, hermiteSimpson) methods, and will not work for
    % chebyshev or rungeKutta methods.
    cc = soln(end).interp.collCst(t);
    
    subplot(2,2,1);
    plot(t,cc(1,:))
    title('Collocation Error:   dx/dt - f(t,x,u)')
    ylabel('d/dt cart position')
    
    subplot(2,2,3);
    plot(t,cc(2,:))
    xlabel('time')
    ylabel('d/dt pole angle')
    
    idx = 1:length(soln(end).info.error);
    subplot(2,2,2); hold on;
    plot(idx,soln(end).info.error(1,:),'ko');
    title('State Error')
    ylabel('cart position')
    
    subplot(2,2,4); hold on;
    plot(idx,soln(end).info.error(2,:),'ko');
    xlabel('segment index')
    ylabel('pole angle');
end

%%%% Plot the state and control against time
figure(1); clf;
plotPendulumCart(t,z,u,p);






