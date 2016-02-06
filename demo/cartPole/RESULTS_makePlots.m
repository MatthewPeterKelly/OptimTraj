% MAIN.m
%
% Solve the cart-pole swing-up problem

clc; clear;

p.m1 = 1.0;  % (kg) Cart mass
p.m2 = 0.3;  % (kg) pole mass
p.g = 9.81;  % (m/s^2) gravity
p.l = 0.5;   % (m) pendulum (pole) length

dist = 1.0;  %How far must the cart translate during its swing-up
maxForce = 20;  %Maximum actuator forces    
duration = 2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( cartPoleDynamics(x,u,p) );
problem.func.pathObj = @(t,x,u)( u.^2 );  %Force-squared cost function

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up problem bounds                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = duration;
problem.bounds.finalTime.upp = duration;

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

problem.guess.time = [0,duration];
problem.guess.state = [problem.bounds.initialState.low, problem.bounds.finalState.low];
problem.guess.control = [0,0];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e5);

% problem.options.method = 'trapezoid'; problem.options.trapezoid.nGrid = 20;
problem.options.method = 'hermiteSimpson'; problem.options.hermiteSimpson.nSegment = 25;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = trajOpt(problem);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Display Solution                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Unpack the simulation
n = length(soln.grid.time);
t = linspace(soln.grid.time(1), soln.grid.time(end), 15*(n-1)+1);
z = soln.interp.state(t);
u = soln.interp.control(t);

%%%% Plots:
figure(1); clf;
plotPendulumCart(t,z,u,p);


%%%% Draw Trajectory:
[p1,p2] = cartPoleKinematics(z,p);

figure(2); clf;
nFrame = 9;  %Number of frames to draw
drawCartPoleTraj(t,p1,p2,nFrame);


%%%% Save an animation:
% % val = [p1,p2];
% % xLow = min(val(1,:));
% % xUpp = max(val(1,:));
% % yLow = min(val(2,:));
% % yUpp = max(val(2,:));
% % drawFun = @(t,p)( drawCartPoleAnim(t,p,xLow, xUpp, yLow, yUpp) );
% % P.plotFunc = drawFun;
% % P.figNum = 7;
% % P.frameRate = 24;
% % P.fileName = 'cartPoleAnimation';
% % saveAnimation(t,[p1;p2],P);


%%%% Show only solution grid:
figure(3); clf;
tGrid = soln.grid.time;
zGrid = soln.grid.state;
uGrid = soln.grid.control;

subplot(3,1,1);
plot(tGrid,zGrid(1,:),'ko')

subplot(3,1,2);
plot(tGrid,zGrid(2,:),'ko')

subplot(3,1,3);
plot(tGrid,uGrid,'ko')

%%%% Show both grids:
figure(4); clf;
tGrid = soln.grid.time;
zGrid = soln.grid.state;
uGrid = soln.grid.control;

idx = 1:2:length(tGrid);  %Only plot knot points

colorState = [0.2,0.2,0.8];
colorControl = [0.6, 0.1, 0.7];

subplot(3,1,1); hold on
plot(t,z(1,:),'Color',colorState,'LineWidth',3)
plot(tGrid(idx),zGrid(1,idx),'ko','MarkerSize',10,'LineWidth',2)

subplot(3,1,2); hold on
plot(t,z(2,:),'Color',colorState,'LineWidth',3)
plot(tGrid(idx),zGrid(2,idx),'ko','MarkerSize',10,'LineWidth',2)

subplot(3,1,3); hold on
plot(t,u,'Color',colorControl,'LineWidth',3)
plot(tGrid(idx),uGrid(idx),'ko','MarkerSize',10,'LineWidth',2)


%%%% Show the error in the collocation constraint between grid points:
figure(5); clf;

idx = 1:2:length(tGrid);  %Only plot knot points
cc = soln.interp.collCst(t);
ccIdx = soln.interp.collCst(tGrid(idx));

subplot(2,2,1); hold on;
plot(tGrid(idx),ccIdx(1,:),'ko','MarkerSize',7,'LineWidth',2);
plot(t,cc(1,:))
title('Collocation Error:   dx/dt - f(t,x,u)');
ylabel('d/dt cart position');

subplot(2,2,3); hold on;
plot(tGrid(idx),ccIdx(2,:),'ko','MarkerSize',7,'LineWidth',2);
plot(t,cc(2,:))
xlabel('time')
ylabel('d/dt pole angle')

idx = 1:length(soln.info.error);
subplot(2,2,2); hold on;
plot(idx,soln.info.error(1,:),'ko','MarkerSize',8,'LineWidth',3);
title('State Error')
ylabel('cart position')

subplot(2,2,4); hold on;
plot(idx,soln.info.error(2,:),'ko','MarkerSize',8,'LineWidth',3);
xlabel('segment index')
ylabel('pole angle');


%%%% Save script for paper:
% save2pdf('cartPole_drawSoln_25.pdf',figure(2));
% save2pdf('cartPole_plotSoln_25.pdf',figure(4));
% save2pdf('cartPole_error_25.pdf',figure(5));
