%TEST  --  Dynamics
%
% This script runs a simple open-loop simulation of the robot to do a
% sanity check on the dynamics, energy, and contact force equtions.
%


clc; clear;

% Loads the struct of physical parameters (masses, lengths, ...)
p = getPhysicalParameters();

% Initial conditions:
q0 = [-0.25; 0.2; -0.1; -0.35; -0.55];
dq0 = [0.2; -0.15; 0.05; 0.4; 0.45];

z0 = [q0;dq0];

tSpan = [0,0.5];
dynFun = @(t,z)( dynamics(z,zeros(5,1),p) );

options = odeset(...
    'AbsTol',1e-10,...
    'RelTol',1e-10);

% Simulation with ode45:
sol = ode45(dynFun,tSpan,z0,options);
t = linspace(tSpan(1),tSpan(2),150);
z = deval(sol,t);
q = z(1:5,:);
dq = z(6:10,:);
u = zeros(5,length(t));
ddq = dynSs(q,dq,u,p);
[Fx, Fy] = contactForces(q,dq,ddq,p);
[KE, PE] = energy(q,dq,p);

% Animate the result:
Anim.figNum = 1;
Anim.plotFunc = @(t,q)( drawRobot(q,p) );
Anim.verbose = true;
Anim.speed = 0.25;  %Playback at half-speed
animate(t,q,Anim);

% Check energy conservation:
sprintf('Energy variation: %4.4g',std(KE+PE))


