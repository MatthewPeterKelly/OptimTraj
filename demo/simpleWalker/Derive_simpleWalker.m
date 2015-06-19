% Derive_simpleWalker.m
%
% This script uses Matlab symbolic toolbox to derive the dynamics and
% kinematics equations for the simple walker model, which is mathematically
% identical to an acrobot.
%
% NOTATION:
% 
%   1 = stance leg (connected to the ground)
%   2 = swing leg (hanging from the hip)
%
%   q = angle
%   dq = dq/dt = angular rate
%   ddq = ddq/ddt = angular acceleration
%
clc; clear;

syms q1 q2 dq1 dq2 ddq1 ddq2 'real'   % states 
syms u 'real' % hip torque
syms d m I g l 'real' % physical parameters

% d = distance along leg from hip to the center of mass of the leg
% m = mass of each leg
% I = moment of inertia of each leg about its center of mass
% g = gravity
% l = leg length

%%%% Unit vectors:
i = sym([1;0]);
j = sym([0;1]);

e1 = cos(q1)*(-j) + sin(q1)*(i);    % hip -> stance foot
e2 = cos(q2)*(-j) + sin(q2)*(i);    % hip -> swing foot

%%%% State vectors:
z = [q1;q2;dq1;dq2];
dz = [dq1;dq2;ddq1;ddq2];

%%%% Kinematics:
pHip = -l*e1;
p1 = pHip +d*e1;   %Center of mass of leg one
p2 = pHip +d*e2;   %Center of mass of leg two

dp1 = jacobian(p1,z)*dz;  %Chain rule to get velocity of hip joint
dp2 = jacobian(p2,z)*dz; 

ddp1 = jacobian(dp1,z)*dz;  
ddp2 = jacobian(dp2,z)*dz; 

%%%% Define a function for doing '2d' cross product: dot(a x b, k)
cross2d = @(a,b)(a(1)*b(2) - a(2)*b(1));

%%%% Angular momentum balance of system about stance foot (origin)
sumTorques1 = cross2d(p1,-m*g*j) + cross2d(p2,-m*g*j);
sumInertial1 = cross2d(p1,m*ddp1) + I*ddq1 + cross2d(p2,m*ddp2) + I*ddq2;
eqn1 = sumTorques1-sumInertial1;

%%%% Angular momentum balance of swing leg about hip joint:
sumTorques2 = cross2d(p2-pHip,-m*g*j) + u;
sumInertial2 = cross2d(p2-pHip,m*ddp2) + I*ddq2;
eqn2 = sumTorques2-sumInertial2;

%%%% Solve dynamics:
ddq = [ddq1;ddq2];
eqns = [eqn1;eqn2];
[MM,ff] = equationsToMatrix(eqns,ddq);
soln.ddq = MM\ff;
soln.ddq1 = simplify(soln.ddq(1));
soln.ddq2 = simplify(soln.ddq(2));

%%%% Generate an optimized matlab function for dynamics:
matlabFunction(soln.ddq1,soln.ddq2,...
    'file','autoGen_dynamics.m',...
    'vars',{q1,q2,dq1,dq2,u,d, m, I, g, l},...
    'outputs',{'ddq1','ddq2'});

%%%% Compute the energy of the system:
U = m*g*dot(p1,j) + m*g*dot(p2,j);   %Potential Energy
T = 0.5*m*dot(dp1,dp1) + 0.5*m*dot(dp2,dp2) + 0.5*I*dq1^2 + 0.5*I*dq2^2;   %Kinetic Energy

%%%% Generate an optimized matlab function for energy:
matlabFunction(U,T,...
    'file','autoGen_energy.m',...
    'vars',{q1,q2,dq1,dq2,d, m, I, g, l},...
    'outputs',{'U','T'});

%%%% Generate a function for computing the kinematics:
matlabFunction(p1,p2,dp1,dp2,...
    'file','autoGen_kinematics.m',...
    'vars',{q1,q2,dq1,dq2,d,l},...
    'outputs',{'p1','p2','dp1','dp2'});



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%         Derive heel-strike map and collision mechanics                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

pFoot = pHip + l*e2;  %Swing foot position

% Angular momentum of the system about the new stance foot (old swing foot)
hSysBefore = ...
    cross2d(p1-pFoot,m*dp1) + I*dq1 + ...   % old stance leg
    cross2d(p2-pFoot,m*dp2) + I*dq2;    % old swing leg

% Angular momentum of the old stance leg about the hip
hLegBefore = cross2d(p1-pHip,m*dp1) + I*dq1;  % old stance leg

% Introduce new variables for the state after the collision:
q1New = q2;
q2New = q1;
syms dq1New dq2New  'real'   % angular rates after collision

% Unit vectors after the collision:   (new naming convention)
e1New = cos(q1New)*(-j) + sin(q1New)*(i);    % hip -> stance foot
e2New = cos(q2New)*(-j) + sin(q2New)*(i);    % hip -> swing foot

% Kinematics:
pHipNew = -l*e1New;
p1New = pHipNew + d*e1New;
p2New = pHipNew + d*e2New;

dp1New = jacobian(p1New,[q1New;q2New])*[dq1New;dq2New];  
dp2New = jacobian(p2New,[q1New;q2New])*[dq1New;dq2New];  

% Angular momentum of the system after collision about new stance foot:
hSysAfter = cross2d(p2New,m*dp2New) + I*dq2New + ...
    cross2d(p1New,m*dp1New) + I*dq1New;

% Angular momentum of the new swing leg about the hip
hLegAfter = cross2d(p2New-pHipNew,m*dp2New) + I*dq2New;

% solve the dynamics:
eqnsHs = [hSysBefore-hSysAfter; hLegBefore-hLegAfter];
varsHs = [dq1New; dq2New];
[AA,bb] = equationsToMatrix(eqnsHs, varsHs);
soln.hs = AA\bb;
soln.dq1New = simplify(soln.hs(1));
soln.dq2New = simplify(soln.hs(2));

% Write the heel-strike map to a file:
matlabFunction(q1New,q2New,soln.dq1New,soln.dq2New,...
    'file','autoGen_heelStrike.m',...
    'vars',{q1,q2,dq1,dq2, m, I, d,l},...
    'outputs',{'q1New','q2New','dq1New','dq2New'});

