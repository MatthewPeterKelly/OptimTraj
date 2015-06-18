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
syms m1 m2 g l 'real' % physical parameters

%%%% Unit vectors:
i = sym([1;0]);
j = sym([0;1]);

e1 = cos(q1)*(-j) + sin(q1)*(i);    % stance foot -> hip
e2 = cos(q2)*(-j) + sin(q2)*(i);    % hip -> swing foot

%%%% State vectors:
z = [q1;q2;dq1;dq2];
dz = [dq1;dq2;ddq1;ddq2];

%%%% Kinematics:
p1 = l*e1;
p2 = p1 + l*e2;

dp1 = jacobian(p1,z)*dz;  %Chain rule to get velocity of hip joint
dp2 = jacobian(p2,z)*dz; 

ddp1 = jacobian(dp1,z)*dz;  
ddp2 = jacobian(dp2,z)*dz; 

%%%% Define a function for doing '2d' cross product: dot(a x b, k)
cross2d = @(a,b)(a(1)*b(2) - a(2)*b(1));

%%%% Angular momentum balance of system about stance foot:
sumTorques1 = cross2d(p1,-m1*g*j) + cross2d(p2,-m2*g*j);
sumInertial1 = cross2d(p1,m1*ddp1) + cross2d(p2,m2*ddp2);
eqn1 = sumTorques1-sumInertial1;

%%%% Angular momentum balance of swing leg about hip joint:
sumTorques2 = cross2d(p2-p1,-m2*g*j) + u;
sumInertial2 = cross2d(p2-p1,m2*ddp2);
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
    'vars',{q1,q2,dq1,dq2,u,m1,m2,g,l},...
    'outputs',{'ddq1','ddq2'});

%%%% Compute the energy of the system:
U = m1*g*dot(p1,j) + m2*g*dot(p2,j);   %Potential Energy
T = 0.5*m1*dot(dp1,dp1) + 0.5*m2*dot(dp2,dp2);   %Kinetic Energy

%%%% Generate an optimized matlab function for energy:
matlabFunction(U,T,...
    'file','autoGen_energy.m',...
    'vars',{q1,q2,dq1,dq2,m1,m2,g,l},...
    'outputs',{'U','T'});

%%%% Generate a function for computing the kinematics:
matlabFunction(p1,p2,dp1,dp2,...
    'file','autoGen_kinematics.m',...
    'vars',{q1,q2,dq1,dq2,l},...
    'outputs',{'p1','p2','dp1','dp2'});



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%         Derive heel-strike map and collision mechanics                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Angular momentum of the system about the new stance foot (old swing foot)
hSysBefore = cross2d(p1-p2,m1*dp1);
% Note: the old swing foot has zero angular momentum because its relative 
% distance is zero, and the old stance foot has no angular momentum because
% its speed is zero.

% Angular momentum of the old stance leg about the hip
hLegBefore = sym(0);  %Since old stance foot has no speed

% Introduce new variables for the state after the collision:
q1New = q2;
q2New = q1;
syms dq1New dq2New  'real'   % angular rates after collision

% Unit vectors after the collision:
e1New = cos(q1New)*(-j) + sin(q1New)*(i);    % stance foot -> hip
e2New = cos(q2New)*(-j) + sin(q2New)*(i);    % hip -> swing foot

% Kinematics:
p1New = l*e1New;
p2New = p1New + l*e2New;

dp1New = jacobian(p1New,[q1New;q2New])*[dq1New;dq2New];  
dp2New = jacobian(p2New,[q1New;q2New])*[dq1New;dq2New];  

% Angular momentum of the system after collision:
hSysAfter = cross2d(p2New,m2*dp2New) + cross2d(p1New,m1*dp1New);

% Angular momentum of the new stance leg about the hip
hLegAfter = cross2d(p2New-p1New,m2*dp2New);

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
    'vars',{q1,q2,dq1,m1,m2},...
    'outputs',{'q1New','q2New','dq1New','dq2New'});

