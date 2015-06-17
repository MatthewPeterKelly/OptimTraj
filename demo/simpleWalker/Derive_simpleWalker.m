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

dp1 = jacobian(p1,z)*dz;  %Chain rule to get velocity of elbow joint
dp2 = jacobian(p2,z)*dz; 

ddp1 = jacobian(dp1,z)*dz;  
ddp2 = jacobian(dp2,z)*dz; 

%%%% Define a function for doing '2d' cross product: dot(a x b, k)
cross2d = @(a,b)(a(1)*b(2) - a(2)*b(1));

%%%% Angular momentum balance of system about shoulder joint:
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


