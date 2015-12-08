% Derive_cartPole.m
%
% This script derives the equations of motion for a simple cart-pole
% system: A cart that travels on horizontal rails, with a pendulum hanging
% from it. A force actuator is used to move the car along the rails. Both
% the cart and the pendulum are point masses.
%

clear; clc;

syms x q dx dq ddx ddq 'real'   % states
syms u 'real' % actuation
syms m1 m2 g l 'real' % physical parameters

%%%% Unit vectors:
i = sym([1;0]);
j = sym([0;1]);

e = cos(q)*(-j) + sin(q)*(i);    % cart -> tip of pendulum

%%%% State vectors:
z = [x;q;dx;dq];
dz = [dx;dq;ddx;ddq];

%%%% Kinematics:
derivative = @(f)( jacobian(f,z)*dz );   %Chain rule!

p1 = x*i;  %Position of the center of the cart
p2 = p1 + l*e;  %Position of the end of the pendulum

dp1 = derivative(p1);  
dp2 = derivative(p2); 

ddp1 = derivative(dp1); 
ddp2 = derivative(dp2); 

%%%% Define a function for doing '2d' cross product: dot(a x b, k)
cross2d = @(a,b)(a(1)*b(2) - a(2)*b(1));

%%%% Horizontal force balance for entire system (+i direction)
sumForces1 = u;
sumInertial1 = m1*dot(ddp1,i) + m2*dot(ddp2,i);
eqn1 = sumForces1 - sumInertial1;

%%%% Angular momentum balance of pendulum about support point:
sumTorques2 = cross2d(p2-p1,-m2*g*j);   %Gravity torque on pendulum
sumInertial2 = cross2d(p2-p1,m2*ddp2);
eqn2 = sumTorques2 - sumInertial2;

%%%% Write out dynamics in matrix form:   MM*ddq = ff
ddq = [ddx;ddq];
eqns = [eqn1;eqn2];
[MM,ff] = equationsToMatrix(eqns,ddq); 
dyn = simplify(MM\ff);

%%%% Generate an optimized matlab function for dynamics:
matlabFunction(dyn(1),dyn(2),...
    'file','autoGen_cartPoleDynamics.m',...
    'vars',{q,dq,u,m1,m2,g,l},...
    'outputs',{'ddx','ddq'});

%%%% Compute the energy of the system:
U = m2*g*dot(p2,j);   %Potential Energy
T = 0.5*m1*dot(dp1,dp1) + 0.5*m2*dot(dp2,dp2);   %Kinetic Energy

%%%% Generate a matlab function for energy:
matlabFunction(U,T,...
    'file','autoGen_cartPoleEnergy.m',...
    'vars',{x,q,dx,dq,m1,m2,g,l},...
    'outputs',{'U','T'});

%%%% Generate a function for computing the kinematics:
syms empty 'real'  %fixes a bug in matlabFunction related to vectorization
p1(2) = p1(2) + empty;
dp1(2) = dp1(2) + empty;
matlabFunction(p1,p2,dp1,dp2,...
    'file','autoGen_cartPoleKinematics.m',...
    'vars',{x,q,dx,dq,l,empty},...
    'outputs',{'p1','p2','dp1','dp2'});

