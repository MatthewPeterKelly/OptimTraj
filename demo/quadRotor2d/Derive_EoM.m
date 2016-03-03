% Derive  -- Equations of Motion - 2D (planar) quad-rotor
%
% Model:
%   We will assume that the quadrotor is modeled by two point-masses that
%   are connected by a rigid mass-less bar. There are two motors, each of
%   which apply a force orthogonal to the connecting bar. Assume that air
%   drag is negligable.
%

clc; clear;

%%%% Create symbolic variables:
syms x y q 'real'  % horizontal and vertical position, angle (+i axis = 0)
syms dx dy dq 'real'  % rates
syms ddx ddy ddq 'real' % accelerations
syms u1 u2 'real'  % force applied by rotors
syms d m g 'real'  % distance from center to rotor, rotor mass, gravity


%%%% Unit vectors:

% Inertial reference frame:
i = sym([1;0;0]);
j = sym([0;1;0]);
k = sym([0;0;1]);

% Body reference frame:
e = cos(q)*i + sin(q)*j;  % Unit vector from center to rotor 2
n = -sin(q)*i + cos(q)*j;  % Unit vector orthogonal to e


%%%% Position vectors:
Pc = x*i + y*j;   % Position of the center of the quad rotor
P1 = Pc - d*e;  % Position of rotor 1
P2 = Pc + d*e;  % Position of rotor 2


%%%% Define the state vector and its derivative
z = [x;y;q;dx;dy;dq];  %State
dz = [dx;dy;dq;ddx;ddy;ddq];  %Derivative of state


%%%% Kinematics:

% Define the derivative operator. CHAIN RULE.  
derivative = @(f)( jacobian(f,z)*dz );

% Velocities
dPc = derivative(Pc);
dP1 = derivative(P1);
dP2 = derivative(P2);

% Accelerations
ddPc = derivative(dPc);
ddP1 = derivative(dP1);
ddP2 = derivative(dP2);


%%%% Define all force vectors:
Fg1 = -m*g*j;  %Force of gravity on mass 1
Fg2 = -m*g*j;  %Force of gravity on mass 2
Fu1 = u1*n;  %Force due to rotor 1
Fu2 = u2*n;  %Force due to rotor 1


%%%% Force Balance:
sumForces = Fg1 + Fg2 + Fu1 + Fu2;
massAccel = m*ddPc;
eqn1 = dot(sumForces-massAccel,i);
eqn2 = dot(sumForces-massAccel,j);


%%%% Angular Momentum Balance:   (about center of quad rotor)
sumTorques = ...
    cross(P1-Pc, Fg1) + ...
    cross(P2-Pc, Fg2) + ...    
    cross(P1-Pc, Fu1) + ...
    cross(P2-Pc, Fu2);
angMomentum = ...
    cross(P1-Pc, m*ddP1) + ...
    cross(P2-Pc, m*ddP2);
eqn3 = simplify(dot(sumTorques-angMomentum, k));


%%%% Collect and solve equations:
vars = [ddx;ddy;ddq];  %This is what we want to find  (the accelerations)
eqns = [eqn1;eqn2;eqn3];  %These are the dynamics equations
[M,f] = equationsToMatrix(eqns, vars); % EoM are linear in acceleration

% In this problem, the equations of motion are simple, so we can solve them
% analytically using the "\" command. For more complicated systems, it is
% typically much faster to do this solve numerically at run-time.
soln = simplify(M\f);  
ddxSoln = soln(1);
ddySoln = soln(2);
ddqSoln = soln(3);


%%%% Automatically write out dynamics function:
% This will automatically optimize the function to minimize computation
% time. Not too important here, but makes a big difference for complicated
% functions. It also prevents copy-paste errors.
matlabFunction(...
    ddxSoln, ddySoln, ddqSoln,...
    'file','autoGen_dynamics.m',...
    'vars',{q,u1,u2,m,g,d},...
    'outputs',{'ddx','ddy','ddq'});


%%%% Display the solution to the user:
disp(['ddx = ' ddxSoln]);
disp(['ddy = ' ddySoln]);
disp(['ddq = ' ddqSoln]);




