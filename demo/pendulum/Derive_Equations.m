% Derive Equations - Pendulum
%
% This script derives the equations of motion for a simple pendulum, as
% well as the gradients of the dynamics.
%

clc; clear;

syms t 'real'
syms q dq ddq 'real'
syms u 'real'
syms m g l c 'real'

%%%% Dynamics

sumTorques = u - m*g*l*sin(q) - c*dq;
sumInertia = m*l*l*ddq;
ddqSoln = solve(sumTorques-sumInertia,ddq);

z = [t;q;dq;u];  %Lifted state space (state + input)
F = [dq; ddqSoln];   %First-order dynamics
Fz = jacobian(F,z);  %Gradient of the dynamics
Fz = reshape(Fz,numel(Fz),1);

%%%% Objective function

C = u^2;
Cz = jacobian(C,z);
Cz = reshape(Cz,numel(Cz),1);

%%%% Write files:
syms empty 'real'   %Stupid hack for vectorization.

matlabFunction(F+empty,Fz+empty,...
    'file','autoGen_dynamics.m',...
    'vars',{t,q dq u m g l c empty},...
    'outputs',{'F','Fz'});

matlabFunction(C+empty,Cz+empty,...
    'file','autoGen_objective.m',...
    'vars',{t,q dq u m g l c empty},...
    'outputs',{'C','Cz'});

