function [u, qRef] = controller(z,p)
% [u, qRef] = controller(z,p)
%
% This function implements a non-linear feedback control law that will
% stabilize the quadrotor to the origin (stationary hover).
%
% INPUTS:
%   z = [6, n] = [x; y; q; dx; dy; dq] = state matrix
%   p = parameter struct:
%       .g = acceleration due to gravity
%       .d = half distance between rotors
%       .m = mass of each rotor (half-mass of the quad-rotor)
%
% OUTPUTS:
%   u = [2, n] = [u1; u2] = control matrix
%
% NOTES:
%
%   System Dynamics:
%       ddx = -(sin(q)*(u1 + u2))/m
%       ddy = (u1*cos(q) + u2*cos(q) - 2*g*m)/m
%       ddq = -(u1 - u2)/(2*d*m)
%
%   How does it work?
%       The basic idea here is that we want to find the controls to drive
%       all components of z->0. The system is underactuated (2 motors, 3
%       degrees of freedom), so we have to be a bit clever about it.
%
% Let's start by rewriting the control inputs:
%   U = u1 + u2   
%   T = u2 - u1  
%
%   ddx = -(sin(q)*U)/m
%   ddy = (U*cos(q) - 2*g*m)/m
%   ddq = T/(2*d*m)
%
% Now, we can divide up the system into two parts:
%
% Position:
%       ddx = -(sin(q)*U)/m
%       ddy = (U*cos(q) - 2*g*m)/m
%
% Orientation:
%       ddq = T/(2*d*m)
%
% Here is the key step. We can put a very fast controller on the angle of
% the quadrotor, which essientially let's us prescribe any angle q. Then,
% we pretend that q is actually a control input. Then we can look at
% position and see that we can control x and y if we can pretend that q is
% an input, rather than a state. 
%
% Simple version: Stabilize q to some arbitrary target qRef very quickly.
% Choose qRef, such that it will stabilize the position (x,y) along with a
% judicious choice of U.
%
% Last piece of the puzzle: Compliant control. Suppose that we have a
% simple second order system:
%   ddz = F 
% We can choose F to make z behave like a simple damped oscillator:
%   ddz = F = -k*x + -c*dx
% Finally, we can select k and c to achieve a desired behavior.
%   ddz = F = -wn^2*x + -2*xi*wn*dx
% Where wn is the natural frequence and xi is the damping
%
% We can then use these control laws to compute a desired acceleration and
% then solve the equations of motion backwards for T,U,and q
%
%

% Controller parameters:
wFast = p.wFast;  % (rad/s)   
wSlowX = p.wSlowX;  % (rad/s)
wSlowY = p.wSlowY;  % (rad/s)
xi = p.xi;  % Slightly over-damped
uMax = p.uMax;  % Maximum force available by rotors

% Unpack the state and parameters
x = z(1,:);
y = z(2,:);
q = z(3,:);
dx = z(4,:);
dy = z(5,:);
dq = z(6,:);
g = p.g;
d = p.d;
m = p.m;

% Position Controller (slow)
%
%   ddx = -(sin(q)*U)/m
%   ddy = (U*cos(q) - 2*g*m)/m
%
%   -m*ddx = U*sin(q) === Fx
%   m*ddy + 2*m*g = U*cos(q) ==== Fy
%
ddx = -wSlowX^2*x + -2*xi*wSlowX*dx;
ddy = -wSlowY^2*y + -2*xi*wSlowY*dy;
Fx = -m*ddx;
Fy = m*ddy + 2*m*g;
[qRef,U] = cart2pol(Fy,Fx);

% Fix angle wrapping for qRef
qErr = pi + (q-qRef);
qShift = floor(qErr/(2*pi));
qRef = qRef + qShift*2*pi;

% Orientation Controller (fast)
%       ddq = T/(2*d*m)
ddq = -wFast^2*(q-qRef) + -2*wFast*xi*dq;
T = ddq*(2*d*m);

% Convert back into u1 and u2
%   U = u1 + u2   
%   T = u2 - u1  
u1 = 0.5*(U - T);
u2 = 0.5*(U + T);
u = [u1;u2];

% Implement controller saturation:
iUpp = u>uMax; u(iUpp) = uMax;
iLow = u< -uMax; u(iLow) = -uMax;

end