function dz = cartPoleDynamics(z,u,p)
% dz = cartPoleDynamics(z,u,p)
%
% This function computes the first-order dynamics of the cart-pole.
%
% INPUTS:
%   z = [4, n] = [x;q;dx;dq] = state of the system
%   u = [1, n] = horizontal force applied to the cart
%   p = parameter struct
%       .g = gravity
%       .m1 = cart mass
%       .m2 = pole mass
%       .l = pendulum length
% OUTPUTS:
%   dz = dz/dt = time derivative of state
%
%

% x = z(1,:);   %Not used in dynamics
q = z(2,:);
dx = z(3,:);
dq = z(4,:);

[ddx,ddq] = autoGen_cartPoleDynamics(q, dq, u, p.m1, p.m2, p.g, p.l);

dz = [dx;dq;ddx;ddq];

end