function [energy, potential, kinetic] = cartPoleEnergy(z,p)
% [energy, potential, kinetic] = cartPoleEnergy(z,p)
%
% This function computes the mechanical energy of the cart-pole.
%
% INPUTS:
%   z = [4, n] = [x;q;dx;dq] = state of the system
%   p = parameter struct
%       .g = gravity
%       .m1 = cart mass
%       .m2 = pole mass
%       .l = pendulum length
% OUTPUTS:
%   energy = total mechanical energy
%   potential = potential energy
%   kinetic = kinetic energy
%
%

x = z(1,:);
q = z(2,:);
dx = z(3,:);
dq = z(4,:);

[potential, kinetic] = autoGen_cartPoleEnergy(x, q, dx, dq, p.m1, p.m2, p.g, p.l);

energy = potential + kinetic;

end