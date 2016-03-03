function dz = dynamics(z, u, p)
% dz = dynamics(z, u, p)
%
% This function computes the dynamics of a simple planar quad-rotor
% helicopter.
%
% INPUTS:
%   z = [6, n] = [x; y; q; dx; dy; dq] = state matrix
%   u = [2, n] = [u1; u2] = control matrix
%   p = parameter struct:
%       .g = acceleration due to gravity
%       .d = half distance between rotors
%       .m = mass of each rotor (half-mass of the quad-rotor)
%
% OUTPUTS:
%   dz = [6, n] = [dx; dy; dq; ddx; ddy; ddq] = derivative of state matrix
%

% Unpack the inputs
q = z(3,:);  % Angle of the quad-rotor
u1 = u(1,:);  % rotor 1 force magnitude
u2 = u(2,:);  % rotor 2 force magnitude

% The actual dynamics are done by autoGen_dynamics.m, which is computed
% automatically using the symbolic math toolbox. See "Derive_EoM.m".
[ddx,ddy,ddq] = autoGen_dynamics(q,u1,u2, p.m, p.g, p.d);

% Pack up the outputs
dz = [z(4:6,:); ddx;ddy;ddq];

end