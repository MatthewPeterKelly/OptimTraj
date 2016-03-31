function dz = dynamics(z, u, p)
% dz = dynamics(z, u, p)
%
% This function computes the dynamics of a simple planar quad-rotor
% helicopter.
%
% INPUTS:
%   z = [6, n] = [X; dX] = state matrix
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
X = z(1:3,:);   % configuration
dX = z(4:6,:);  % rates

% Call to the actual dynamics function
ddX = dynQuadRotor(X, u, p);

% Pack up the outputs
dz = [dX; ddX];

end