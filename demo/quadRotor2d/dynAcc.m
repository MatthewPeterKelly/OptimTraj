function dz = dynAcc(z, u, p)
% dz = dynAcc(z, u, p)
%
% This function computes the dynamics of a simple planar quad-rotor
% helicopter, including chain integrator for acceleration cost function.
%
% INPUTS:
%   z = [9, n] = [X; V1; V2] = state matrix
%   u = [5, n] = [U1; U2] = control matrix
%   p = parameter struct:
%       .g = acceleration due to gravity
%       .d = half distance between rotors
%       .m = mass of each rotor (half-mass of the quad-rotor)
%
% OUTPUTS:
%   dz = derivative of state matrix
%

% Unpack the inputs
X = z(1:3,:);   % configuration
V1 = z(4:6,:);  % rates (dynamics)
% V2 = z(7:9,:);  % rates (integrator)    %unused
U1 = u(1:2,:);  % actuators
U2 = u(3:5,:);  % kinematics

% Call to the actual dynamics function
dV1 = dynQuadRotor(X, U1, p);

% Chain integrators
dX = V1;
dV2 = U2;

% Pack up the outputs
dz = [dX; dV1; dV2];

end