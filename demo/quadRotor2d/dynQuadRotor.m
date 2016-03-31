function ddX = dynQuadRotor(X, u, p)
% ddX = dynamics(X, u, p)
%
% This function computes the second-order form of the quad-rotor dynamics
%
% INPUTS:
%   X = [3, n] = [x; y; q] = configuration matrix
%   u = [2, n] = [u1; u2] = control matrix
%   p = parameter struct:
%       .g = acceleration due to gravity
%       .d = half distance between rotors
%       .m = mass of each rotor (half-mass of the quad-rotor)
%
% OUTPUTS:
%   ddX = [3, n] = [ddx; ddy; ddq] = second derivative of configuration
%

% Unpack the inputs
q = X(3,:);  % Angle of the quad-rotor
u1 = u(1,:);  % rotor 1 force magnitude
u2 = u(2,:);  % rotor 2 force magnitude

% The actual dynamics are done by autoGen_dynamics.m, which is computed
% automatically using the symbolic math toolbox. See "Derive_EoM.m".
[ddx,ddy,ddq] = autoGen_dynamics(q,u1,u2, p.m, p.g, p.d);

% Pack up the outputs
ddX = [ddx;ddy;ddq];

end