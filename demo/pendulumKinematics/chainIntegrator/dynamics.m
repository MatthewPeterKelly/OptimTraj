function dz = dynamics(z,u)
% dz = dynamics(z,u)
%
% Computes the dynamics of a quadruple chain integrator
%
% INPUTS: 
%   z = [4*nx,nt] = [position;velocity;acceleration;jerk]
%   u = [nx,nt] = snap
%
% OUTPUTS:
%   dz = dz/dt
%

nx = size(u,1);    % dimension of the position space
idx = (nx+1):(4*nx);   % index of velocity, acceleration, and jerk
dz = [z(idx,:);u];  % derivatives

end