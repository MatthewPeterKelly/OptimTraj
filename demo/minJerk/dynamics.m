function dz = dynamics(z,u)
% dz = dynamics(z,u)
%
% Computes the first-order form of the dynamics for the combined chain
% integrator and pendulum system
%
% INPUTS:
%   z = [x;v1;v2;a2];
%   u = [u1;u2];
%
% OUTPUTS:
%   dz = dz/dt
%

x = z(1,:);
v1 = z(2,:);
% v2 = z(3,:);   %Unused
a2 = z(4,:);   
u1 = u(1,:);
u2 = u(2,:);

% Pendulum physics
dv1 = pendulum(x,v1,u1);

% Integrator chain physics:
dx = v1;
dv2 = a2;
da2 = u2;

% Combine:
dz = [dx;dv1;dv2;da2];

end