function [c, ceq] = pathConstraint(z)
% [c, ceq] = pathConstraint(z)
%
% Computes the path constraint so that the chain integrator dynamics match
% the pendulum dynamics.
%
% INPUTS:
%   z = [x;v1;v2];
%
% OUTPUTS:
%   dz = dz/dt
%
% NOTES:
%   ddx = dv1   % definition
%   v1 = v2     % path constraint
%   dv2 = u2    % dynamics
%

% x = z(1,:); %Unused
v1 = z(2,:);
v2 = z(3,:);  

c = [];
ceq = v1-v2;  

end