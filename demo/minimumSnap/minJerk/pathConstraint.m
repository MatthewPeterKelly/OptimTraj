function [c, ceq] = pathConstraint(z)
% [c, ceq] = pathConstraint(z)
%
% Computes the path constraint so that the chain integrator dynamics match
% the pendulum dynamics.
%
% INPUTS:
%   z = [x;v1;v2;a2];
%
% OUTPUTS:
%   dz = dz/dt
%
% NOTES:
%   dddx = da2   % definition
%   da2 = u2    % dynamics
%   dv2 = a2    % dynamics
%   v1 = v2     % path constraint
%   

% x = z(1,:);  %Unused
v1 = z(2,:);
v2 = z(3,:); 
% a2 = z(4,:); %unused

c = [];
ceq = v1-v2;  

end