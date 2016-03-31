function dObj = pathObjective(u)
% dObj = pathObjective(u)
%
% Computes the integrand of the objective function, in this case,
% minimizing the integral of acceleration squared
%
% INPUTS:
%   u = [u1;u2];
%
% OUTPUTS:
%   dOBj = ddx.^2
%
% NOTES:
%   ddx = dv1   % definition
%   v1 = v2     % path constraint
%   dv2 = u2    % dynamics
%

% u1 = u(1,:);  % Unused
u2 = u(2,:);

dObj = u2.^2;

end