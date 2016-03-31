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

% u1 = u(1,:);  % Unused
u2 = u(2,:);

dObj = u2.^2;

end