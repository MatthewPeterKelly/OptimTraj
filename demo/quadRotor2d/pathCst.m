function [c, ceq] = pathCst(z)
% [c, ceq] = pathCst(z)
%
% Computes a velocity-matching path constraint
%
% INPUTS:
%   z = [6, n] = [X; V1; V2; ... ] = state matrix
%
% OUTPUTS:
%   c = []
%   ceq = V1 - V2
%

V1 = z(4:6,:);
V2 = z(7:9,:);

c = [];
ceq = V1-V2;

end