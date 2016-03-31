function [dObj, uStar] = pathObj(u,w)
% [dObj, uStar] = pathObj(u,w)
%
% This function computes the cost function integrand for the trajectory
%
% INPUTS:
%   u = [U1;U2] = [actuator; derivative] = [2+3, n]
%
% OUTPUTS:
%   dObj = [1,n] = path objective
%   uStar = [3,n] = normalized derivatives
%

if nargin == 1
    w = [1,1,1];
end

n = size(u,2);

uStar = zeros(3,n);
uStar(1,:) = w(1)*u(3,:);
uStar(2,:) = w(2)*u(4,:);
uStar(3,:) = w(3)*u(5,:);

dObj = sum(uStar.^2,1);

end