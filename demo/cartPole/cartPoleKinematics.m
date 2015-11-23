function [p1,p2,dp1,dp2] = cartPoleKinematics(z,p)
%  [p1,p2,dp1,dp2] = cartPoleKinematics(z,p)
%
% This function computes the mechanical energy of the cart-pole.
%
% INPUTS:
%   z = [4, n] = [x;q;dx;dq] = state of the system
%   p = parameter struct
%       .g = gravity
%       .m1 = cart mass
%       .m2 = pole mass
%       .l = pendulum length
% OUTPUTS:
%   p1 = position of the cart
%   p2 = position of the tip of the pole
%   dp1 = velocity of the cart
%   dp2 = velocity of the tip of the pole
%
%

x = z(1,:);
q = z(2,:);
dx = z(3,:);
dq = z(4,:);

empty = zeros(size(x));  %See note below
[p1,p2,dp1,dp2] = autoGen_cartPoleKinematics(x,q,dx,dq,p.l,empty);

end


% NOTE: Matlab symbolic toolbox does not properly vectorize constant
% expressions in the functions that are created by matlabFunction. This is
% a hack that fixes the problem. The idea is to create a new "variable"
% that is added to the constant expression. The result is a line like:
%   tmpVar = const + empty;
% where const is a scalar, and empty is a vector of zeros that matches the
% desired size of tmpVar. Matlab then correctly parses this by copying the
% const expression onto each element of empty (which is zero)
