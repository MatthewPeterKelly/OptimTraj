function [f, fx, fy] = hills(x,y)
% [f, fx, fy] = hills(x,y)
%
% This function computes the integrand of the objective function for the
% toy car problem. 
%
% INPUTS:
%   x = x position of car
%   y = y position of car
%
% OUTPUTS:
%  f = f(x,y) = integrand of the objective function
%
% DERIVATION:
%   syms x y 'real';
%   f = sin(3*x) + cos(2*y) + cos(x.*y);
%   fx = diff(f,x);
%   fy = diff(f,y);
%   disp(['f = ' vectorize(f) ';']);
%   disp(['fx = ' vectorize(fx) ';']);
%   disp(['fy = ' vectorize(fy) ';']);

f = cos(2.*y) + sin(3.*x) + cos(x.*y);
fx = 3.*cos(3.*x) - y.*sin(x.*y);
fy = - 2.*sin(2.*y) - x.*sin(x.*y);

end