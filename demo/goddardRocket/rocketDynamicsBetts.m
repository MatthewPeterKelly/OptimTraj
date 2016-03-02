function dz = rocketDynamicsBetts(z,u)
% dz = rocketDynamicsBetts(z,u)
%
% This problem is taken from the book by John T. Betts:
% "Practical Methods for Optimal Control and Estimation using Nonlinear
% Programming"  2nd Edition. 2010. SIAM.
%
% INPUTS:
%   z = [3,n] = [h; v; m] = state vector
%   u = [1,n] = [T] = control = thrust
%

sigma = 5.4915*10^-5;
h0 = 23800;
g = 32.174;
c = 1580.9425;

h = z(1,:);   %Height
v = z(2,:);   %Velocity
m = z(3,:);   %Mass
T = u;        %Thrust

D = sigma*v.^2.*exp(-h/h0);

dh = v;  %vertical velocity
dv = (T-D)./m - g;   %vertical acceleration
dm = -T/c;   %mass rate

dz = [dh;dv;dm];

end