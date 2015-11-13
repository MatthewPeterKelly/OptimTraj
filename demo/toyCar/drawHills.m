function drawHills(xBnd,yBnd)
% drawHills(xBnd,yBnd)
%
% This function draws the integrand of the objective function as a surface
% function.
%
% INPUTS:
%   xBnd = [xLow, xUpp];
%   yBnd = [yLow, yUpp];
%

nx = 150;
ny = 150;

x = linspace(xBnd(1), xBnd(2), nx);
y = linspace(yBnd(1), yBnd(2), ny);

[xx,yy] = meshgrid(x,y);
zz = hills(xx,yy);

N = 12;
contour(xx,yy,zz,N,'LineWidth',2);

end