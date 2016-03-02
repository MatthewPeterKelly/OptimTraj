function dz = rocketDynamics(z,u)
% dz = rocketDynamics(z,u)
%
% The basic dynamics and drag coefficient data are from the paper:
%
%   "Drag-law Effects in the Goddard Problem"
%   P. Tsiotras, H. Kelley, H.Kelley    1991
%
% INPUTS:
%   z = [3,n] = [h; v; m] = state vector
%   u = [1,n] = [T] = control = thrust
%

h = z(1,:);   %Height
v = z(2,:);   %Velocity
m = z(3,:);   %Mass
T = u;        %Thrust

%%%% Density of air:
% altitude = [0, 1e4, 2e4, 3e4, 4e4];  %(m)  %height above the ground
% density = [1.23, 0.41, 0.089, 0.018, 0.004];   %(kg/m^3)  density of air
density = 1.474085291*(0.9998541833.^h);  %Data fit off of wolfram alpha

%%%% Drag coefficient, calculated from paper:
A1 = 0.0095;
A2 = 25;
A3 = 0.953;
A4 = 0.036;
speedOfSound = 280;  %(m/s)  %At 10 km altitude
mach = v/speedOfSound;
Cd = A1*atan(A2*(mach-A3))+A4;

%%%% Compute the drag:
Area = pi*3.66;  %(m^2) cross-sectional area (SpaceX F9 Falcon)
D = 0.5*Cd.*Area.*density.*v.^2;

%%%% Compute gravity from inverse-square law:
rEarth = 6.3674447e6;  %(m) radius of earth
mEarth = 5.9721986e24;  %(kg) mass of earth
G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
g = G*mEarth./((h+rEarth).^2);

%%%% Complete the calculation:
c = 4500;  %(m/s) rocket exhaust velocity - chemical rocket
dh = v;  %vertical velocity
dv = (T-D)./m - g;   %vertical acceleration
dm = -T/c;   %mass rate

dz = [dh;dv;dm];

end