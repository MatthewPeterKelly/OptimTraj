function [h] = plotPropLoc(d_prop, loc, ax)
% [h] = plotPropLoc(d_prop, loc, ax)
%   
% Plots 3D propeller location and thrust axis, in bodyframe, for visualization.
%
% Inputs:
%   d_prop = [scalar] [m] propeller diameter
%   loc = [3x1] [m] XYZ location of propeller
%   ax  = [3x1] [unit vector] XYZ axis of propeller
%
% Outputs:
%   h = figure axis handle
%
% Written by Conrad McGreal 2020-02-10

%% Prepare input parameters
r = d_prop/2 ; % radius required for plotting

% Want these to be column vectors, not row vectors.
if size(ax,2)>size(ax,1)
    ax=ax';
end

if size(loc,2)>size(loc,1)
    loc=loc';
end

%% Define unit vectors u and v
% u and v define a new coordinate system in a plane perpendicular to ax
a=[1;0;0];
b=[0;1;0];

if isempty(find(cross(a,ax), 1))==1
    a=[0;0;1];
elseif isempty(find(cross(b,ax), 1))==1
    b=[0;0;1];
end
alpha=dot(ax,a)/dot(ax,ax);
u=a-alpha*ax;
v=cross(u,ax);

u=u/sqrt(sum(u.*u)); 
v=v/sqrt(sum(v.*v));

%% Plot the circle
phi = 0:pi/100:2*pi;
 
% transformation
x = loc(1) + r*cos(phi)*u(1)+r*sin(phi)*v(1) ; 
y = loc(2) + r*cos(phi)*u(2)+r*sin(phi)*v(2) ; 
z = loc(3) + r*cos(phi)*u(3)+r*sin(phi)*v(3) ; 

% plot propeller outer disk
h = plot3(x, y, z,'k--'); hold on; grid on; 

% plot thrust axis
avs = 0.25*d_prop ; % axis vector scaling
tip = loc + avs*(ax) ; % tip of axis vector
ax_plot = [loc, tip] ; % plot vector for axis
plot3(ax_plot(1,:),ax_plot(2,:),ax_plot(3,:),'k'); 

% labeling
xlabel('x position [m]'); ylabel('y position [m]'); zlabel('z position [m]') ;
title('Propeller Locations') 
