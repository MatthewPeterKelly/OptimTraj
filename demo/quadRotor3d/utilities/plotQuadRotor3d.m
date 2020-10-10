function plotQuadRotor3d(soln,i)
%
% Plots a solution to quadRotor3d
%
% Input:
%   soln = OptimTraj soln struct
%   i = selects which solution to plot. Defaults to last solution in
%           solution struct.
% 
% Written by Conrad McGreal 2020-02-21

% handle optional input arg
if ~exist('i,','var')
    i = numel(soln) ; 
end

% Create labels
ylabels_states = {'x position [m]','y position [m]','z position [m]',...
    'pitch [deg]','roll [deg]','yaw [deg]'...
    'x velocity [m/s]','y velocity [m/s]',' z velocity [m/s]',...
    'pitch rate [deg/s]','roll rate [deg/s]','yaw rate [deg/s]'} ; 
ylabels_controls = {'throttle 1 []','throttle 2 [N]','throttle 3 []','throttle 4 []'} ; 

%% plot solution in three figures
% Unpack solution vectors
plotPoints = 150 ; 
t = linspace(soln(i).grid.time(1), soln(i).grid.time(end), plotPoints);
x = soln(i).interp.state(t);
x(4:6,:) = rad2deg(x(4:6,:)) ; % convert to deg
x(10:12,:) = rad2deg(x(10:12,:)) ; % convert to deg
u = soln(i).interp.control(t);

% grab for title 
method_name = soln(i).problem.options.method ;  % for plot titles
if strcmp(method_name,'trapezoid') 
    nGrid = soln(i).problem.options.(method_name).nGrid ; 
elseif strcmp(method_name,'hermiteSimpson')
    nGrid = soln(i).problem.options.(method_name).nSegment ; 
end 

% Plot positions 
figure
for j=1:6
    subplot(6,1,j) % x position
    plot(t,x(j,:)); hold on; grid on; 
    xlabel('time [sec]'); ylabel(ylabels_states{j}) ; 
end
subplot(6,1,1) % add title
title(strcat(method_name,', n = ',num2str(nGrid)),'Interpreter','none') ;

% Plot velocities
figure
for k=7:12
    subplot(6,1,k-6) % x position
    plot(t,x(k,:)); hold on; grid on ; 
    xlabel('time [sec]'); ylabel(ylabels_states{k}) ; 
end
subplot(6,1,1) % add title
title(strcat(method_name,', n = ',num2str(nGrid)),'Interpreter','none') ;

% Plot control vector
figure
for jj=1:4
    subplot(4,1,jj) % x position
    plot(t,u(jj,:)); hold on; grid on ; 
    xlabel('time [sec]'); ylabel(ylabels_controls{jj}) ; 
    if jj==3 
        ylim([0, 1.2*max(u(jj,:))]) 
    end
end
subplot(4,1,1)
title(strcat(method_name,', n = ',num2str(nGrid)),'Interpreter','none') ;