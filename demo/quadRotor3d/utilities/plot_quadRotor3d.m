function plot_quadRotor3d(soln)
%
% Plots a solution to quadRotor3d
%
% Input:
%   soln = OptimTraj soln struct
%
% 

ylabels_states = {'x position [m]','y position [m]','z position [m]',...
    'pitch [deg]','roll [deg]','yaw [deg]'...
    'x velocity [m/s]','y velocity [m/s]',' z velocity [m/s]',...
    'pitch rate [deg/s]','roll rate [deg/s]','yaw rate [deg/s]'} ; 
ylabels_controls = {'throttle 1 []','throttle 2 [N]','throttle 3 []','throttle 4 []'} ; 

%
for i = 1:numel(soln)
    % Unpack solution vectors
    t = soln(i).grid.time ; 
    x = soln(i).grid.state ; 
    x(3:6,:) = rad2deg(x(3:6,:)) ; % convert to deg
    x(10:12,:) = rad2deg(x(10:12,:)) ; % convert to deg
    u = soln(i).grid.control ;

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
        plot(t,x(j,:),'bo'); hold on; grid on ; 
        plot(t,x(j,:),'k--');
        xlabel('time [sec]'); ylabel(ylabels_states{j}) ; 
    end
    subplot(6,1,1) % add title
    title(strcat(method_name,', n = ',num2str(nGrid)),'Interpreter','none') ;
    
    % Plot velocities
    figure
    for k=7:12
        subplot(6,1,k-6) % x position
        plot(t,x(k,:),'bo'); hold on; grid on ; 
        plot(t,x(k,:),'k--');
        xlabel('time [sec]'); ylabel(ylabels_states{k}) ; 
    end
    subplot(6,1,1) % add title
    title(strcat(method_name,', n = ',num2str(nGrid)),'Interpreter','none') ;

    % Plot control vector
    figure
    for jj=1:4
        subplot(4,1,jj) % x position
        plot(t,u(jj,:),'bo'); hold on; grid on ; 
        plot(t,u(jj,:),'k--');
        xlabel('time [sec]'); ylabel(ylabels_controls{jj}) ; 
        if jj==3 
            ylim([0, 1.2*max(u(jj,:))]) 
        end
    end
    subplot(4,1,1)
    title(strcat(method_name,', n = ',num2str(nGrid)),'Interpreter','none') ;
end 