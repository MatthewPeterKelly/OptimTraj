function showPropulsionModel(propulsionModel)
%
% Shows propulsion model details: (See function 'TEST_definePropulsionModel.m')
%   - Location of propellers in aircraft body frame
%   - Plots for: 
%       RPM vs throttle
%       Thrust vs RPM
%       Torque vs RPM
%   
% Input:
%   propulsionModel (See function 'definePropulsionModel.m')
%   
% Written by Conrad McGreal

%% plot propellers
n_motors = numel(propulsionModel) ; 
legendcell = cell(n_motors,1) ; 

%% Plot propulsion system geometry
figure 
for i=1:n_motors
    d_prop = propulsionModel(i).d_prop ; 
    loc = propulsionModel(i).thrustLocation ; 
    ax = propulsionModel(i).thrustAxis ; 
    plotPropLoc(d_prop, loc, ax) ; 
    
    % labeling this motor
    x = loc(1) ; 
    y = loc(2) ; 
    z = loc(3) ; 
    dx = 0.1 ; % label offset (from datapoint)
    dy = 0.1 ;
    
    this_label = strcat(['Motor ',num2str(i)]); 
    legendcell{i} = this_label ; 
    if propulsionModel(i).isSpinDirectionCCW
        dir_label = 'Spins CCW' ;
    else
        dir_label = 'Spins CW' ; 
    end
    c1 = cellstr(this_label) ; % make it a cell
    c2 = cellstr(dir_label) ; 
    text(x+dx, y+dy,z, c1) ; % print label to figure
    text(x+dx, y, z, c2) ; % print label to figure
end
title('Propeller Locations - Top view')
xlabel('x position [m]'); ylabel('y position [m]'); zlabel('z position [m]'); 
axis equal
drawnow

%% Plot performance; i.e. Thrust, Torque, & RPM vs throttle
figure
rho = 1.225 ; % set default air density 

for i = 1:n_motors
motorid = i ; % select motor to plot for; at the moment only plots figure for one motor
RPM = 0:100:propulsionModel(motorid).maxRPM ; 
throttle = linspace(0,1,numel(RPM)) ; 
d_prop = propulsionModel(motorid).d_prop ;
C_t = propulsionModel(motorid).C_t ;
C_q = propulsionModel(motorid).C_q ;
[thrust, torque] = computePropOpPoint(RPM, rho, d_prop, C_t, C_q) ; 

% plot it
ax1 = subplot(3,1,1) ; 
plot(throttle,RPM); grid on; hold on; 
xlabel('throttle []'); ylabel('RPM');

subplot(3,1,2) 
plot(RPM,thrust); grid on; hold on; 
xlabel('RPM'); ylabel('thrust [N]');

subplot(3,1,3)
plot(RPM,torque); grid on; hold on; 
xlabel('RPM'); ylabel('torque [Nm]');
end 
titleCell = 'motor model performance details' ; 
title(ax1, titleCell) 
legend(ax1, legendcell,'location','southeast')
drawnow
