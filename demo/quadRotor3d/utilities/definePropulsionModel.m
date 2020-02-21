function [propulsionModel] = definePropulsionModel(quadRotorParams,plotflag)
% [propulsionModel] = definePropulsionModel(quadrotorParams)
%
% Defines a propulsion parameter plant model for use with dynBodyFrame.
% Supports any number of thrusters, n.
%
% First row of 'thrustLocations','thrustAxes','isSpinDirectionCCW' are
% associated with the first motor. 
% Second row of 'thrustLocations','thrustAxes','isSpinDirectionCCW', are associate with the second
% motor, and so forth. 
% 
% Following parameters are the same for each thruster
%   (i.e.assumes same properties of prop-motor combination at all locations):
%       maxThrust, maxRPM, maxTorque, d_prop
%
% Inputs: 
%   quadRotorParams [struct] parameter struct with the following fields
%       .thrustLocations = [n x 3] [m] location of motors 
%       .thrustAxes = [n x 3] [unitvectors] direction of thrust axis
%       .isSpinDirectionCCW = [n x 1] [bool] if set 
%       .maxThrust [n x 1] thrust at 100% throttle (N)
%       .maxRPM [n x 1] RPM at 100% throttle (RPM)
%       .maxTorque [n x 1] torque at 100% throttle (Nm)
%       .d_prop [n x 1] propeller diameter (m)
%
% Output:
%   quadRotorModel = [struct] with n elements, one for each motor.
%
% Written by Conrad McGreal 2020

% give shorter name
p = quadRotorParams ; 

%% Compute propulsion system coefficients
% See: https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node86.html
% air density during propulsion data collection (kg/m^3)
% (used to determine propulsion system coefficients).
rho = 1.225 ;     

%% Assign to struct
propulsionModel = struct() ; % initialize output struct

for i=1:size(p.thrustLocations,1) 
    propulsionModel(i).thrustAxis = p.thrustAxes(i,:) ;      % [port, nose, top] % nose faces north when body and world axis are aligned (world is "East North Up")
    propulsionModel(i).thrustLocation = p.thrustLocations(i,:) ; 
    propulsionModel(i).isSpinDirectionCCW = p.isSpinDirectionCCW(i,:) ; 
    propulsionModel(i).maxRPM = p.maxRPM(i) ; 
    propulsionModel(i).maxTorque = p.maxTorque(i) ; 
    propulsionModel(i).d_prop = p.d_prop(i) ; 

    % compute coefficients.
    C_t = p.maxThrust(i) / (rho * (p.maxRPM(i)/60)^2 * p.d_prop(i)^4) ; % thrust coefficient
    propulsionModel(i).C_t = C_t ; 

    C_q = p.maxTorque(i) / (rho * (p.maxRPM(i)/60)^2 * p.d_prop(i)^5) ; % torque coefficient
    propulsionModel(i).C_q = C_q ; 

end

% plot propulsion model, if desired
if ~exist('plotflag','var')
else
    if plotflag
        showPropulsionModel(propulsionModel) ; 
    end
end