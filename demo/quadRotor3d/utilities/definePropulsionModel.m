function [quadRotorModel] = definePropulsionModel(quadRotorParams)
% [quadrotorModel] = definePropulsionModel(quadrotorParams)
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
%       .maxThrust [scalar] thrust at 100% throttle (N)
%       .maxRPM [scalar] RPM at 100% throttle (RPM)
%       .maxTorque [scalar] torque at 100% throttle (Nm)
%       .d_prop [scalar] propeller diameter (m)
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

% compute coefficients.
C_t = p.maxThrust / (rho * (p.maxRPM/60)^2 * p.d_prop^4) ; % thrust coefficient
C_q = p.maxTorque / (rho * (p.maxRPM/60)^2 * p.d_prop^5) ; % torque coefficient

%% Assign to struct
quadRotorModel = struct() ; % initialize output struct

for i=1:size(p.thrustLocations,1) 
quadRotorModel(i).thrustAxis = p.thrustAxes(i,:) ;      % [port, nose, top] % nose faces north when body and world axis are aligned (world is "East North Up")
quadRotorModel(i).thrustLocation = p.thrustLocations(i,:) ; 
quadRotorModel(i).isSpinDirectionCCW = p.isSpinDirectionCCW(i,:) ; 
quadRotorModel(i).maxRPM = p.maxRPM ; 
quadRotorModel(i).maxTorque = p.maxTorque ; 
quadRotorModel(i).d_prop = p.d_prop ; 
quadRotorModel(i).C_t = C_t ; 
quadRotorModel(i).C_q = C_q ; 
end