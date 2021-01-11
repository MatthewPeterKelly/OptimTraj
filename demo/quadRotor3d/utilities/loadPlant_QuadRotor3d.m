function [p] = loadPlant_QuadRotor3d()
%
% Convenience function to load 3d quadcopter.
% Provided to simplify main function and modularize plant model definition.
%
% Depends:
%   - definePropulsionModel.m

% Enviromental params
p.g = -9.81 ; % World Coords is XYZ = [East, North, Up], i.e. gravity is a negative number
p.rho = 1.225 ; % air density during flight (kg/m^3) 

% Inertial params
p.m = 5 ; 
p.I = [0.625 0 0; 0 0.625 0; 0 0 1.25] ; % inertia tensor coords: 
p.cg = [0 0 0] ; % (m) location of center of gravity

% control params
p.uMax = 1 ; % maximum throttle setting

% Propulsion system params - shared for all motors:
qRP.d_prop = 0.305*ones(4,1) ; % propeller diameter (m)
qRP.maxThrust = 25*ones(4,1) ; % thrust at 100% throttle (N)
qRP.maxRPM = 10000*ones(4,1) ; % RPM at 100% throttle (RPM)
qRP.maxTorque = ones(4,1) ;  % torque at 100% throttle (Nm)
qRP.thrustLocations = [0.5 0 0; 0 0.5 0; -0.5 0 0; 0 -0.5 0]; % motor locations (each row one motor in coords: [port, nose, top] 
qRP.thrustAxes = repmat([0 0 1],4,1) ; % thrust axes of each motor in coords port, nose, top.
qRP.isSpinDirectionCCW = [1; 0; 1; 0] ; % bool to reverse motor spin direction around 'thrustAxes'.
plotflag = 0 ; 
[p.propulsion] = definePropulsionModel(qRP,plotflag); 