function [out] = dynBodyFrame(omega, u, p)
%
% Computes combined (from all motors and props) forces and moments, and resultant acceleration, on a vehicle.
% Does not include effects from gravity.
%
% INPUTS: 
%   omega = [3,n] = [rad/s] angular velocity. 
%   u = [N,n] (0-1) = "throttle" vector where N is number of motors.  
%                       All elements should be: 0 < u(i) < 1 
%                       N = number of motors
%                       n = number of time steps
%   rho = [scalar] (kg-m^3) air density
%   p = parameter struct: 
%       .g = [scalar] (m/s/s) acceleration due to gravity
%       .rho = [scalar] (kg/m^3) air density
%       .m = [scalar] (kg) total quadcopter mass 
%       .I = [3x3] (kg-m^2) inertia matrix kg-m^2]
%       .cg = [1,3] (m) center-of-gravity location in body coords.
%       .propulsion(i) = [struct] propulsion parameter structwith fields
%           .thrustAxis = [3x1] (unit vector) in body coordinates
%           .thrustlocation = [3x1] (m) XYZ in body coordinates.
%           .isSpinDirectionCCW = [bool] if true, propeller spin direction 
%                                    is counterclockwise around thrust axis.
%           .C_t = [scalar] () thrust coefficient see https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node86.html
%           .C_q = [scalar] () torque coefficient
%           .maxRPM = [scalar] (RPM) maximum propeller RPM 
%           .maxTorque = [scalar] (Nm) maximum propeller torque
%           .d_prop = [scalar] (m) propeller diameter
%
% OUTPUTS:
%   out = [6, n] = [linearAccelerations; angularAccelerations] where
%               linearAccelerations = [x_accel; y_accel; z_accel] (m/s/s) forces due to thrust
%               angularAccelerations = [pitch_accel; roll_accel; yaw_accel] (rad/s/s) moment due to countertorque.
%
% Written by Conrad McGreal 2020/01/24  

%% prep
n_motors = numel(p.propulsion) ; 
n_time = size(u,2) ; 

force_totals = zeros(3,n_time,n_motors) ; 
moment_totals = zeros(3,n_time,n_motors) ; 

%% body frame: forces and moments due to propulsion.
for i=1:n_motors % actuators
    this_motor = p.propulsion(i) ; 
    this_RPM = u(i,:) * this_motor.maxRPM ; 

    % Compute wrenches on prop (prop frame)
    [thrust, torque] = computePropOpPoint(this_RPM, p.rho, this_motor.d_prop, ...
        this_motor.C_t, this_motor.C_q) ; 

    % Compute prop wrenches in body frame
    forces = this_motor.thrustAxis' * thrust  ; 

    % Compute moments due to thrust 
    thrustArm = this_motor.thrustLocation - p.cg ;  % moment arm of motor to CG
    thrustArms = repmat(thrustArm,size(forces,2),1) ; 
    thrustMoments = cross(thrustArms,forces')  ; % moment, due to thrust, on vehicle

    % Reverse torque direction if required.
    if this_motor.isSpinDirectionCCW == 1 
       torqueAxis = -this_motor.thrustAxis ; 
    else
       torqueAxis = this_motor.thrustAxis ; 
    end

    % Compute moments due to countertorque. 
    torMoments = torque' * torqueAxis ; 

    % Add moments due to countertorque and moments due to thrust
    moments = torMoments + thrustMoments ; 
    
    % save
    force_totals(:,:,i) = forces ; 
    moment_totals(:,:,i) = moments' ; 
end

%% Calculate accelerations (in body frame) based on forces and moments
% linear accelerations
resultantForce = sum(force_totals,3) ; % add together, for each time step, total forces from all motors.
linAccel = resultantForce./p.m ; % acceleration in each linear direction for every timestep

% angular accelerations
resultantMoment = sum(moment_totals,3) ; % add together, for each time step, total moment from all motors.
bias = cross(omega,p.I*omega); % Eulerian bias acceleration term
angAccel = p.I\(resultantMoment-bias) ; % solve Euler's equation for ang eccleration

% Prepare output
out = [linAccel ; angAccel] ; 
