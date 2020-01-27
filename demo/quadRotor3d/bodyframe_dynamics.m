function [out] = bodyframe_dynamics(u, p)
%
% Computes combined (from all motors and props) forces and moments, and resultant acceleration, on a vehicle.
% Does not include effects from gravity.
%
% INPUTS: 
%   u = [N,n] (0-1) = "throttle" vector where N is number of motors.  
%                       All elements should be: 0 < u(i) < 1 
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
n_time = size(u,2) ; 
n_motors = numel(p.propulsion) ; 

% containers to aggregate wrenches from all motors at each timestep
forces = zeros(n_motors,3) ; 
torMoments = zeros(n_motors,3) ; 
thrustMoments = zeros(n_motors,3) ; 
moments = zeros(n_motors,3) ; 
out = zeros(6,n_time) ; % output container

%% body frame: forces and moments due to propulsion.
for j=1:n_time % timesteps
    for i=1:n_motors % actuators
        this_motor = p.propulsion(i) ; 
        this_RPM = u(i,j) * this_motor.maxRPM ; 

        % Compute wrenches on prop (prop frame)
        [thrust, torque] = computePropOpPoint(this_RPM, p.rho, this_motor.d_prop, ...
            this_motor.C_t, this_motor.C_q) ; 

        % Compute prop wrenches in body frame
        forces(i,:) = thrust * this_motor.thrustAxis ; 
        
        % Compute moments due to thrust 
        thrustArm = this_motor.thrustLocation - p.cg ;  % moment arm of motor to CG
        thrustMoments(i,:) = cross(thrustArm,forces(i,:))  ; % moment, due to thrust, on vehicle

        % Compute moments due to countertorque. 
        if this_motor.isSpinDirectionCCW == 1 
           torMoments(i,:) = torque * -this_motor.thrustAxis ; 
        else
           torMoments(i,:) = torque * this_motor.thrustAxis ; 
        end
        
        % Add moments due to countertorque and moments due to thrust
        moments(i,:) = torMoments(i,:) + thrustMoments(i,:) ; 
    end
    
    %% Calculate accelerations (in body frame) based on forces and moments
linAccel = sum(forces)/p.m ; % total linear acceleration from all motors
angAccel = sum(moments) / p.I ;    % total angular acceleration from all motors

    % Prepare output
    out(:,j) = [linAccel' ; angAccel'] ; 
end 
