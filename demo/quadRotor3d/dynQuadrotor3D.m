function dz = dynQuadRotor3d(z, u, p)
%
% This function computes the dynamics of a 6 DOF quadcopter.
%
% INPUTS:
%   z = [12, n] = [X; dX] = [x, y, z, pitch, roll, yaw, dx, dy, dz, dpitch, droll, dyaw] = state matrix
%   u = [4, n] = [u1; u2; u3; u4] = control matrix
%   p = parameter struct: 
%       .g = [scalar] (m/s/s) acceleration due to gravity
%       .rho = [scalar] (kg/m^3) air density
%       .m = [scalar] (kg) total quadcopter mass 
%       .I = [3x3] (kg-m^2) inertia matrix kg-m^2]
%       .cg = [1,3] (m) center-of-gravity location in body coords.
%       .propulsion(i) = [struct] with fields
%           .thrustAxis = [3x1] (unit vector) in body coordinates
%           .thrustLocation = [3x1] (m) XYZ in body coordinates.
%           .isSpinDirectionCCW = [bool] if true, propeller spin direction 
%                                    is counterclockwise around thrust axis.
%           .C_t = [scalar] () thrust coefficient, see https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node86.html
%           .C_q = [scalar] () torque coefficient
%           .maxRPM = [scalar] (RPM) maximum propeller RPM 
%           .maxTorque = [scalar] (Nm) maximum propeller torque
%           .d_prop = [scalar] (m) propeller diameter
%
% OUTPUTS:
%   dz = [12, n] = [dx; dy; dz; da; db; dg; ddx; ddy; ddz; dda; ddb; ddg] = derivative of state matrix
%

% Unpack the inputs
eul = z(4:6,:)' ; % unpack euler angles
linVel = z(7:9,:) ; % linear velocities
omega = z(10:12,:) ; % unpack angular rates 
dX = z(7:12,:);  % rates
ddX = zeros(size(dX)) ; % container

% Compute bodyframe dynamics
[ddX_body] = dynBodyFrame(omega, u, p) ; 

% unpack bodyframe dynamics 
linAccel_body = ddX_body(1:3,:)' ; % linear acceleration in body frame
angAccel_body = ddX_body(4:6,:)' ; % angular acceleration

% Convert bodyframe acceleration, based on world pose, into world frame accelerations.
BodyRotMat = Euler2RotMat(eul) ;  % create body rotation matrix

for i=1:size(BodyRotMat,3)
    linAccel_world = linAccel_body(i,:) * BodyRotMat(:,:,i)' ; % transform body frame linear acceleration into world frame 
    angAccel_world = angAccel_body(i,:) * BodyRotMat(:,:,i)' ; % transform body frame angular acceleration into world frame 

    % Add acceleration due to gravity
    linAccel_world(3) = linAccel_world(3) + p.g ; 

    % Insert this time step into vector
    ddX(:,i) = [linAccel_world'; angAccel_world'] ; 
end 

% Pack up the outputs
dX = [linVel; omega] ; 
dz = [dX; ddX] ; 