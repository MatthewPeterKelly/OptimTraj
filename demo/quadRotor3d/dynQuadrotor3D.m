function dz = dynQuadrotor3D(z, u, p)
%
% This function computes the dynamics of a 6 DOF quadcopter.
%
% INPUTS:
%   z = [12, n] = [X; dX] = state matrix
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
n = size(z,2) ; 
dX = z(7:12,:);  % rates

% Compute bodyframe dynamics
[ddX_body] = bodyframe_dynamics(u, p) ; 

% Initialize output container
dz = zeros(size(z)) ;

for i=1:n % for each time vector
    % unpack bodyframe dynamics 
    linAccel_body = ddX_body(1:3,i) ; % linear acceleration in body frame
    angAccel_body = ddX_body(4:6,i) ; % angular acceleration

    % Convert bodyframe acceleration, based on world pose, into world frame accelerations.
    BodyRotMat = Euler2RotMat(z(4,i), z(5,i), z(6,i)) ;  % create body rotation matrix
    linAccel_world = linAccel_body' * BodyRotMat' ; % transform body frame linear acceleration into world frame 
    angAccel_world = angAccel_body' * BodyRotMat' ; % transform body frame angular acceleration into world frame 

    % Add acceleration due to gravity
    linAccel_world(3) = linAccel_world(3) + p.g ; 

    % Pack up the outputs
    ddX = [linAccel_world'; angAccel_world'] ; 
    dz(:,i) = [dX(:,n); ddX] ; 
end