function [thrust, torque] = computePropOpPoint(RPM, rho, d_prop, C_t, C_q)
% [thrust, torque] = computePropOpPoint(RPM, rho, d_prop, C_t, C_q)
% 
% Calculates thrust and torque for a series of propeller RPM inputs.
%   Note: Other parameters are considered constant 
%   See: https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node86.html
%
%
% INPUTS: 
%   RPM = [1 x n] (RPM) propeller RPM
%   rho = [scalar] (kg/m^3) air density
%   d_prop = [scalar] (m) propeller diameter
%   C_t = [scalar] () propeller thrust coefficient.
%   C_q = [scalar] () propeller torque coefficient.
%
% OUTPUTS:
%   thrust = [scalar] (N) propeller thrust
%   torque = [scalar] (Nm) propeller torque
%
% Written by Conrad McGreal 2020-1-25

revs = RPM/60 ; % convert to revolution per second

% Calculate thrust
thrust = C_t*rho*revs.^2*d_prop^4 ; 

% Calculate torque 
torque = C_q*rho*revs.^2*d_prop^5 ; 