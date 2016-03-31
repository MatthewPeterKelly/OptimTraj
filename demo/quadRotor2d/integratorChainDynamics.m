function dz = integratorChainDynamics(z,u)
% dz = integratorChainDynamics(z,u)
% 
% This function computes the dynamics of a chain of first-order integrators
% for a 3-dimensional system.
%
% INPUTS:
%   z = [3*n, m] = [pos; vel; accel; jerk; snap; ...] =  state matrix
%   u = [3, m] = derivative of highest derivative in state.
%
% OUTPUTS:
%   dz = derivative of z
%

dz = [z(4:end, :); u];

end