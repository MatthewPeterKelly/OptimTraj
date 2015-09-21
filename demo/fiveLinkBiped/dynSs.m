function ddq = dynSs(q,dq,u,p)
% ddq = dynSs(q,dq,u,p)
%
% This function computes the dynamics of the five-link biped during single
% stance phase of motion (one foot on ground, one in the air).
%
% INPUTS:
%   q = [5,n] = configuration
%   dq = [5,n] = rates
%   u = [5,n] = torque inputs
%   p = paramter struct
%
% OUTPUTS:
%   ddq = [5,n] = accelerations
%

M = zeros(5,5);  %Temporary mass matrix

% Vectorized calculation for mass matrix and generalized forcing terms
[MM,Idx,F] = autoGen_dynSs(...
    q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),...
    dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
    u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),...
    p.m1, p.m2, p.m3, p.m4, p.m5, p.I1, p.I2, p.I3, p.I4, p.I5, p.l1, p.l2, p.l3, p.l4, p.c1, p.c2, p.c3, p.c4, p.c5, p.g);

% Loop through and numerically invert the mass matrix
ddq = zeros(size(q));
nGrid = size(q,2);
for i=1:nGrid
   M(Idx) = MM(:,i);
   f = F(:,i);
   ddq(:,i) = M\f;
end

end