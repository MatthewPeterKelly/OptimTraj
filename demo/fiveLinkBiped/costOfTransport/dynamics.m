function [dState, dStateGrad] = dynamics(time,state,control,p)
% [dState, dStateGrad] = dynamics(time,state,control,p)
%
% Computes the first-order dynamics of the five-link biped, wrapper for
% dynSs.
%
% INPUTS:
%   time = 
%   state = [q; dq; u] = [angles, rates, torques];
%   control = [du;sn;sp] = [torque rate; slack neg; slack pos];
%   p = parameter struct

nt = length(time);
empty = zeros(1,nt);  % For vectorization
q = state(1:5,:);
dq = state(6:10,:);
u = state(11:15,:);
du = control(1:5,:);
ddq = zeros(size(q));

if nargout == 1   % Using numerical gradients
    [m,mi,f,fi] = autoGen_dynSs(...
        q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),...
        dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),...
        p.m1, p.m2, p.m3, p.m4, p.m5, p.I1, p.I2, p.I3, p.I4, p.I5, p.l1, p.l2, p.l3, p.l4, p.c1, p.c2, p.c3, p.c4, p.c5, p.g, empty);
    
    M = zeros(5,5);  %Mass matrix
    F = zeros(5,1);
    for i=1:nt
        M(mi) = m(:,i);
        F(fi) = f(:,i);
        ddq(:,i) = M\F;  %Numerically invert the mass matrix
    end

else  % Using analytic gradients
    [m,mi,f,fi,mz,mzi,mzd,fz,fzi,fzd] = autoGen_dynSs(...
        q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),...
        dq(1,:),dq(2,:),dq(3,:),dq(4,:),dq(5,:),...
        u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),...
        p.m1, p.m2, p.m3, p.m4, p.m5, p.I1, p.I2, p.I3, p.I4, p.I5, p.l1, p.l2, p.l3, p.l4, p.c1, p.c2, p.c3, p.c4, p.c5, p.g, empty);
    
    M = zeros(5,5);  %Mass matrix
    F = zeros(5,1);
    nz = 31;   %Number of dimensions for gradients [t;q;dq;u;du;sn;sp]
    ddqGrad = zeros(5,nz,nt);
    Mz = zeros(mzd); 
    Fz = zeros(fzd);
    for i=1:nt
        M(mi) = m(:,i);
        F(fi) = f(:,i);
        Mz(mzi) = mz(:,i);
        Fz(fzi) = fz(:,i);
        ddq(:,i) = M\F;  %Numerically invert the mass matrix
        for j=1:nz
            ddqGrad(:,j,i) = M\( -Mz(:,:,j)*ddq(:,i) + Fz(:,:,j) );  % Derivative of a matrix inverse http://www.atmos.washington.edu/~dennis/MatrixCalculus.pdf
        end
    end
    
    dqGrad = zeros(5,nz,nt); 
    idxBase_dq = 2 + 5 -1; % time + angle - one index
    for i=1:5
       dqGrad(i,idxBase_dq+i,:) = 1; 
    end
    
    duGrad = zeros(5,nz,nt);
    idxBase_du = 2 + (5 + 5 + 5) -1; % time + angle + rate + torque - one index
    for i=1:5
       duGrad(i,  idxBase_du + i,:) = 1; 
    end
    
    dStateGrad = cat(1,dqGrad,ddqGrad,duGrad);
    
end

dState= [dq;ddq;du];

end