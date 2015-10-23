function [ceq, ceqGrad] = cst_heelStrike(x0,xF,p)
% [ceq, ceqGrad] = cst_heelStrike(x0,xF,p)
%
% This function computes the heel-strike constraint as well as gradients.
%
% INPUTS:
%   t0 = time at the start of the trajectory
%   x0 = state at the start of the trajectory
%   tF = time at the end of the trajectory
%   xF = state at the end of the trajectory
%
% OUTPUTS:
%   ceq = equality constraint
%   ceqGrad = gradient of equality constraints
%

qm = xF(1:5);
qp = x0(1:5);
dqm = xF(6:10);
dqp = x0(6:10);

if nargout == 1   %numerical gradients
    [m,mi,f,fi] = autoGen_cst_heelStrike(...
        qp(1),qp(2),qp(3),qp(4),qp(5),...  %Angles "plus" - immediately after heel strike
        qm(1),qm(2),qm(3),qm(4),qm(5),...  %Angles "minus" - immediately before heel-strike
        dqm(1),dqm(2),dqm(3),dqm(4),dqm(5),...   %Rates "minus" - immediately before heel-strike
        p.m1, p.m2, p.m3, p.m4, p.m5, p.I1, p.I2, p.I3, p.I4, p.I5, p.l1, p.l2, p.l3, p.l4, p.l5, p.c1, p.c2, p.c3, p.c4, p.c5, 0);
    
    M = zeros(5,5);  %Mass matrix
    F = zeros(5,1);
    M(mi) = m(:,1);
    F(fi) = f(:,1);
    dqpDyn = M\F;  %Numerically invert the mass matrix
    
    ceqPos = qp - qm([5;4;3;2;1]);
    ceqVel = dqpDyn - dqp;
    ceq = [ceqPos;ceqVel];
    
else %Analytic gradients
    
    [m,mi,f,fi,mz,mzi,mzd,fz,fzi,fzd] = autoGen_cst_heelStrike(...
        qp(1),qp(2),qp(3),qp(4),qp(5),...  %Angles "plus" - immediately after heel strike
        qm(1),qm(2),qm(3),qm(4),qm(5),...  %Angles "minus" - immediately before heel-strike
        dqm(1),dqm(2),dqm(3),dqm(4),dqm(5),...   %Rates "minus" - immediately before heel-strike
        p.m1, p.m2, p.m3, p.m4, p.m5, p.I1, p.I2, p.I3, p.I4, p.I5, p.l1, p.l2, p.l3, p.l4, p.l5, p.c1, p.c2, p.c3, p.c4, p.c5, 0);
    
    M = zeros(5,5);  %Mass matrix
    F = zeros(5,1);
    M(mi) = m(:,1);
    F(fi) = f(:,1);
    dqpDyn = M\F;  %Numerically invert the mass matrix
    
    Mz = zeros(mzd);
    Fz = zeros(fzd);
    Mz(mzi) = mz;
    Fz(fzi) = fz;
    
    nz = 22;   %Number of dimensions for gradients [t0;x0;tF;xF]
    dqpDynGrad = zeros(5,nz);
    for j=1:nz
        dqpDynGrad(:,j) = M\( -Mz(:,:,j)*dqpDyn + Fz(:,:,j) );  % Derivative of a matrix inverse http://www.atmos.washington.edu/~dennis/MatrixCalculus.pdf
    end
    
    dqpGrad = zeros(5,nz);
    qpGrad = zeros(5,nz);
    qmGrad = zeros(5,nz);
    for i=1:5
        qpGrad(i,1+i) = 1;
        dqpGrad(i,6+i) = 1;
        qmGrad(i,12+i) = 1;
    end
    
    ceqPos = qp - qm([5;4;3;2;1]);
    ceqVel = dqpDyn - dqp;
    ceq = [ceqPos;ceqVel];
    
    ceqPosGrad = qpGrad - qmGrad([5;4;3;2;1],:);
    ceqVelGrad = dqpDynGrad - dqpGrad;
    ceqGrad = [ceqPosGrad;ceqVelGrad];
    
end

end