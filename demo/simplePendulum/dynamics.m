function dx = dynamics(x,u,p)
% dx = dynamics(x,u,p)
%
% Computes the dynamics for the simple pendulum
%

q = x(1,:);
dq = x(2,:);

k = p.k;    c = p.c;
ddq = -c*dq - k*sin(q) + u;
dx = [dq;ddq];

end