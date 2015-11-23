function data = plotPendulumCart(t,z,u,p)

x = z(1,:);
q = z(2,:);
dx = z(3,:);
dq = z(4,:);

[energy, potential, kinetic] = cartPoleEnergy(z,p);

subplot(3,2,1);
plot(t,x)
ylabel('x')

subplot(3,2,2);
plot(t,q)
ylabel('q')

subplot(3,2,3);
plot(t,dx)
ylabel('dx')

subplot(3,2,4);
plot(t,dq)
ylabel('dq')

subplot(3,2,5)
plot(t,u)
ylabel('u')
xlabel('t')

subplot(3,2,6); hold on;
plot(t,potential,'b-','LineWidth',1)
plot(t,kinetic,'r-','LineWidth',1)
plot(t,energy,'k-','LineWidth',2)
ylabel('energy')
xlabel('t')
legend('potential','kinetic','total');

data.energy = energy;
data.potential = potential;
data.kinetic = kinetic;

data.x = x;
data.q = q;
data.dx = dx;
data.dq = dq;

end