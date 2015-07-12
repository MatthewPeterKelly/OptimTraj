function plotAcrobot(t,z,u,p)

[energy.total, energy.potential, energy.kinetic] = acrobotEnergy(z,p);

subplot(2,3,1);
plot(t,z(1,:))
xlabel('t')
ylabel('q1')
title('link one angle')
subplot(2,3,2);
plot(t,z(2,:))
xlabel('t')
ylabel('q2')
title('link two angle')
subplot(2,3,4);
plot(t,z(3,:))
xlabel('t')
ylabel('dq1')
title('link one rate')
subplot(2,3,5);
plot(t,z(4,:))
xlabel('t')
ylabel('dq2')
title('link two rate')

subplot(2,3,3); hold on
plot(t,energy.total,'k')
plot(t,energy.potential,'r')
plot(t,energy.kinetic,'b')
legend('total','potential','kinetic')
xlabel('t')
ylabel('e');
title('mechanical energy')

subplot(2,3,6)
plot(t,u)
xlabel('t')
ylabel('u')
title('torque between links')



end