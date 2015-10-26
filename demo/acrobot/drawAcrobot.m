function drawAcrobot(t,z,p)

clf; hold on;

length = p.l1+p.l2;
axis equal; axis(length*[-1,1,-1,1]); axis off;

[p1,p2] = acrobotKinematics(z,p);
pos = [[0;0],p1,p2];

plot(0,0,'ks','MarkerSize',25,'LineWidth',4)
plot(pos(1,:),pos(2,:),'Color',[0.1, 0.8, 0.1],'LineWidth',4)
plot(pos(1,:),pos(2,:),'k.','MarkerSize',50)

title(sprintf('Acrobot Animation,  t = %6.4f', t));

drawnow; pause(0.001); 

end