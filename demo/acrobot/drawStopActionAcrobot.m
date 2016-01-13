function drawStopActionAcrobot(soln,p)

clf; hold on;

length = p.l1+p.l2;
axis equal; axis(length*[-1,1,-1,1]); axis off;

nFrame = 10;
t = linspace(soln.grid.time(1), soln.grid.time(end), nFrame);
Z = soln.interp.state(t);

for i=1:nFrame
    z = Z(:,i);
    
    [p1,p2] = acrobotKinematics(z,p);
    pos = [[0;0],p1,p2];
    
    val = 0.1 + 0.8*(i-1)/(nFrame-1);
    color = hsv2rgb([0.3333, 0.875, val]);
    
    plot(0,0,'ks','MarkerSize',25,'LineWidth',4)
    plot(pos(1,:),pos(2,:),'Color',color,'LineWidth',4)
    plot(pos(1,1),pos(2,1),'k.','MarkerSize',30)
    plot(pos(1,2),pos(2,2),'k.','MarkerSize',30)
    plot(pos(1,3),pos(2,3),'k.','MarkerSize',50)
    
end

end