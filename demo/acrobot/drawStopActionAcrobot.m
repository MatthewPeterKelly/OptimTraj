function drawStopActionAcrobot(soln,p)

clf; hold on;

length = p.l1+p.l2;
axis equal; axis(length*[-1,1,-1,1]); axis off;

% Start by drawing the trajectory of the tip
nCurve = 250;
t = linspace(soln.grid.time(1), soln.grid.time(end), nCurve);
Z = soln.interp.state(t);
[~,p2] = acrobotKinematics(Z,p);
plot(p2(1,:),p2(2,:),'LineWidth',2,'Color',[0.2,0.2,0.8]);

% Now draw the stop-action frames
nFrame = 10;
t = linspace(soln.grid.time(1), soln.grid.time(end), nFrame);
Z = soln.interp.state(t);

for i=1:nFrame
    z = Z(:,i);
    
    [p1,p2] = acrobotKinematics(z,p);
    pos = [[0;0],p1,p2];
    
    val = 0.3 + 0.7*(i-1)/(nFrame-1);
    color = hsv2rgb([0.0, 0.8, val]);
    
    plot(0,0,'ks','MarkerSize',25,'LineWidth',3)
    plot(pos(1,:),pos(2,:),'Color',color,'LineWidth',3)
    plot(pos(1,1),pos(2,1),'k.','MarkerSize',20)
    plot(pos(1,2),pos(2,2),'k.','MarkerSize',20)
    plot(pos(1,3),pos(2,3),'k.','MarkerSize',40)
    
end

end