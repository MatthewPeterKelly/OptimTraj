function drawRobotStopAction(q,p)
% drawRobotStopAction(q,p)
%
% This function draws the robot with configuration q and parameters p over
% a series of configurations, leaving hold on so that the motion can be
% seen in a single (static) frame.
%
% INPUTS:
%   q = [5, 1] = column vector of a single robot configuration
%   p = parameter struct
%


nTime = size(q,2);

% Compute the points that will be used for plotting
[P, G] = getPoints(q,p);

x = [zeros(1,nTime); P(1:2:end,:)];
y = [zeros(1,nTime); P(2:2:end,:)];

P1 = P(1:2,:);
P2 = P(3:4,:);
P3 = P(5:6,:);
P4 = P(7:8,:);
P5 = P(9:10,:);

% G1 = G(1:2,:);
% G2 = G(3:4,:);
% G3 = G(5:6,:);
% G4 = G(7:8,:);
% G5 = G(9:10,:);

% Heuristics:
L = (p.l1 + p.l2);  % Maximum extended leg length
xBnd = L*[-1.2,1.2];
yBnd = [-0.2*L, L + p.l3];

% Colors:
colorGround = [118,62,12]/255;
colorStance = [200,60,60]/255;
colorSwing = [60,60,200]/255;
colorTorso = [160, 80, 160]/255;

% Plot parameters:
legWidth = 2;
jointSize = 20;

% Set up the figure
hold off;

% Plot the ground:
plot(xBnd,[0,0],'LineWidth',6,'Color',colorGround);

hold on;

% Plot the links:
for i=1:nTime
plot(x(1:3,i),y(1:3,i),'--','LineWidth',legWidth,'Color',colorStance);
plot(0, 0,'k.','MarkerSize',jointSize,'Color',colorStance);
plot(P1(1,i), P1(2,i),'k.','MarkerSize',jointSize,'Color',colorStance);
plot(P2(1,i), P2(2,i),'k.','MarkerSize',jointSize,'Color',colorStance);
% plot(G1(1,i), G1(2,i),'ko','MarkerSize',8,'LineWidth',2);
% plot(G2(1,i), G2(2,i),'ko','MarkerSize',8,'LineWidth',2);
end
for i=1:nTime
plot(x(3:4,i),y(3:4,i),'LineWidth',legWidth+1,'Color',colorTorso);
plot(P3(1,i), P3(2,i),'k.','MarkerSize',jointSize,'Color',colorTorso);
% plot(G3(1,i), G3(2,i),'ko','MarkerSize',8,'LineWidth',2,'Color',colorTorso);
end
for i=1:nTime
plot(x([3,5,6],i),y([3,5,6],i),'LineWidth',legWidth,'Color',colorSwing);
plot(P4(1,i), P4(2,i),'k.','MarkerSize',jointSize,'Color',colorSwing);
plot(P5(1,i), P5(2,i),'k.','MarkerSize',jointSize,'Color',colorSwing);
% plot(G4(1,i), G4(2,i),'ko','MarkerSize',8,'LineWidth',2,'Color',colorSwing);
% plot(G5(1,i), G5(2,i),'ko','MarkerSize',8,'LineWidth',2,'Color',colorSwing);
end

% Format the axis:
axis([xBnd,yBnd]); axis equal; axis off;

end