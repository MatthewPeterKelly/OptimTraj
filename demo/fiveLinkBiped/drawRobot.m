function drawRobot(q,p)
% drawRobot(q,p)
%
% This function draws the robot with configuration q and parameters p
%
% INPUTS:
%   q = [5, 1] = column vector of a single robot configuration
%   p = parameter struct
%


% Compute the points that will be used for plotting
[P, G] = getPoints(q,p);

x = [0; P(1:2:end,1)];
y = [0; P(2:2:end,1)];

P1 = P(1:2,1);
P2 = P(3:4,1);
P3 = P(5:6,1);
P4 = P(7:8,1);
P5 = P(9:10,1);

G1 = G(1:2,1);
G2 = G(3:4,1);
G3 = G(5:6,1);
G4 = G(7:8,1);
G5 = G(9:10,1);

% Heuristics:
L = (p.l1 + p.l2);  % Maximum extended leg length
xBnd = L*[-1.2,1.2];
yBnd = [-0.2*L, L + p.l3];

% Colors:
colorGround = [118,62,12]/255;
colorStance = [200,60,60]/255;
colorSwing = [60,60,200]/255;
colorTorso = [160, 80, 160]/255;

% Set up the figure
hold off;

% Plot the ground:
plot(xBnd,[0,0],'LineWidth',6,'Color',colorGround);

hold on;

% Plot the links:
plot(x(1:3),y(1:3),'LineWidth',4,'Color',colorStance);
plot(x(3:4),y(3:4),'LineWidth',4,'Color',colorTorso);
plot(x([3,5,6]),y([3,5,6]),'LineWidth',4,'Color',colorSwing);

% Plot the joints:
plot(0, 0,'k.','MarkerSize',30);
plot(P1(1), P1(2),'k.','MarkerSize',30);
plot(P2(1), P2(2),'k.','MarkerSize',30);
plot(P3(1), P3(2),'k.','MarkerSize',30);
plot(P4(1), P4(2),'k.','MarkerSize',30);
plot(P5(1), P5(2),'k.','MarkerSize',30);

% Plot the CoM:
plot(G1(1), G1(2),'ko','MarkerSize',8,'LineWidth',2);
plot(G2(1), G2(2),'ko','MarkerSize',8,'LineWidth',2);
plot(G3(1), G3(2),'ko','MarkerSize',8,'LineWidth',2);
plot(G4(1), G4(2),'ko','MarkerSize',8,'LineWidth',2);
plot(G5(1), G5(2),'ko','MarkerSize',8,'LineWidth',2);

% Format the axis:
axis([xBnd,yBnd]); axis equal; axis off;

end