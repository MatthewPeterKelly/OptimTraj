function drawCartPoleAnim(~,p,xLow, xUpp, yLow, yUpp)
% drawCartPoleTraj(t,p,xLow, xUpp, yLow, yUpp)
%
% INPUTS:
%   t =  [1,n] = time stamp for the data in p1 and p2
%   p =  [4,n] = [p1;p2];
%

clf; hold on;

Cart_Width = 0.15;
Cart_Height = 0.05;

p1 = p(1:2,:);
p2 = p(3:4,:);


Pole_Width = 4;  %pixels

%%%% Figure out the window size:

xLow = xLow - 0.7*Cart_Width;
xUpp = xUpp + 0.7*Cart_Width;

yLow = yLow - 0.7*Cart_Height;
yUpp = yUpp + 0.7*Cart_Height;

Limits = [xLow,xUpp,yLow,yUpp];

%%%% Get color map for the figure
% map = colormap;
% tMap = linspace(t(1),t(end),size(map,1))';

%%%% Plot Rails
plot([Limits(1) Limits(2)],-0.5*Cart_Height*[1,1],'k-','LineWidth',2)

%%%% Draw the trace of the pendulum tip  (continuously vary color)
% nTime = length(t);
% for i=1:(nTime-1)
%     idx = i:(i+1);
%     x = p2(1,idx);
%     y = p2(2,idx);
%     c = interp1(tMap,map,mean(t(idx)));
%     plot(x,y,'Color',c);
% end

%%%% Compute the frames for plotting:
% tFrame = linspace(t(1), t(end), nFrame);
% cart = interp1(t',p1',tFrame')';
% pole = interp1(t',p2',tFrame')';

cart = p1;
pole = p2;

% for i = 1:nFrame
    
    % Compute color:
    color = [0.2,0.7,0.1];  %interp1(tMap,map,tFrame(i));
    
    %Plot Cart
    x = cart(1) - 0.5*Cart_Width;
    y = -0.5*Cart_Height;
    w = Cart_Width;
    h = Cart_Height;
    hCart = rectangle('Position',[x,y,w,h],'LineWidth',2);
    set(hCart,'FaceColor',color);
    set(hCart,'EdgeColor',0.8*color);
    
    %Plot Pendulum
    Rod_X = [cart(1), pole(1)];
    Rod_Y = [cart(2), pole(2)];
    plot(Rod_X,Rod_Y,'k-','LineWidth',Pole_Width,'Color',color)
    
    %Plot Bob and hinge
    plot(pole(1),pole(2),'k.','MarkerSize',40,'Color',color)
    plot(cart(1),cart(2),'k.','MarkerSize',60,'Color',color)
    
% end

%These commands keep the window from automatically rescaling in funny ways.
axis(Limits);
axis('equal');
axis manual;
axis off;

end



function [xLow, xUpp, yLow, yUpp] = getBounds(p1,p2)
%
% Returns the upper and lower bound on the data in val
%

val = [p1,p2];
xLow = min(val(1,:));
xUpp = max(val(1,:));
yLow = min(val(2,:));
yUpp = max(val(2,:));

end