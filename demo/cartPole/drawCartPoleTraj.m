function drawCartPoleTraj(t,p1,p2,nFrame)
% drawCartPoleTraj(t,p1,p2,nFrame)
%
% INPUTS:
%   t =  [1,n] = time stamp for the data in p1 and p2
%   p1 = [2,n] = [x;y] = position of center of the cart
%   p2 = [2,n] = [x;y] = position of tip of the pendulum
%   nFrame = scalar integer = number of "freeze" frames to display
%

clf; hold on;

Cart_Width = 0.15;
Cart_Height = 0.05;

Pole_Width = 4;  %pixels

%%%% Figure out the window size:

[xLow, xUpp, yLow, yUpp] = getBounds(p1,p2);

xLow = xLow - 0.7*Cart_Width;
xUpp = xUpp + 0.7*Cart_Width;

yLow = yLow - 0.7*Cart_Height;
yUpp = yUpp + 0.7*Cart_Height;

Limits = [xLow,xUpp,yLow,yUpp];

%%%% Plot Rails
plot([Limits(1) Limits(2)],-0.5*Cart_Height*[1,1],'k-','LineWidth',2)

%%%% Draw the trace of the pendulum tip
x = p2(1,:);
y = p2(2,:);
z = zeros(size(x));
col = t;  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'FaceColor','no',...
        'EdgeColor','interp',...
        'LineWidth',2);


map = colormap;    
tMap = linspace(t(1),t(end),size(map,1))';

%%%% Compute the frames for plotting:
tFrame = linspace(t(1), t(end), nFrame);
cart = interp1(t',p1',tFrame')';
pole = interp1(t',p2',tFrame')';

for i = 1:nFrame
    
    % Compute color:
    color = interp1(tMap,map,tFrame(i));
    
    %Plot Cart
    x = cart(1,i) - 0.5*Cart_Width;
    y = -0.5*Cart_Height;
    w = Cart_Width;
    h = Cart_Height;
    rectangle('Position',[x,y,w,h],'LineWidth',2);
    
    %Plot Pendulum
    Rod_X = [cart(1,i), pole(1,i)];
    Rod_Y = [cart(2,i), pole(2,i)];
    plot(Rod_X,Rod_Y,'k-','LineWidth',Pole_Width,'Color',color)
    
    %Plot Bob and hinge
    plot(pole(1,i),pole(2,i),'k.','MarkerSize',40,'Color',color)
    plot(cart(1,i),cart(2,i),'k.','MarkerSize',60,'Color',color)
    
end

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