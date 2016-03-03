% MAIN  --  Quad-Rotor Simulation
%
% Simulates a quad-rotor using ode45, running a controller that will
% stabilize it to the origin.
%
%

% Dynamics paramters
p.g = 9.81; % (m/s^2) gravity
p.d = 0.3;  % (m) half-width of quad rotor
p.m = 0.2;  % (m) half-mass of the quad rotor

% Controller parameters:
p.wFast = 20;  % (rad/s) - char. freq. of orientation controller    
p.wSlowX = 2;  % (rad/s) - char. freq. of horizontal controller    
p.wSlowY = 5;  % (rad/s) - char. freq. of vertical controller 
p.xi = 1.0;  % (1/1)  -  effective damping ratio in the controller
p.uMax = 5*(p.m*p.g);  % Maximum force available by each rotor

% Initial state and simulation duration
z0 = 2.0*randn(6,1);
tSpan = [0,5];

% Function handles for simulation
ctrlFun = @(z)(  controller(z, p)  );
dynFun = @(t,z)(  dynamics(z, ctrlFun(z), p)  );

% Run the simulation
soln = ode45(dynFun,tSpan,z0);

% Unpack the solution:
t = linspace(tSpan(1), tSpan(2), 150);
z = deval(soln,t);
[u, qRef] = ctrlFun(z);
x = z(1,:);
y = z(2,:);
q = z(3,:);
dx = z(4,:);
dy = z(5,:);
dq = z(6,:);
u1 = u(1,:);
u2 = u(2,:);

% Plot:
figure(1); clf;

subplot(2,2,1); hold on;
plot(tSpan,[0,0],'k--');
plot(t,x);
xlabel('t')
ylabel('x')

subplot(2,2,2); hold on;
plot(tSpan,[0,0],'k--');
plot(t,y);
xlabel('t')
ylabel('y')

subplot(2,2,3); hold on;
plot(t,qRef,'k--');
plot(t,q);
xlabel('t')
ylabel('q')

subplot(2,2,4); hold on;
plot(tSpan,p.uMax*[1,1],'k--');
plot(tSpan,-p.uMax*[1,1],'k--');
plot(t,u1);  plot(t,u2);
xlabel('t')
ylabel('u')
legend('u1','u2');


