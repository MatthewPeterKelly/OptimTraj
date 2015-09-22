% TEST_heelStrike.m
%
% This script is used to test some simple cases of heel-strike, to make
% sure that it produces at least reasonable results.
%

clc; clear;

% Loads the struct of physical parameters (masses, lengths, ...)
p = getPhysicalParameters();

% Pick a test configuration
q = [-0.8; -0.6; -0.2; 0.5; -0.85];
dq = [0;0;0;0;0.5];

% Compute heel-strike
[qNext, dqNext] = heelStrikeMap(q,dq,p);

% Draw the robot to check configuration:
figure(1); 
subplot(1,2,1);
drawRobot(q,p);
title('Before');
subplot(1,2,2);
drawRobot(qNext,p);
title('After');

% Print out velocities:
for i=1:5
    fprintf('dq(%d) = [%4.4f  -->  %4.4f]\n',i,dq(i),dqNext(i));
end




