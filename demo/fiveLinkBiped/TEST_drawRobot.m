% TEST_drawRobot.m
%
% This script is used to test the robot and understand the sign conventions
% for the angles of each link.
%

clc; clear;

% Loads the struct of physical parameters (masses, lengths, ...)
p = getPhysicalParameters();

% Pick a test configuration
q = [-0.25; 0.2; -0.1; -0.35; -0.55];

% Draw the robot to check configuration:
figure(1); 
drawRobot(q,p);




