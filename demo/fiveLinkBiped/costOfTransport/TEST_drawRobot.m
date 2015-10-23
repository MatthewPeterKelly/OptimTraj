% TEST_drawRobot.m
%
% This script is used to test the robot and understand the sign conventions
% for the angles of each link.
%

clc; clear;

% Loads the struct of physical parameters (masses, lengths, ...)
p = getPhysicalParameters();

% Pick a test configuration
q = [...
    -0.3; % stance leg tibia angle
    0.7; % stance leg femur angle
    0.0; % torso angle
    -0.5; % swing leg femur angle
    -0.6]; % swing leg tibia angle

% Draw the robot to check configuration:
figure(1); clf;
drawRobot(q,p);




