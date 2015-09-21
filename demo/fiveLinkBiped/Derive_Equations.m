%%%% Derive Equations - Five Link Biped Model %%%%
%
% This script derives the equations of motion, as well as some other useful
% equations (kinematics, contact forces, ...) for the five-link biped
% model. 
%
%
% Nomenclature:
%
% - There are five links, which will be numbered starting with "1" for the
% stance leg tibia, increasing as the links are father from the base joint,
% and ending with "5" for the swing leg tibia.
%   1 - stance leg tibia (lower leg)
%   2 - stance leg femur  (upper leg)
%   3 - torso
%   4 - swing leg femur
%   5 - swing leg tibia
%
% - This script uses absolute angles, which are represented with "q"
%
% - Derivatives with respect to time are notated by prepending a "d". For
% example the rate of change in an absolute angle is "dq" and angular
% acceleration would be "ddq"
%
% - Joint positions are given with "P", center of mass positions are "G"
%

clc; clear;

%%%% Link state and derivatives
syms q1 q2 q3 q4 q5 'real' % Absolute link orientations
syms dq1 dq2 dq3 dq4 dq5 'real' % Absolute link angular rates
syms ddq1 ddq2 ddq3 ddq4 ddq5 'real' % Absolute link angular accelerations

%%%% System inputs (controls)
syms u1 'real' % Torque acting on stance leg tibia from ground
syms u2 'real' % torque acting on stance leg femur from stance leg tibia
syms u3 'real' % torque acting on torso from the stance leg femur
syms u4 'real' % torque acting on swing leg femur from torso
syms u5 'real' % torque acting on swing leg tibia from swing leg femur

%%%% Physical paramters 
syms m1 m2 m3 m4 m5 'real' % Link masses
syms c1 c2 c3 c4 c5 'real' % center of mass distance from joint
syms l1 l2 l3 l4 l5 'real' % link length
syms I1 I2 I3 I4 I4 'real' % link moment of inertia about center of mass


