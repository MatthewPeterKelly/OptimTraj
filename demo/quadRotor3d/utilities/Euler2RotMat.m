function [R] = Euler2RotMat(eul)
% [R] = Euler2RotMat(theta_x, theta_y, theta_z)
%
% Stand-in for eul2rotm (built-in function MATLAB Robotics System Toolbox)
% Converts Euler angles [rad] to 3x3 Rotation Matrix
%
% This implementation based on:
% S.M.LaValle, Planning Algorithms, Cambridge University Press, 2006, p.98 
% (PDF available here: http://planning.cs.uiuc.edu/ch3.pdf)
% as cited in: http://nghiaho.com/?page_id=846 
%
% Inputs:
%   eul = [n x 3] - Euler Rotation Angles
%
% Outputs:
%   R = [3 x 3 x n] - rotation matrix. 
%
% Written by Conrad McGreal @ 2019
R = zeros(9,size(eul,1)) ; 

for i=1:size(eul,1) 
    theta_x = eul(i,1) ; 
    theta_y = eul(i,2) ; 
    theta_z = eul(i,3) ; 

    cx = cos(theta_x) ; 
    sx = sin(theta_x) ; 
    cy = cos(theta_y) ; 
    sy = sin(theta_y) ; 
    cz = cos(theta_z) ; 
    sz = sin(theta_z) ; 

    X = [1 0 0
        0 cx -sx
        0 sx cx];

    Y = [cy 0 sy
        0 1 0
        -sy 0 cy];

    Z = [cz -sz 0
        sz cz 0
        0 0 1];

    ZY = matrixMultiply(Z,Y) ; 
    R_square = matrixMultiply(ZY,X) ; 
    R(:,i) = [R_square(1)
        R_square(2)
        R_square(3)
        R_square(4)
        R_square(5)
        R_square(6)
        R_square(7)
        R_square(8)
        R_square(9)] ;
end

function [C] = matrixMultiply(A,B)
% Helper function for matrix multiplication
a = A(1) ; 
d = A(2) ; 
g = A(3) ; 
b = A(4) ; 
e = A(5) ; 
h = A(6) ; 
c = A(7) ; 
f = A(8) ; 
i = A(9) ; 
j = B(1) ; 
m = B(2) ; 
p = B(3) ; 
k = B(4) ; 
n = B(5) ; 
q = B(6) ; 
l = B(7) ; 
o = B(8) ; 
r = B(9) ; 

C = [a*j+b*m+c*p, a*k+b*n+c*q, a*l+b*o+c*r
     d*j+e*m+f*p, d*k+e*n+f*q, d*l+e*o+f*r
     g*j+h*m+i*p, g-k+h-n+i-q, g*l+h*o+i*r];
