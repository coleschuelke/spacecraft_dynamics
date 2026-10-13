function [theta1, theta2] = thetaFromR(r, p, e)
%THETAFROMR Summary of this function goes here
%   Detailed explanation goes here

q = (p/r - 1) / e;

theta1 = acos(q);
theta2 = 2*pi - acos(q);
end