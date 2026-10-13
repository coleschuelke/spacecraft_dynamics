function [E] = EFromTheta(theta, e)
%EFROMTHETA Summary of this function goes here
%   Detailed explanation goes here


q = (e + cos(theta)) / (1 + e*cos(theta));
if q >= 0
    E = acos(q);
else
    E = 2*pi - acos(q);
end