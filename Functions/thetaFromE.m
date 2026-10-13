function [theta] = thetaFromE(E, e)
%THETAFROME Summary of this function goes here
%   Detailed explanation goes here
q = (cos(E) -e) / (1 - e*cos(E));

if q >= 0
    theta = acos(q);
else
    theta = 2*pi - acos(q);
end