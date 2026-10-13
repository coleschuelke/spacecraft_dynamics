function [lat, lon, alt] = ecef2lla(r_ecef)
% ECEF2LLA Converts Earth-Centered, Earth-Fixed (ECEF) position to 
% geodetic Latitude, Longitude, and Altitude (LLA) using WGS-84.
% Uses Bowring's 1976 closed-form algorithm for efficiency.
%
% Inputs:
%   r_ecef - 3x1 or 1x3 Position vector in ECEF frame (km)
%
% Outputs:
%   lat    - Geodetic latitude (radians)
%   lon    - Geodetic longitude (radians)
%   alt    - Altitude above WGS-84 ellipsoid (km)

    % WGS-84 defining parameters
    a = 6378.1370;             % Semi-major axis (km)
    f = 1 / 298.257223563;     % Flattening
    b = a * (1 - f);           % Semi-minor axis
    e2 = 2*f - f^2;            % First eccentricity squared
    ep2 = (a^2 - b^2) / b^2;   % Second eccentricity squared

    x = r_ecef(1);
    y = r_ecef(2);
    z = r_ecef(3);

    % Distance from spin axis
    p = sqrt(x^2 + y^2);

    % Longitude is exact
    lon = atan2(y, x);

    % Bowring's parametric latitude approximation
    theta = atan2(z * a, p * b);

    % Geodetic latitude
    lat = atan2(z + ep2 * b * sin(theta)^3, ...
                p - e2  * a * cos(theta)^3);

    % Radius of curvature in the prime vertical
    N = a / sqrt(1 - e2 * sin(lat)^2);

    % Altitude above ellipsoid
    % Avoid division by zero near the poles by using z for high latitudes
    if abs(lat) < pi/4
        alt = p / cos(lat) - N;
    else
        alt = z / sin(lat) - N * (1 - e2);
    end
end