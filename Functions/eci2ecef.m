function [r_ecef, v_ecef] = eci2ecef(r_eci, v_eci, gmst_rad)
% ECI2ECEF Converts state vectors from Earth-Centered Inertial (ECI) 
% to Earth-Centered, Earth-Fixed (ECEF) frame.
%
% Inputs:
%   r_eci    - 3x1 Position vector in ECI frame (km or m)
%   v_eci    - 3x1 Velocity vector in ECI frame (km/s or m/s)
%   gmst_rad - Greenwich Mean Sidereal Time in radians
%
% Outputs:
%   r_ecef   - 3x1 Position vector in ECEF frame
%   v_ecef   - 3x1 Velocity vector in ECEF frame

    % Ensure inputs are column vectors
    r_eci = r_eci(:);
    v_eci = v_eci(:);

    % WGS-84 Earth rotation rate (rad/s)
    omega_e = 7.2921151467e-5; 
    omega_vec = [0; 0; omega_e];

    % Z-axis rotation matrix for the given sidereal time
    c = cos(gmst_rad);
    s = sin(gmst_rad);
    R3 = [ c,  s,  0;
          -s,  c,  0;
           0,  0,  1];

    % Transform position
    r_ecef = R3 * r_eci;

    % Transform velocity (accounting for rotating reference frame)
    % v_ecef = R3 * (v_eci - omega x r_eci)
    v_ecef = R3 * (v_eci - cross(omega_vec, r_eci));
end