function [r, v] = coe2rv(a, e_mag, i, RAAN, omega, theta, mu)
% COE2RV Converts Classical Orbital Elements to Position and Velocity.
    %
    % Inputs:
    %   a     - Semi-major axis (same units as r)
    %   e_mag - Eccentricity (dimensionless)
    %   i     - Inclination (radians)
    %   RAAN  - Right ascension of the ascending node (radians)
    %   omega - Argument of perigee (radians)
    %   theta - True anomaly (radians)
    %   mu    - Gravitational parameter (units must match r and v)
    %           Defaults to Earth (398600.4418 km^3/s^2) if not provided.
    % Outputs:
    %   r  - 3x1 position vector (km or m)
    %   v  - 3x1 velocity vector (km/s or m/s)
    
    if nargin < 7
        mu = 398600.4418; % Default Earth mu in km^3/s^2
    end
   
    % 1. Compute semi parameter (p)
    p = a * (1 - e_mag^2);
   
    % 2. Calculate position and velocity in the Perifocal (PQW) frame
    r_mag = p / (1 + e_mag * cos(theta));
   
    r_pqw = [r_mag * cos(theta);
             r_mag * sin(theta);
             0];
             
    v_pqw = sqrt(mu / p) * [-sin(theta);
                            e_mag + cos(theta);
                            0];
                            
    % 3. Construct PQW -> ECI rotation matrix using a 3-1-3 Euler sequence
    % Standard ECI to PQW rotation is C = R3(omega) * R1(i) * R3(RAAN).
    
    R3_RAAN = [ cos(RAAN),  sin(RAAN), 0;
               -sin(RAAN),  cos(RAAN), 0;
                0,          0,         1];
                
    R1_i =    [ 1,  0,       0;
                0,  cos(i),  sin(i);
                0, -sin(i),  cos(i)];
                
    R3_omega =[ cos(omega),  sin(omega), 0;
               -sin(omega),  cos(omega), 0;
                0,           0,          1];
    
    % Using R^-1 = R^T
    R_pqw2eci = R3_RAAN' * R1_i' * R3_omega';
    
    % 4. Rotate PQW vectors into the ECI frame
    r = R_pqw2eci * r_pqw;
    v = R_pqw2eci * v_pqw;
    
end