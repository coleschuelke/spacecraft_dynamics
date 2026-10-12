function [a, e_mag, i, RAAN, omega, theta] = rv2coe(r, v, mu)
    % RV2COE Converts ECI position and velocity vectors to classical orbital elements.
    %
    % Inputs:
    %   r  - 3x1 position vector (km or m)
    %   v  - 3x1 velocity vector (km/s or m/s)
    %   mu - Gravitational parameter (units must match r and v)
    %        Defaults to Earth (398600.4418 km^3/s^2) if not provided.
    %
    % Outputs:
    %   a     - Semi-major axis (same units as r)
    %   e_mag - Eccentricity (dimensionless)
    %   i     - Inclination (radians)
    %   RAAN  - Right ascension of the ascending node (radians)
    %   omega - Argument of perigee (radians)
    %   theta - True anomaly (radians)

    if nargin < 3
        mu = 398600.4418; % Default Earth mu in km^3/s^2
    end

    % Ensure column vectors
    r = r(:);
    v = v(:);

    r_mag = norm(r);
    v_mag = norm(v);

    % Specific angular momentum
    h = cross(r, v);
    h_mag = norm(h);

    % Node vector
    K = [0; 0; 1];
    n = cross(K, h);
    n_mag = norm(n);

    % Eccentricity vector
    e = (1/mu) * ((v_mag^2 - mu/r_mag)*r - dot(r,v)*v);
    e_mag = norm(e);

    % Specific orbital energy and semi-major axis
    energy = (v_mag^2)/2 - mu/r_mag;
    
    if abs(e_mag - 1) > 1e-10
        a = -mu / (2 * energy);
    else
        a = inf; % Parabolic orbit
    end

    % Inclination
    i = acos(h(3) / h_mag);

    % RAAN
    if n_mag > 1e-10
        RAAN = acos(n(1) / n_mag);
        if n(2) < 0
            RAAN = 2*pi - RAAN;
        end
    else
        RAAN = 0; % Equatorial orbit singularity
    end

    % Argument of Perigee
    if n_mag > 1e-10 && e_mag > 1e-10
        omega = acos(dot(n, e) / (n_mag * e_mag));
        if e(3) < 0
            omega = 2*pi - omega;
        end
    else
        omega = 0; % Circular or equatorial singularity
    end

    % True Anomaly
    if e_mag > 1e-10
        theta = acos(dot(e, r) / (e_mag * r_mag));
        if dot(r, v) < 0
            theta = 2*pi - theta;
        end
    else
        % Handle circular orbit singularities
        if n_mag > 1e-10
            % Circular inclined: theta becomes argument of latitude
            theta = acos(dot(n, r) / (n_mag * r_mag));
            if r(3) < 0
                theta = 2*pi - theta;
            end
        else
            % Circular equatorial: theta becomes true longitude
            theta = acos(r(1) / r_mag);
            if r(2) < 0
                theta = 2*pi - theta;
            end
        end
    end
end