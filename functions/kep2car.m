function [r_ECI, v_ECI] = kep2car(kep, mu)
% kep2car Transformation from Keplerian elements to cartesian coordinates
%
% Function to transform the data expressed in Keplerian elements to
% cartesian coordinates in Earth Centered Inertial reference frame
%
% INPUT:
% kep [1x6] Vector of the Keplerian elements [a,e,i,Om,om,theta]
% mu  [1x1] Planetary gravitational constant [km^3/s^2]
%
% OUTPUTS:
% r_ECI [nx3] Position vector in ECI frame [km]
% v_ECI [nx3] Velocity vector in ECI frame [km/s]
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

a = kep(1); % Semi-major axis [km]
e = kep(2); % Eccentricity [-]
i = kep(3); % Inclination [-]
Om = kep(4); % RAAN [rad]
om = kep(5); % Argument of perigee [rad]
theta = kep(6); % True anomaly [rad]

p = a*(1-e^2); % Semi-latus rectum [km]
h = sqrt(p*mu); % Angular momentum [km^2/s]
r = p / (1+e*cos(theta)); % Position radius [km]

r_PF = r*[cos(theta), sin(theta), 0]'; % Position vector in perifocal frame [km]
v_PF = (mu/h) * [-sin(theta), (e+cos(theta)), 0]'; % Velocity vector in perifocal frame [km/s]

% Rotation matrices:
R_om = [cos(om)  sin(om)    0   ;
        -sin(om) cos(om)    0   ;
           0        0       1   ];
       
R_i  = [   1        0       0   ;
           0     cos(i)   sin(i);
           0    -sin(i)   cos(i)];
       
R_Om = [cos(Om)  sin(Om)    0   ;
        -sin(Om) cos(Om)    0   ;
           0        0       1   ];
    
R313 = R_om * R_i * R_Om; % ECI -> PF


% PF -> ECI
r_ECI = R313'* r_PF; % Position vector in ECI frame [km]
v_ECI = R313'* v_PF; % Velocity vector in ECI frame [km/s]


end

