function [a,e,i,Om,om,theta] = car2kep(rr,vv, mu)
% car2kep Transformation from cartesian coordinates to Keplerian elements
%
% Function to transform the data expressed in cartesian coordinates to 
% the orbital Keplerian elements
%
% INPUT:
% rr  [3x1] position vector [km]
% vv  [3x1] velocity vector [km/s]
% mu  [1x1] gravitational parameter [km^3/s^2]
%
% OUTPUTS:
% a     [1x1] semi-major axis [km]
% e     [1x1] eccentricity [-]
% i     [1x1] inclination [rad]
% Om    [1x1] RAAN [rad]
% om    [1x1] argument of periapsis [rad]
% theta [1x1] true anomaly [rad]
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

r = norm(rr); % Norm of the position vector [km]
v = norm(vv); % Norm of the velocity vector [km/s]
k = [0 0 1]'; % Z axis versor [-]

a = -mu / (v^2 - 2*mu/r); % Semi-major axis [km]

h = cross(rr,vv); % Angular momentum vector [km^2/s]
i = acos(h(3)/norm(h)); % Inclination [rad]
N = cross(k,h) / norm(cross(h,k)); % Line of nodes versor [km^2/s]
e_vect = cross(vv,h)/mu - rr/r; % Eccentricity vector [-]

Om = acos(N(1)); % Right Ascension of the Ascending Node [rad]
if N(2) < 0
    Om = 2*pi - Om;
end

om = acos(dot(N,e_vect)/norm(e_vect)); % Argument of perigee [rad]
if e_vect(3) < 0
    om = 2*pi - om;
end

theta = acos(dot(e_vect,rr)/(r*norm(e_vect))); % True anomaly [rad]
if dot(rr,vv) < 0
    theta = 2*pi - theta;
end

e = norm(e_vect); % Norm of the eccentricity vector [-]

end

