function [a_J2_rsw, a_J2_xyz] = acc_J2(kep, mu)
% acc_J2 calculate the accelerations due to J2 effect
%
% Function that calculates the perturbing acceleration due to J2 effect as
% function of cartesian position (x,y,z frame) and as function of keplerian
% elements, in RSW frame.
% 
% INPUTS:
%  kep      [1x6] Vector of the Keplerian elements [a,e,i,Om,om,theta]
%  mu       [1x1] Planetary gravitational constant [km^3/s^2]
% 
% OUTPUTS:
%  a_J2_rsw [1x3] Perturbing acceleration due to J2 in RSW frame [km/s^2]
%  a_J2_xyz [1x3] Perturbing acceleration due to J2 in xyz frame [km/s^2]
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

J2 = astroConstants(9); % Gravitatonal field constant of the Earth
R_e = 6378.1363; % Radius of the Earth [km]

a = kep(1); e = kep(2); 
i = kep(3); OM = kep(4);
om = kep(5); th = kep(6);

p = a*(1-e^2); % Semi-latus rectum [km]
r = p/(1+e*cos(th)); % Position radius [km]

% RSW 
k = -3/2*J2*mu*R_e^2/r^4;
a_J2_rsw(1) = k*(1-3*((sin(i))^2*(sin(th+om))^2)); % Acceleration along x [km/s^2]
a_J2_rsw(2) = k*((sin(i))^2*sin(2*(th+om))); % Acceleration along y [km/s^2]
a_J2_rsw(3) = k*(sin(2*i)*sin(th+om)); % Acceleration along z [km/s^2]

% Cartesian coordinates 
[rr,~] = kep2car(kep,mu); % Position vector [km]
r = norm(rr); % Norm of the position vector [km]
a_J2_xyz(1) = -k * rr(1)/r*(5*rr(3)^2/r^2-1); % Acceleration along x [km/s^2]
a_J2_xyz(2) = -k * rr(2)/r*(5*rr(3)^2/r^2-1); % Acceleration along y [km/s^2]
a_J2_xyz(3) = -k * rr(3)/r*(5*rr(3)^2/r^2-3); % Acceleration along z [km/s^2]
 
end