function [a_SRP_rsw, a_SRP_xyz] = acc_SRP(kep_sc, t, data)
% acc_SRP calculate the accelerations due to SRP effect
%
% Function that calculates the perturbing acceleration due to SRP effect as
% function of cartesian position (x,y,z frame) and as function of keplerian
% elements, in RSW frame.
%
% INPUT 
% kep_sc    [1x6]    Vector of the Keplerian elements [a,e,i,Om,om,theta]
% t         [1x1]    Time [s]
% data      [struct] Data of the mission
%
% OUTPUT
% a_SRP_rsw [1x3]    SRP perturbing acceleration in RSW frame [km/s^2]
% a_SRP_xyz [1x3]    SRP perturbing acceleration in xyz frame [km/s^2]
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

Cr = data.Cr;    % Reflectivity
A2m = data.A2m;  % A/m ratio [m^2/kg]
srp = data.srp;  % Solar radiation pressure  N/(m)^2
eps = deg2rad(23.4);     % Obliquity of the ecliptic [rad]
AU = astroConstants(2);  % Astronomical unit [km]
mu_s = astroConstants(4);  % Sun’s gravitational parameter [km^3/s^2]
Re = data.Re; % Radius of the Earth [km]

t_MJD = data.epoch_MJD + t/86400; % Time [s]
kep_e = uplanet(t_MJD, 3); % Kep of Earth @ Sun (sun centered ecliptic inertial frame)
[rr_e, ~] = kep2car(kep_e, mu_s); % Earth's position vector [km]

R1_eps = [1   0   0; % Rotation matrix from equatorial to ecliptic frame
    0  cos(eps)   sin(eps);
    0  -sin(eps)  cos(eps)];
rr_s2e = R1_eps' * rr_e(:); % Position vector from sun to earth wrt equatorial frame [km]
r = norm(rr_s2e); % Norm of the position vector [km]

[rr_e2sc, ~] = kep2car(kep_sc, data.mu); % Position vector from earth to sc (equatorial) [km]
rr_s2sc = rr_s2e + rr_e2sc; % Position vector from sun 2 sc (equatorial) [km]

% Check on eclipse --------------------------------------------------------
alpha = acos(dot(rr_s2e, rr_s2sc)/norm(rr_s2e)/norm(rr_s2sc));  % Angle between sun-earth and sun-sc [rad]
beta = asin(data.Re/norm(rr_s2e));    % Angle between sun-earth and tangent to earth (control angle) [rad]

if beta < alpha && norm(rr_s2sc) > norm(rr_s2e)
    isSRP = 0;  % Eclipse
else 
    isSRP = 1;  % No eclipse
end

% r = norm(rr_e2sc);
% cosphi = dot(rr_e2sc, -rr_s2e)/norm(rr_e2sc)/norm(rr_s2e);
%     if cosphi<0 && r*sqrt(1-cosphi^2) < Re
%         isSRP = 0;    % eclipse
%     else
%         isSRP = 1;    % no eclipse
%     end

%---------------------------------------------

r = norm(rr_s2sc);
a_SRP_xyz = isSRP*srp*10^(-3)*(AU/r)^2*Cr*A2m* (rr_s2sc./r); % No minus sign 
% since rr_s2sc is already directed from sun to s/c (not from s/c to sun as in the formula)


% RSW FRAME
i = kep_sc(3); OM = kep_sc(4); om = kep_sc(5); theta = kep_sc(6);

% Rotation matrix
R_om = [cos(om+theta)  sin(om+theta)    0   ;
        -sin(om+theta) cos(om+theta)    0   ;
           0        0       1   ];
       
R_i  = [   1        0       0   ;
           0     cos(i)   sin(i);
           0    -sin(i)   cos(i)];
       
R_Om = [cos(OM)  sin(OM)    0   ;
        -sin(OM) cos(OM)    0   ;
           0        0       1   ];

R313 = R_om * R_i * R_Om; % ECI -> RSW
a_SRP_rsw = R313*a_SRP_xyz;

end

