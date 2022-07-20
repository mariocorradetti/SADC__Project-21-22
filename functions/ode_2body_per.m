function dy = ode_2body_per(t, y, data, code)
% ode_2body_per ODE system for the perturbed two-body problem
%
% Function to solve the ODE system for the perturbed two-body problem for
% the cases of no perturbations, J2 effect, SRP effect and both effects.
%
% INPUTS:
% t    [1x1]       Time [s]
% y    [6x1]       State vector of the problem
%                  1:3 Position vector [km]
%                  4:6 Velocity vector [km/s]
% data [struct]    Data of the mission
%
% OUTPUT:
% dy   [6x1]       Derivative of the state vector
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

mu = data.mu; % Planetary gravitational constant [km^3/s^2]

rr = y(1:3); % Position vector [km]
vv = y(4:6); % Velocity vector [km/s]

r = norm(rr); % Norm of the position vector [km]
v = norm(vv); % Norm of the velocity vector [km]

[a,e,i,Om,om,theta] = car2kep(rr,vv, mu);
kep = [a,e,i,Om,om,theta];

if code == 0 % no perturbations
    a_per = [0;0;0];
elseif code == 1 % only J2
    [~, a_J2_xyz] = acc_J2(kep, mu);
    a_per = a_J2_xyz;
elseif code == 2 % only SRP
    [~, a_SRP_xyz] = acc_SRP(kep, t, data);
    a_per = a_SRP_xyz';
elseif code == 3 % J2 + SRP
    [~, a_SRP_xyz] = acc_SRP(kep, t, data);
    [~, a_J2_xyz] = acc_J2(kep, mu);
    a_per = a_J2_xyz + a_SRP_xyz';
end

a_x = a_per(1); % Acceleration component along x [km/s^2]
a_y = a_per(2); % Acceleration component along y [km/s^2]
a_z = a_per(3); % Acceleration component along z [km/s^2]

% Derivative of the state vector (r' and r'')
dy = [vv(1); vv(2); vv(3); -mu/r^3*rr(1)+a_x; -mu/r^3*rr(2)+a_y; -mu/r^3*rr(3)+a_z];

end