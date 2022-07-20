function dy = ode_2body_oblat(t,y,mu)
% ode_2body_oblat ODE system for the perturbed two-body problem (only J2)
%
% Function to solve the ODE system for the perturbed two-body problem
%
% INPUTS:
% t  [1x1]        Time [s]
% y  [6x1]        State vector of the problem
%                 1:3 Position vector [km]
%                 4:6 Velocity vector [km/s]
% mu [1x1]        Planetary gravitational constant [km^3/s^2]
%
% OUTPUT:
% dy [6x1]        Derivative of the state vector
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

J2 = astroConstants(9); % Gravitatonal field constant of the Earth
R_e = 6378.1363; % Radius of the Earth [km]

rr = y(1:3); % Position vector [km]
vv = y(4:6); % Velocity vector [km/s]

r = norm(rr); % Norm of the position vector [km]
v = norm(vv); % Norm of the velocity vector [km/s]

% Perturbations due to second zonal harmonic J2
aJ_x = (3/2)*(J2*mu*R_e^2)/(r^4)*((rr(1)/r)*(5*(rr(3)^2/r^2)-1)); % [km/s^2]
aJ_y = (3/2)*(J2*mu*R_e^2)/(r^4)*((rr(2)/r)*(5*(rr(3)^2/r^2)-1)); % [km/s^2]
aJ_z = (3/2)*(J2*mu*R_e^2)/(r^4)*((rr(3)/r)*(5*(rr(3)^2/r^2)-3)); % [km/s^2]

% Derivative of the state vector (r' and r'')
dy = [vv(1);vv(2);vv(3);-mu/r^3*rr(1)+aJ_x;-mu/r^3*rr(2)+aJ_y;-mu/r^3*rr(3)+aJ_z];

end
