function dy = ode_2body(t, y, mu)
% ode_2body ODE system for the two-body problem
% 
% Function to solve the ODE system for the two-body problem
%
% INPUTS:
% t  [1x1]        Time vector [s]
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

rr = y(1:3); % Position vector [km]
vv = y(4:6); % Velocity vector [km/s]

r = norm(rr); % Norm of the position vector [km]
v = norm(vv); % Norm of the velocity vector [km/s]

% Derivative of the state vector (r' and r'')
dy = [vv(1); vv(2); vv(3); -mu/r^3*rr(1); -mu/r^3*rr(2); -mu/r^3*rr(3)];

end


