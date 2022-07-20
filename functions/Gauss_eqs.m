function [kep_dot] = Gauss_eqs(t, kep, data, code)
% Gauss_eqs Gauss planetary equations in RSW frame
%
% INPUTS:
%  t      [1x1]    Time [s]
%  kep    [1x6]    Vector of the Keplerian elements [a,e,i,Om,om,theta]
%  data   [struct] Data of the mission
%  code   [1x1]    Number that indicates the perturbation involved
%
% OUTPUTS:
% kep_dot [1x6]    Derivatives of the keplerian elements
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

mu=data.mu; % Planetary gravitational constant [km^3/s^2]

if code == 0 % No perturbations
    a_per = [0;0;0];
elseif code == 1 % Only J2
    [a_J2_rsw, ~] = acc_J2(kep, mu);
    a_per = a_J2_rsw;
elseif code == 2 % Only SRP
    [a_SRP_rsw, ~] = acc_SRP(kep, t, data);
    a_per = a_SRP_rsw;
elseif code == 3 % J2 + SRP
    [a_SRP_rsw, ~] = acc_SRP(kep, t, data);
    [a_J2_rsw, ~] = acc_J2(kep, mu);
    a_per = a_J2_rsw + a_SRP_rsw';
end

a_r = a_per(1); % Radial component of the acceleration [km/s^2]
a_s = a_per(2); % Transversal component of the acceleration [km/s^2]
a_w = a_per(3); % Out-of-plane component of the acceleration [km/s^2]

a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
th = kep(6);

% GAUSS PLANETARY EQS IN RSW FRAME
p = a*(1-e^2); % Semi-latus rectum [km]
h = sqrt(p*mu); % Angular momentum [km^2/s]
r = p/(1+e*cos(th)); % Position radius [km]

kep_dot = zeros(length(kep),1);
kep_dot(1) = 2*a^2/h*(e*sin(th)*a_r + p/r*a_s);
kep_dot(2) = (p*sin(th)*a_r + ((p+r)*cos(th) + r*e)*a_s)/h;
kep_dot(3) = r*cos(th+om)*a_w/h;
kep_dot(4) = r*sin(th+om)*a_w/(h*sin(i));
kep_dot(5) = (-p*cos(th)*a_r + (p+r)*sin(th)*a_s)/(h*e)- kep_dot(4)*cos(i);
kep_dot(6) = h/r^2 + (p*cos(th)*a_r - (p+r)*sin(th)*a_s)/(h*e);

end