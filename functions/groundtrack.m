function [lonN,latN,alpha,delta]  = groundtrack(r0,v0,time,data,code)
% groundtrack Compute the satellite's ground track
%
% INPUTS:
%  r0    [3x1]    Initial space vector [km]
%  v0    [3x1]    Initial velocity vector [km]
%  time  [1xn]    Time vector to be considered [s]
%  data  [struct] Data of the mission
%  code  [1x1]    Number that indicates if J2 perturbation is considered or
%                 not
%
% OUTPUTS:
%  lonN  [nx1]    Vector of longitude [deg]
%  latN  [nx1]    Vector of latitude  [deg]
%  alpha [nx1]    Vector of right ascension [rad]
%  delta [nx1]    Vector of declinations [rad]
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

% Set data
mu = data.mu; % Earth’s gravitational parameter [km^3/s^2]
we = data.we; % Earth’s rotation velocity [rad/s]
th_G = data.th_g0; % Longitude of Greenwich meridian at initial time [rad]

% Orbit propagation
[a,e,i,OM,om,th] = car2kep(r0,v0,mu);
kep0 = [a e i OM om th];
options = odeset('RelTol',2.22045e-14,'AbsTol',1e-15);

if code == 0 % Unperturbed
    [T,kep] = ode113(@(t,kep) Gauss_eqs(t,kep,data,code),time,kep0,options);
    r=zeros(length(T),3);
    for i=1:length(T)
        kp=kep(i,:);
        r1 = kep2car(kp, mu); % Returns column vector
        r(i,:) = r1';
    end
       [alpha,delta]=cart2RA_dec(r);
end

if code == 1 % Perturbed with J2
    [T,kep] = ode113(@(t,kep) Gauss_eqs(t,kep,data,code),time,kep0,options);
    r = zeros(length(T),3);
    for i=1:length(T)
        kp = kep(i,:);
        r1 = kep2car(kp, mu); % Returns column vector
        r(i,:) = r1';
    end
    [alpha,delta] = cart2RA_dec(r);
end

lat = delta; % Latitude
lon = wrapToPi(alpha-(th_G+we*T)); % Longitude

% Interpolate btw 2 points for the ground track
latN = lat(1);
lonN = lon(1);
ind = 1;
% I have to insert NaN everytime I pass form lat = pi to lat = -pi
for i=2:length(lon)-2
    % CASE IN WHICH I PASS FROM PI TO -PI
    if lon(i)>lon(i+1) && lon(i)>lon (i-1) && lon(i+1)<lon(i+2)
        val = interp1([lon(i) lon(i+1)+2*pi],[lat(i) lat(i+1)],pi); % Interpolate between the 2 successive points and evaluate at pi before shift to -pi
        latN = [ latN ; lat(ind+1:i) ; val ; NaN ; val];  % Last element val because when I pass from pi to -pi the latitude is the same
        lonN = [ lonN ; lon(ind+1:i) ; pi ; NaN ; -pi]; % Last element is -pi because I pass from pi to -pi
        ind = i;
        % CASE WHEN I PASS FROM -PI TO PI
    elseif lon(i)<lon(i+1) && lon(i)<lon(i-1) && lon(i+1)>lon(i+2)
        val = interp1([lon(i) lon(i+1)-2*pi],[lat(i) lat(i+1)],-pi); % Interpolate between the 2 successive points and evaluate at -pi before shift to pi
        latN = [latN ; lat(ind+1 : i) ; val ; NaN ; val]; % Last element val because when I pass from pi to -pi the latitude is the same
        lonN = [lonN ; lon(ind+1 : i) ; -pi ; NaN ;pi]; % Last element is pi because I pass from -pi to pi
        ind = i;
    end
end
latN = [latN ; lat((ind+1):end)]*180/pi; % Complete the vector of latN with the remaining element of the initial vector lat
lonN = [lonN ; lon((ind+1):end)]*180/pi; % Complete the vector of lonN with the remaining element of the initial vector lon
end


