function [alpha_vec,delta_vec] = cart2RA_dec(r)
% cart2RA_dec Transform from cartesian coordinates to RA and declination
%
% Function to transform the data expressed in cartesian coordinates to
% Right Ascension and declination vectors
%
% INPUT:
% r [nx3] position vector [km]
%
% OUTPUTS:
% alpha_vec [nx1] Right Ascension vector [rad]
% delta_vec [nx1] Declination vector [rad] 
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

rx=r(:,1); % x component of position vector [km]
ry=r(:,2); % y component of position vector [km]
rz=r(:,3); % z component of position vector [km]

alpha_vec=[];
delta_vec=[];
for i=1:length(rx)
    x=rx(i);
    y=ry(i);
    z=rz(i);
    r_norm=sqrt(x^2+y^2+z^2); % Norm of the position vector [km]
    
    delta=asin(z/r_norm); 
    if y/r_norm>0
    alpha=acos((x/r_norm)/cos(delta));
    else
     alpha=2*pi-acos((x/r_norm)/cos(delta));
    end
    alpha_vec=[alpha_vec;alpha]; % Right Ascension vector [rad]
    delta_vec=[delta_vec;delta]; % Declination vector [rad] 
end
end

