function [] = Plot_Sun
% Plot_Sun plot the Sun as an ellipsoid centered in the orgin of the axis
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

C = imread('Sun.jpg'); 
theta=0;
R=astroConstants(3); % Sun mean radius [km]
[x, y, z] = ellipsoid(0,0,0,R*30,R*30,R*30,1E2);
p=surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]), 'FaceColor', 'texturemap','EdgeColor','none');
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
axis equal;
hold on;
end