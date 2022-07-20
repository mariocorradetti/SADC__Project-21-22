function [] = Plot_Earth(x,y,z)
% Plot_Earth plot the Earth as an ellipsoid centered in (x,y,z)
%
% INPUTS 
% x [1x1] x coordinate of the center
% y [1x1] y coordinate of the center
% z [1x1] z coordinate of the center
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

C = imread('EarthTexture.jpg'); 
theta=0;
[x, y, z] = ellipsoid(x,y,x, 6378.135/1e3, 6378.135/1e3,6356.750/1e3,1E2);
p=surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]), 'FaceColor', 'texturemap','EdgeColor','none');
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
axis equal;
hold on;
end