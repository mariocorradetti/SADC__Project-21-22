function print_GT(lonN, latN, name)
% print_GT Plot the groundtrack on the Earth surface
%
% INPUTS: 
%  lonN  [nx1]   Vector of longitude [deg]
%  latN  [nx1]   Vector of latitude  [deg]
%  name  [str]   Title of the plot
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

x=[-180,180];
y=[-90,90];
A=imread('EarthTexture.jpg');
A=flipud(A);
image(x,y,A);
set(gca,'YDir','normal')
hold on

axis([-180 180 -90 90]);
title('Ground track');
xlabel('Longitude [rad]');
ylabel('Latitude [rad]');
title(name)

hold on
plot(lonN(1:end),latN(1:end),'y','LineWidth',1.5);
scatter(lonN(1),latN(1),30,'k','filled')
text(lonN(1),latN(1),'Start','Color','white','Fontsize',14,'horizontalAlignment','left','verticalAlignment','top');
scatter(lonN(end),latN(end),30,'k','filled')
text(lonN(end),latN(end),'Finish','Color','white','Fontsize',14,'horizontalAlignment','right','verticalAlignment','bottom');
end