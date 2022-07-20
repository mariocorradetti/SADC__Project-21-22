function str = date2string(dateVec)
% date2string Transorm the vector of date in a string vector
% 
% Function to transform the vector indicating the date into a string
% vector
% 
% INPUT:
% dateVec [1x6] Date a 6-element vector [year, month, day, hour, minute, second]
%
% OUTPUT:
% str           Date written in a string vector
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

months = {'January', 'February', 'March', 'April', 'May', 'June', ...
	'July', 'August', 'September', 'October', 'November', 'December'};
str = sprintf('%i %s %i at %02i:%02i:%02i', dateVec(1), months{dateVec(2)}, dateVec(3), dateVec(4), dateVec(5), round(dateVec(6)));

end

