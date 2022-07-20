function [err_kep] = plot_kep_variation(time, vect1, vect2, title_string)
% plot_kep_variation Plot the errors between the different propagation
%                    methods
%
% INPUTS:
%  time         [nx1]  Vector of time output of the propagation method [s]
%  vect1        [nx6]  Vector of keplerian elements output of the first
%                      propagation method 
%  vect2        [nx6]  Vector of keplerian elements output of the second
%                      propagation method 
%  title_string [str]  Title of the plot
%
% OUTPUT:
%  err_kep      [nx6]  Error between vect1 and vect2
%
% AUTHORS:
%  Balossi
%  Corradetti
%  Donato
%  Gelosa

err_kep = abs(vect1 - vect2);

figure()
sgtitle(title_string)
subplot(2,3,1)
plot(time./86400, err_kep(:,1))
title({'a'})
ylabel('err_a [km]'); xlabel('time [days]')
subplot(2,3,2)
plot(time./86400, err_kep(:,2))
title({'e'})
ylabel('err_e [-]'); xlabel('time [days]')
subplot(2,3,3)
plot(time./86400, err_kep(:,3))
title({'i'})
ylabel('err_i [rad]'); xlabel('time [days]')
subplot(2,3,4)
plot(time./86400, err_kep(:,4))
title({'\Omega'})
ylabel('err_\Omega [rad]'); xlabel('time [days]')
subplot(2,3,5)
plot(time./86400, err_kep(:,5))
title({'\omega'})
ylabel('err_\omega [rad]'); xlabel('time [days]')
subplot(2,3,6)
plot(time./86400, err_kep(:,6))
title({'\theta'})
ylabel('err_\theta [rad]'); xlabel('time [days]')

end