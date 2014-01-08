function y = timeFilter(a_one, a_two, n_one, n_two, tau_one, tau_two, x)

% Time filter function to fit time courses: tau_n are the time constants of
% the filters, a_n are amplitudes, n_n are the number of stages of the
% filters - from Petrusca et. al, 2007.

y = a_one*((x/tau_one).*exp(1-x/tau_one)).^n_one + a_two*((x/tau_two).*exp(1-x/tau_two)).^n_two;






