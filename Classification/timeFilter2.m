function y = timeFilter2(a_one, a_two, n_one, n_two, tau_one, tau_two, x)


% a_one = 1;
% a_two = 1;
% 
% tau_one = 15;
% tau_two = 15;
% 
% n_one = 5;
% n_two = 10;

% x = [-400:10:0];
% x = (1:25);

y = a_one*(x/tau_one).^(n_one).*exp(-n_one*(x/tau_one-1)) - a_two*(x/tau_two).^(n_two).*exp(-n_two*(x/tau_two-1));


