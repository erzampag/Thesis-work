function [y_new, RL, DOT, TTP] = fit_time_course(x, y)

% This function fits the time course to timeFilter.m, and calculates
% response latency and degree of transiency. 
%
% ENZ, Winter 2013

printfig = 0;

%% Fitting and extracting coefficients

f = fittype('timeFilter(a_one, a_two, n_one, n_two, tau_one, tau_two, x)');
% f = fittype('timeFilter2(a_one, a_two, n_one, n_two, tau_one, tau_two, x)');

options = fitoptions('Method','NonlinearLeastSquares');

% % % % options.StartPoint = [1 1 15 15 5 5]; can fit 1 2 3

a_one_est = max(y);
a_two_est = min(y);
tau_one_est = round(find(a_one_est == y(:))*25/length(y));
tau_two_est = round(find(a_two_est == y(:))*25/length(y));

if tau_one_est > tau_two_est
    t_one = tau_one_est + 1;
    t_two = tau_two_est - 1;
else
    t_one = tau_one_est - 1;
    t_two = tau_two_est + 1;
end

if t_one < 1
    t_one = 1;
end
if t_two < 1
t_two = 1;
end

new_index_1 = round(t_one*length(y)/25);
new_index_2 = round(t_two*length(y)/25);

n_one_est = log(y(new_index_1)/a_one_est) / (log(t_one*length(y)/25/tau_one_est) + 1 - (new_index_1)/tau_one_est);
n_two_est =  log(y(new_index_2)/a_two_est) / (log(t_two*length(y)/25/tau_two_est) + 1 - (new_index_2)/tau_two_est);

if n_one_est > 500
    n_one_est = 500;
end
if n_two_est > 500
    n_two_est = 500;
end

[a_one_est a_two_est real(n_one_est) real(n_two_est) tau_one_est tau_two_est];

options.StartPoint = [a_one_est a_two_est real(abs(n_one_est)) real(abs(n_two_est)) tau_one_est tau_two_est]; %

% options.Lower = [-500 -500 -500 -500 -500 -500];
% options.Upper = [500 500 500 500 500 500];

options.Lower = [-1 -1 0 0 0 0];
options.Upper = [1 1 500 500 30 30];

f = setoptions(f, options);

[fo, gof] = fit(x(find(x > 0)), y(find(x >0)), f); %,'Weights', (1/(error.^2))

MyCoeffs = coeffvalues(fo);

a_one = MyCoeffs(1);
a_two = MyCoeffs(2);
n_one = MyCoeffs(3);
n_two = MyCoeffs(4);
tau_one = MyCoeffs(5);
tau_two = MyCoeffs(6);

%% New y values of time course

x_new = (1:0.1:25);  %%% new timecourse is always 241 values long

y_new = a_one.*((x_new./tau_one).*exp(1-x_new./tau_one)).^n_one + a_two.*((x_new./tau_two).*exp(1-x_new./tau_two)).^n_two;

if printfig == 1
    figure
    ylim([-0.65 0.65]);
    plot(linspace(-400,0,length(x_new)), y_new, 'k', 'linewidth', 5);
    hold on
    plot(x*(400/24) - 400*25/24, y, 'r.', 'markersize', 30);
    
    % % % plot(x_new, y_new, 'k', 'linewidth', 5);
    % % % hold on
    % % % plot(1:25, y, 'r.', 'markersize', 30);
end

%% New zero crossing
    function y = g(x)
        y = a_one*((x./tau_one).*exp(1-x./tau_one)).^n_one + a_two*((x./tau_two).*exp(1-x./tau_two)).^n_two;
    end
fun = @g;

% % if tau_one < tau_two
% %     x0 = [tau_one tau_two];
% % else
% %     x0 = [tau_two tau_one];
% % end

if t_one < t_two
    x0 = [t_one t_two];
else
    x0 = [t_two t_one];
end

RL_unscaled = fzero(fun, x0);
RL = RL_unscaled*(400/24) - 400*25/24; % convert to ms, assuming [-400:0]

if printfig == 1
    plot([RL RL], [-0.65 0.65]);
    plot((-400:-1),zeros(1,400),':');
    % plot((1:25), zeros(1,25), ':');
end

%% New DOT

area1 = quad(fun,1,RL_unscaled);
area2 = quad(fun,RL_unscaled,25);

DOT = 1 - abs((area1 + area2)/(abs(area1) + abs(area2)));

if area2 < 0
    DOT = -DOT;
end
%% TTP

TTP_unscaled = [];
if area2 > 0
    TTP_unscaled = x_new(find(y_new == max(y_new)));
else
    TTP_unscaled = x_new(find(y_new == min(y_new)));
end

TTP = TTP_unscaled*(400/24) - 400*25/24;


end
