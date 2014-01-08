% function y = vonMises(k,x)
% function y = vonMises(Rmax,k,x)
 function y = vonMises(Rmax,mu,k,x)



% y = response to motion to given stimulus (R)
% x = given stimulus in degrees
% Rmax = max response (pref_Orate)
% mu = preferred direction in radians (pref_O)
% k = concentration parameter accounting for tuning width

% global Rmax mu
% global mu


y = (Rmax*exp(k*cos(x-mu)))/exp(k);  % from Oesch,et al, 2005, Elstrot et al, 2008

% y = (Rmax*exp(k*cos(x-mu)*pi/180))/exp(k);  % from Oesch,et al, 2005, Elstrot et al, 2008


% y = (pref_Orate*exp(k*cos(x-pref_O)*pi/180))/exp(k);  % from Oesch,et al, 2005, Elstrot et al, 2008


