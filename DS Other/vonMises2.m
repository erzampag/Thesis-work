 function y = vonMises2(Rmax,mu,k,Rmax2,mu2,k2,x)



% y = response to motion to given stimulus (R)
% x = given stimulus in degrees
% Rmax = max response (pref_Orate)
% mu = preferred direction in radians (pref_O)
% k = concentration parameter accounting for tuning width

% global Rmax mu
% global mu


y = (Rmax*exp(k*cos(x-mu)))/exp(k) + (Rmax2*exp(k2*cos(x-mu2)))/exp(k2);  % from Oesch,et al, 2005, Elstrot et al, 2008