function [mu_largest] = fit_tuning_curve2(pref_Orate, pref_O, temp_theta_degrees, true_spike_rate, orientation, NeuronID, q, printfigures, AllCells, Ftwo_Fzero, Fone_Fzero)


f = fittype('vonMises2(Rmax,mu,k,Rmax2,mu2,k2,x)'); % function in /Sherlab

options = fitoptions('Method','NonlinearLeastSquares');
options.StartPoint = [max(true_spike_rate)*50 max(true_spike_rate)*50 1 1 0 pi];
options.Lower = [-max(true_spike_rate)*50 -max(true_spike_rate)*50 0 0 0 0];
options.Upper = [max(true_spike_rate)*100 max(true_spike_rate)*100 inf inf 2*pi 2*pi];

f = setoptions(f, options);
% keyboard
try
[fo, gof] = fit((orientation.*pi/180)',true_spike_rate.*50,f); %,'Weights', (1/(error.^2))
catch
    q
    figure
    polar(orientation(:)*(pi/180), true_spike_rate(:)*50, 'ok'); % data points
        set(gcf, 'MarkerFaceColor', 'r');
end
%     [f, MSGID] = lastwarn();
%     warning('off', MSGID)

goodness(:,q) = struct2cell(gof);

MyCoeffs = coeffvalues(fo);


% coeffnames(f);

%     Rmax = MyCoeffs(1);
%     mu = MyCoeffs(3);
%     k = MyCoeffs(2);

Rmax = MyCoeffs(1);
Rmax2 = MyCoeffs(2);
k = MyCoeffs(3);
k2 = MyCoeffs(4);
mu = MyCoeffs(5);
mu2 = MyCoeffs(6);

if Rmax > Rmax2
    mu_largest = mu;
else
    mu_largest = mu2;
end


% % % if printfigures == 7 && Ftwo_Fzero(q) > sqrt(0.03) && Ftwo_Fzero(q) > Fone_Fzero(q) && max(true_spike_rate)*50 > 100
if printfigures == 7 

    
    figure %von Mises fits
    subplot(1,2,1);
        
        h = polar(orientation(:)*(pi/180), true_spike_rate(:)*50, 'ok'); % data points
        set(h, 'MarkerFaceColor', 'r');
        
        hold on
        
        h = polar(0:0.01:2*pi,(Rmax*exp(k*cos((0:0.01:2*pi)-mu)))/exp(k) + (Rmax2*exp(k2*cos((0:0.01:2*pi)-mu2)))/exp(k2), '-'); % fitted curve
        set(h,'linewidth',4)

title({strcat('ID-', num2str(NeuronID))})



end
