function [fwhh, reduced_chi_square] = fit_tuning_curve(pref_Orate, pref_O, temp_theta_degrees, true_spike_rate, orientation, NeuronID, q, printfigures, AllCells)

%% Tuning Curve Fitting
    %
    %     global Rmax mu
    %     global mu
    %
    %     Rmax = pref_Orate*50;
    %     mu = pref_O*pi/180; let mu be the direction closest to the vector sum
    %
    %     mu = temp_theta_degrees*pi/180; % let mu be the vector sum
    %
    %     f = fittype('vonMises(k,x)'); % function in /Sherlab
    %     f = fittype('vonMises(Rmax,k,x)'); % function in /Sherlab
    f = fittype('vonMises(Rmax,mu,k,x)'); % function in /Sherlab
    
    
    options = fitoptions('Method','NonlinearLeastSquares');
    options.StartPoint = [pref_Orate*50  1 temp_theta_degrees*pi/180];
    options.Lower = [0 0 0];
    
    f = setoptions(f, options);
    
    [fo, gof] = fit((orientation.*pi/180)',true_spike_rate.*50,f); %,'Weights', (1/(error.^2))
    
    %     [f, MSGID] = lastwarn();
    %     warning('off', MSGID)
    
    goodness(:,q) = struct2cell(gof);
    
    MyCoeffs = coeffvalues(fo);
    
    Rmax = MyCoeffs(1);
    mu = MyCoeffs(3);
    MyCoeffs = MyCoeffs(2);
    
    %     fprintf('Rmax= %f, mu = %f, k = %f \n', Rmax, mu, MyCoeffs)
    
    %     test = temp_theta_degrees - mu*180/pi
    
    
    chi_square = 0;
    binomial_error = [];
    predicted_spike_rate = [];
    norm_predicted_spike_rate = [];
    
    ci = [];
    ci = confint(fo);
    
 
    predicted_spike_rate = (Rmax*exp(MyCoeffs*cos((orientation.*(pi/180))-mu)))/exp(MyCoeffs);
   
    norm_predicted_spike_rate = predicted_spike_rate ./ sum(predicted_spike_rate);

    
    % for binomial error, SEM = sqrt(p*(1-p)), estimate p as normalized predicted spike rate,
    %

    binomial_error = sqrt(norm_predicted_spike_rate.*(1-norm_predicted_spike_rate)./5);
    
    
    for i = 1:16
        chi_square = chi_square + ...
            ((predicted_spike_rate(i) - true_spike_rate(i)*50)^2)/predicted_spike_rate(i);
    end
    
    
    for i = 1:16
        chi_square_binom = chi_square + ...
            ((predicted_spike_rate(i) - true_spike_rate(i))^2)/(binomial_error(i)^2);
    end
    
    %      for i = 1:16
    %         chi_square = chi_square + ...
    %             ((norm_predicted_spike_rate(i) - norm_true_spike_rate(i))^2)/(norm_predicted_spike_rate(i));
    %      end
    
    reduced_chi_square = chi_square/15;
    %      reduced_chi_square(q) = chi_square/13;
    
    
    %     reduced_chi_square_binom(q) = chi_square_binom/15;
    
    %      test= 0;
    %     test = sum((predicted_spike_rate' - true_spike_rate).^2);
    %     test
    
    
    fwhh = 2*acos(log((1/2)*exp(MyCoeffs)+(1/2)*exp(-MyCoeffs))/MyCoeffs)  ; % in radians
    
    lowerconf = 2*acos(log((1/2)*exp(ci(1))+(1/2)*exp(-ci(1)))/ci(1));
    upperconf = 2*acos(log((1/2)*exp(ci(2))+(1/2)*exp(-ci(2)))/ci(2));
    
    size_confint(q) = lowerconf-upperconf;
    
    fwhh_fractional_error(q) = size_confint(q)/fwhh; % smaller numbers are good
    
    if printfigures == 5
        
        figure %von Mises fits
        
        h = polar(orientation(:)*(pi/180), true_spike_rate(:)*50, 'ok'); % data points
        set(h, 'MarkerFaceColor', 'r');
        
        hold on
        
        h = polar(0:0.01:2*pi,(Rmax*exp(MyCoeffs*cos((0:0.01:2*pi)-mu)))/exp(MyCoeffs), '-'); % fitted curve
        set(h,'linewidth',4)
        
        
        h = polar([AllCells(q,1), AllCells(q,1)], [0, 100],'k-'); % DS vector arbitrary length
        set(h, 'linewidth', 2)
        
        h = polar(pref_O*pi/180, pref_Orate*50, 'ok'); % determined maximum point
        set(h, 'MarkerFaceColor', 'y');
        
        title({strcat('ID-', num2str(NeuronID), ', X^2=', num2str(reduced_chi_square)),...
            strcat('Pref #Spikes=', num2str(pref_Orate*50), ', CI/FWHH=', num2str(fwhh_fractional_error(q)))},...
            'FontSize', 25);
        
    end
end