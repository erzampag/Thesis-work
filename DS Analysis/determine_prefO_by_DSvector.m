function [DSI, MaxSpikes, index_of_min, index_of_opp, pref_O, DSI_error] = determine_prefO_by_DSvector(theta, orientation, true_spike_rate, q)

%% New Determination of preferred orientation!
    
     %Convert this theta to degrees temporarily
    temp_theta_degrees = [];
    
    if theta < 0
        temp_theta_degrees = (theta + 2*pi) * 180/pi;
    else
        temp_theta_degrees = theta * 180/pi;
    end
    
% %     keyboard
    
    difference_vector = zeros(1,length(orientation));
    for i = 1 : length(orientation);
        difference_vector(i) = abs(orientation(i) - temp_theta_degrees);
    end
    
    min_diff_vect = min(difference_vector);
    index_of_min = find(difference_vector <= min_diff_vect);
    
    if length(index_of_min) > 1 %%%% ACCOUNT FOR MORE THAN ONE PREF_O HERE
        fprintf('more than one preferred orientation\n')
%         keyboard
        index_of_min;
        index_of_min = index_of_min(1);
        index_of_min;
        
    end
    
    pref_O = orientation(index_of_min); % the pref_O is based off the Elstrott paper
    

    
%     index_of_opp = mod((index_of_min+7), 16) + 1;

        index_of_opp = mod((index_of_min + length(orientation)/2 - 1), length(orientation)) + 1;

    
    
    orth_O=pref_O+90;
    if orth_O > 337.5;
        orth_O = orth_O-360;
    end
    
    orth_O_two = pref_O+270; % for ASI
    if orth_O_two > 337.5
        orth_O_two = orth_O_two - 360;
    end
    
    opp_O=pref_O+180;
    if opp_O > 337.5;
        opp_O = opp_O-360;
    end
    
    
    try
    idpref_O= orientation==pref_O;
    catch
        keyboard
    end
    
    
    
    
    
    
    pref_Orate=true_spike_rate(idpref_O);
    idorth_O= orientation==orth_O;
    orth_Orate=true_spike_rate(idorth_O);
    idopp_O= orientation==opp_O;
    opp_Orate=true_spike_rate(idopp_O);
    
    idorth_O_two= orientation==orth_O_two;
    orth_Orate_two=true_spike_rate(idorth_O_two);
    
    %Calculate orientation selectivity and direction selectivity
    %(basically just the percent error between 2 values for each one)
    OSI = (pref_Orate - orth_Orate)/(pref_Orate + orth_Orate);
    OSI = round(OSI/0.01)*.01;
    
    
    DSI = (pref_Orate - opp_Orate)/(pref_Orate + opp_Orate);
    %     DSI = round(DSI/0.01)*.01;
    
    % Appears to be (pref+opp) / (orth1+orth2) : an attempt at OSI
    pref_Arate = pref_Orate + opp_Orate;
    orth_Arate = orth_Orate + orth_Orate_two;
    ASI = (pref_Arate - orth_Arate)/(pref_Arate + orth_Arate);
    
    MaxSpikes = pref_Orate*50;
    
    
        %% Calculation of DSI error
    % I did propagation of errors on the DSI by assuming a gaussian and poisson
    % distribution of the total number of spikes (summed)
    
    %     p_spikes = sum(combined_spikes{index_of_min}.*5);
    %     n_spikes = sum(combined_spikes{index_of_opp}.*5);
    
    %     p_spikes = unave_true_spike_rate(index_of_min)*10;
    %     n_spikes = unave_true_spike_rate(index_of_opp)*10;
    
    p_spikes = true_spike_rate(index_of_min)*50;
    n_spikes = true_spike_rate(index_of_opp)*50;
    
    DSI_error = (2/(p_spikes+n_spikes)^2)*(n_spikes*p_spikes*(n_spikes+p_spikes))^(1/2);
    
end
    
    