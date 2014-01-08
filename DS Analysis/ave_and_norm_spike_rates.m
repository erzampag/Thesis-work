function [true_spike_rate, norm_true_spike_rate, Max_Spikes_absolute, max_index] =...
    ave_and_norm_spike_rates(orientation, all_orientations, num_trials, start_times, end_times, total_spikes,q)   
%% Average spike rates
    % calculate average spike rate of cell during the display of each of 16 orientations of a moving bar
% % % %     keyboard
% % %     
% % %     if ndims(total_spikes) > 2
% % %         temp = zeros(length(total_spikes(1,1,:,1)), length(total_spikes(1,1,1,:)));
% % %         
% % %         for i  = 1: length(total_spikes(1,1,:,1))
% % %             for j = 1:length(total_spikes(1,1,1,:))
% % %                 temp(i,j) = total_spikes(1,1,i,j);
% % %             end
% % %         end
% % %         total_spikes = temp;
% % %     end
    total_spikes = squeeze(total_spikes);
    
    trial_number=0;
    trial_period=zeros(length(orientation),num_trials); %was (16,10)
    
    spike_rate=zeros(length(orientation),num_trials); %was (16,10)
    unave_true_spike_rate = zeros(length(orientation),1);
    true_spike_rate=zeros(length(orientation),1); % this was here first!
    norm_true_spike_rate=zeros(length(orientation),1); %
    
    x=zeros(length(orientation),1);
    y=zeros(length(orientation),1);
    
    error = zeros(length(orientation),1);
    unave_error = zeros(length(orientation),1);
    norm_error = zeros(length(orientation),1);
    
    for i=1:length(orientation); %%% can i do this with allspikes?!
        
        for j=1:length(all_orientations);
            if all_orientations(j)==orientation(i);
                trial_number = trial_number+1;
                trial_period(i,trial_number)=(trial_period(i,trial_number) + end_times(j)-start_times(j)) / 20000;
                spike_rate(i,trial_number)=total_spikes(i,trial_number)/trial_period(i,trial_number);
                % spikes per second for a given trial for a given
                % orientation 
            end
        end
        
        
        % Here we account for the real number of stimulus display epochs during
        % the recording. If there were some display epochs missed we will still
        % calculate a true average.
        for l = 1:trial_number;
            % l = 3:trial_number;  %%%%%THIS IS FOR THE ANOMALOUS RUN!!!
            % l = 1:trial_number;   %%%%% THIS IS WHAT YOU NORMALLY USE
            unave_true_spike_rate(i)=unave_true_spike_rate(i)+spike_rate(i,l);            
        end
        
        %%%%%% FOR NORMAL ANALYSIS
        true_spike_rate(i)=unave_true_spike_rate(i)/trial_number; % ave spike rate for a given orientation
        % % % % %         spike_rate_all_orientation = true_spike_rate; % saving all true spike rates for later analysis
        
        
        
        
        % Need Max spikes also for OS analysis
        max_index = find(true_spike_rate == max(true_spike_rate));
        if length(max_index) >1
            max_index = max_index(1);
%             fprintf('more than one orientations with the highest spike rate for idlist index %d - just choosing the first one\n', q)
        end
        
        Max_Spikes_absolute = true_spike_rate(max_index)*50;
       
        %%%%% FOR WRONG ANALYSIS
        %            true_spike_rate(i)=unave_true_spike_rate(i)/2; % average spike rate for a given orientation
        %           true_spike_rate(i)=unave_true_spike_rate(i)/3; % average spike rate for a given orientation
        
        
        %create matrix recording # of display epochs during an experiment for
        %each orientation
        x(i)=x(i)+trial_number;
        
        %calculate  std error (on the mean!!)
        unave_error(i)=std(spike_rate(i,1:trial_number)');
        error(i)=unave_error(i)/sqrt(trial_number);
        
        trial_number=0;
    end
 
    % (At this point we could plot polar graph of neuron - this has been moved
    % to stand_alone_printed_figures.m)
    
    %% Getting the normalized spike
    %For each orientation, divide by the sum of the rates of all
    % orientations, and divide errorbars by the same(?)
    
    sum_of_rates = 0;
    for i = 1:length(orientation);
        sum_of_rates = sum_of_rates + unave_true_spike_rate(i);
    end
    
    
    for i = 1:length(orientation);
        norm_true_spike_rate(i) = unave_true_spike_rate(i) / sum_of_rates;
        norm_error(i) = unave_error(i) / sum_of_rates;
    end