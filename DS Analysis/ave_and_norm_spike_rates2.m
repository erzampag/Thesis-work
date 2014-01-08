function [true_spike_rate, norm_true_spike_rate, Max_Spikes_absolute, max_index] =...
    ave_and_norm_spike_rates2(spatial, temporal, orientation, num_trials, num_trials_NEW, frames, total_spikes, printfigures)


spike_rate_each_trial = zeros(length(spatial), length(temporal), length(orientation), num_trials);
true_spike_rate = zeros(length(spatial), length(temporal), length(orientation));
norm_true_spike_rate = zeros(length(spatial), length(temporal), length(orientation));

Max_Spikes_absolute = zeros(length(spatial),length(temporal));
max_index = zeros(length(spatial),length(temporal));

% % % keyboard

% Get all spikes rates
for i  = 1: length(spatial)
    for j = 1 : length(temporal)
        for k = 1: length(orientation)
            for m = 1: num_trials
                spike_rate_each_trial(i,j,k,m) = total_spikes(i,j,k,m)/(frames/120); % the spike rate for each stimuls display
            end
            true_spike_rate_OLD(i,j,k) = sum(total_spikes(i,j,k,:))/(frames*num_trials/120); % the spike rate for each spatial/temporal/orientation
            true_spike_rate(i,j,k) = sum(total_spikes(i,j,k,:))/(frames*num_trials_NEW(i,j,k)/120); % spike rate for each sp/temp/orientation
            
            spike_rate_per_period(i,j,k) = sum(total_spikes(i,j,k,:))*temporal(j)/(frames*num_trials_NEW(i,j,k)); %
        end
        
        %%%%% need to account for zero spikes - normalization can't be done
        %%%%% by dividing by the sum because you get an indeterminant
        if sum(true_spike_rate(i,j,:)) ~= 0
            norm_true_spike_rate(i,j,:) = true_spike_rate(i,j,:)/sum(true_spike_rate(i,j,:)); % normalized rate for each spatial/temporal/orientation
        else
            norm_true_spike_rate(i,j,:) = true_spike_rate(i,j,:);
            fprintf('look''s like there''s zero spikes for all orientations of this sp/temp, normalized spike rate set to zero\n')
            
        end
    end
end

keyboard

% Calculate error on the means for the true and norm spike rates assuming
% N distibuted as a Poisson process

for i  = 1: length(spatial)
    for j = 1 : length(temporal)
        for k = 1: length(orientation)
            for m = 1: num_trials
            end
            error(i,j,k) = sqrt(sum(total_spikes(i,j,k,:)))/(frames*num_trials_NEW(i,j,k)/120); % spike rate for each spatial/temporal/orientation
        end
        norm_error(i,j,:) = error(i,j,:)/sum(true_spike_rate(i,j,:)); % normalized rate for each spatial/temporal/orientation
    end
end

%%% old error assuming trials are a gaussian process
% % % % % % for i  = 1: length(spatial)
% % % % % %     for j = 1 : length(temporal)
% % % % % %         for k = 1: length(orientation)
% % % % % %
% % % % % %
% % % % % %             error(i,j,k) = std(spike_rate_each_trial(i,j,k,:)) / sqrt(num_trials);
% % % % % %
% % % % % %
% % % % % %         end
% % % % % %
% % % % % %         sum_of_rates = sum(norm_true_
% % % % % %
% % % % % %         norm_error = std(spike_rate_each_trial(i,j,k,:)/sum(true_spike_rate(i,j,:)) / sqrt(num_trials);
% % % % % %
% % % % % %     end
% % % % % % end



% Determine the maximum spikes and the index of the maximum spikes
for i = 1: length(spatial)
    for j = 1: length(temporal)
        %         if max(true_spike_rate(i,j,:)) > Max_Spikes_absolute(i,j)
        
        Max_Spikes_absolute(i,j) = max(true_spike_rate(i,j,:));
        temp_index = find(true_spike_rate(i,j,:) == Max_Spikes_absolute(i,j));
        max_index(i,j) = temp_index(1); % for if there are more than one orientations where the spike rate is highest.
        
        Max_Spikes_absolute(i,j) = max(true_spike_rate(i,j,:))*(frames*num_trials_NEW(i,j,max_index(i,j))/120);
  
        %         end
    end
end

% Max_Spikes_absolute(i,j) = max(true_spike_rate(i,j,:))*(5*5);

% Max_Spikes_absolute(i,j) = max(true_spike_rate(i,j,:))*(frames*num_trials_NEW(i,j,temp_index)/120);

% %  keyboard
% print N figures that show polar plot AND plot with poisson error bars
if printfigures == 1
    for i  = 1:length(spatial)
        for j = 1:length(temporal)
            figure
            subplot(2,1,1)
            polar([orientation(:)*pi/180'; orientation(1)], [squeeze(true_spike_rate(i,j,:)); true_spike_rate(i,j,1)]);
            title(['Spatial period: ', num2str(spatial(i)), '  Temporal period: ' num2str(temporal(j))])
            subplot(2,1,2)
            plot(orientation(:), squeeze(true_spike_rate(i,j,:)));
            errorbar(orientation(:), true_spike_rate(i,j,:), squeeze(error(i,j,:)));
            
        end
    end
end

% heat map for parameter scan

if printfigures == 8
    
    spat_temp_rates = zeros(length(spatial), length(temporal));
    
    for i  = 1:length(spatial)
        for j = 1:length(temporal)
            spat_temp_rates(i,j) = sum(true_spike_rate(i,j,:))/length(true_spike_rate(i,j,:));
        end
    end
    
    %     for i = 1 : length(spatial)
    %         for j = 1:length(orientation)
    
    
    % let's put these here
    
    
    figure
    imagesc([spatial_conversion(spatial(1)) spatial_conversion(spatial(5))],...
        [temporal_conversion(temporal(1)) temporal_conversion(temporal(5))], spat_temp_rates)
    
    imagesc([spatial(1) spatial(5)], [temporal(1) temporal(5)], spat_temp_rates)
    
    ylabel('Spatial periods'); xlabel('Temporal Periods'); title('Average Spike Rate over all Orientations');
    colormap(gray); colorbar
    
    
    
    % attempt at surface plot
    
    figure
    %     subplot(1,3,1)
    [X,Y] = meshgrid(log(spatial), log(temporal));
    surf(X, Y, spat_temp_rates);
    colormap(gray); colorbar
    
    %     subplot(1,3,2)
    %     [X,Y] = meshgrid(log(spatial), );
    %     surf(X, Y, spat_temp_rates);
    %     colormap(gray); colorbar
    %
    %     subplot(1,3,3)
    %     [X,Y] = meshgrid(log(spatial), log(temporal));
    %     surf(X, Y, spat_temp_rates);
    %     colormap(gray); colorbar
    
end

end


function spatial_conv = spatial_conversion(spatial)
% converts pixels to cyc/mm
spatial_conv = 1/(double(spatial)*9); % 9um = 1 pixel, 31um = 1 degree


end

function temporal_conv = temporal_conversion(temporal)
% converts frames to cyc/sec
temporal_conv = 1/(double(temporal)/120); % 120 frames = 1 sec

end
