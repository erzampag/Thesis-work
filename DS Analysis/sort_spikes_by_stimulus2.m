function [total_spikes, combined_spikes, segments] = sort_spikes_by_stimulus2(orientation, temporal, spatial, frames, fourD_in,...
    num_trials, num_trials_NEW, spikeTimes, start_times, end_times, printfigures)


all_spikes = cell(length(spatial), length(temporal), length(orientation), num_trials); % all spikes in secs timed from start
graph_spikes = cell(length(spatial), length(temporal), length(orientation), num_trials); % same all all_spikes, but with edges for graphing
total_spikes = zeros(length(spatial), length(temporal), length(orientation), num_trials); % total spikes

combined_spikes = cell(length(spatial), length(temporal), length(orientation));

%   keyboard


% ALL_SPIKES: place spikes into 4D array based on start and end times
% TOTAL_SPIKES: the number of spikes in each run
for i = 1 : length(fourD_in)
    
    spike_indices = spikeTimes(:) > start_times(i) & spikeTimes(:) < end_times(i);
    
    all_spikes{fourD_in(i,1), fourD_in(i,2), fourD_in(i,3), fourD_in(i, 4)} ...
        = (spikeTimes(spike_indices)-start_times(i))/20000;
    
    total_spikes(fourD_in(i,1), fourD_in(i,2), fourD_in(i,3), fourD_in(i, 4)) = ...
        sum(length(all_spikes{fourD_in(i,1), fourD_in(i,2), fourD_in(i,3), fourD_in(i, 4)}));
    
end


% GRAPH_SPIKES: the spikes for making raster plots (includes edges)
for i  = 1: length(fourD_in)
    if i == 1 % first one is beginning to start
        
        spike_indices = spikeTimes(:) > 0 & spikeTimes(:) < start_times(i+1);
        graph_spikes{fourD_in(i,1), fourD_in(i,2), fourD_in(i,3), fourD_in(i, 4)} ...
            = (spikeTimes(spike_indices)-start_times(i))/20000;
        
    elseif i ~= 1 && i < length(fourD_in) % middle ones are end to start
        
        spike_indices = spikeTimes(:) > end_times(i - 1) & spikeTimes(:) < start_times(i+1);
        graph_spikes{fourD_in(i,1), fourD_in(i,2), fourD_in(i,3), fourD_in(i, 4)} ...
            = (spikeTimes(spike_indices)-start_times(i))/20000;
        
    elseif i == length(fourD_in) % last one is end to final
        
        spike_indices = spikeTimes(:) > end_times(i - 1) & spikeTimes(:) < spikeTimes(end) +1;
        graph_spikes{fourD_in(i,1), fourD_in(i,2), fourD_in(i,3), fourD_in(i, 4)} ...
            = (spikeTimes(spike_indices)-start_times(i))/20000;
        
    end
end

% keyboard
% COMBINED_SPIKES: histogram of spikes for each stimulus (combined over runs)

% % % % % % segments = .05:.1: frames/120 - 0.05; % 100 segments for 10 seconds, 50 segments for 5 seconds (each segment is 0.1 second)

% may have to create segments for each temporal period (ex, make each
% segment 1/4 of the temporal period
segments = cell(1,length(temporal));

for i = 1: length(temporal)
    segments{i} = temporal(i)/(16*120) : temporal(i)/(8*120) : frames/120 - temporal(i)/(16*120);
end

% segments = .25:.05: frames/120 - 0.25; % 100 segments for 10 seconds, 50 segments for 5 seconds (each segment is 0.1 second)


for i = 1 : length(spatial)
    for j = 1 : length(temporal)
        for k = 1: length(orientation)
            
% % %             combined_spikes_OLD{i,j,k} = 1/num_trials*(hist(cat(1, all_spikes{i,j,k,:}), segments));

            combined_spikes{i,j,k} = 1/num_trials_NEW(i,j,k)*(hist(cat(1, all_spikes{i,j,k,:}), segments{j}));
            
        end
    end
end

%   keyboard

if printfigures == 1
    % GRAPHING RASTER PLOTS
    for i = 1 : length(spatial)
        for j = 1 : length(temporal)
            figure % a figure for each spatial and temporal period
            for k = 1 : length(orientation)
                % subplot(4,4,k)
                axes('position', [(mod((k-1),4)/4)+.05, ((4-ceil((k*4)/16))/4)+.05, .15, .15 ]); % a subplot for each orientation
                axis([-3,frames/120 + 3,0,num_trials]);
                title(orientation(k),'FontSize',9,'FontWeight','bold');
                grid ON; hold on;
                
                for m = 1: num_trials % for each trial, plot the spikes
                    
                    for n = 1: length(graph_spikes{i,j,k,m}) % loop over spikes in one trial
                        x = [graph_spikes{i,j,k,m}(n), graph_spikes{i,j,k,m}(n)];
                        y = [m, m-1];
                        line(x,y);
                    end
                end
            end
            mtit(['Spatial: ', num2str(spatial(i)), '  Temporal: ' num2str(temporal(j))]) % after the figure was made
        end
    end
    
    % GRAPHING HISTOGRAMS
    for i = 1 : length(spatial)
        for j = 1 : length(temporal)
            figure % a figure for each spatial and temporal period
            for k = 1 : length(orientation)
                %                subplot(4,4,k)
                axes('position', [(mod((k-1),4)/4)+.05, ((4-ceil((k*4)/16))/4)+.05, .15, .15 ]); % a subplot for each orientation
                title(orientation(k),'FontSize',9,'FontWeight','bold');
                hold on;
                
                % % % %                 plot(segments, combined_spikes{i,j,k})
                plot(segments{j}, combined_spikes{i,j,k})
                
            end
            mtit(['Spatial: ', num2str(spatial(i)), '  Temporal: ' num2str(temporal(j))]) % after the figure was made
        end
    end
    
end



end



% % % % % % trial_counter = 1;
% % % % % % temp_spike_vector = [];
% % % % % % for i = 1 : length(spikeTimes)
% % % % % %
% % % % % %     if spikeTimes(i) > start_times(trial_counter) && spikeTimes(i) < end_times(trial_counter)
% % % % % %
% % % % % %         temp_spike_vector = [temp_spike_vector (spikeTimes(i)-start_times(trial_counter))/20000];
% % % % % %
% % % % % % %     elseif spikeTimes(i) >end_times(trial_counter) && spikeTimes(i) < start_times(trial_counter + 1)
% % % % % % %         continue
% % % % % %     else
% % % % % %         all_spikes{fourD_in(trial_counter,1), fourD_in(trial_counter,2), fourD_in(trial_counter,3), fourD_in(trial_counter, 4)} ...
% % % % % %             = temp_spike_vector;
% % % % % %         trial_counter = trial_counter +1;
% % % % % %
% % % % % %         temp_spike_vector = [];
% % % % % %         temp_spike_vector = [temp_spike_vector (spikeTimes(i)-start_times(trial_counter))/20000];
% % % % % %
% % % % % %
% % % % % %
% % % % % %     end
% % % % % %
% % % % % % end
% % % % % %
% % % % % %
% % % % % % trial_counter = 1;
% % % % % % for i = 1 : length(spikeTimes)
% % % % % %     if spikeTimes(i) > start_times(trial_counter) && spikeTimes(i) < end_times(trial_counter)
% % % % % %
% % % % % %         all_spikes{fourD_in(trial_counter,1), fourD_in(trial_counter,2), fourD_in(trial_counter,3), fourD_in(trial_counter, 4)} =...
% % % % % %             [all_spikes{fourD_in(trial_counter,:)} (spikeTimes(i)-start_times(trial_counter))/20000];
% % % % % %     else
% % % % % %         trial_counter = trial_counter +1;
% % % % % %
% % % % % %     end
% % % % % %
% % % % % % end
% % % % % %
% % % % % %
% % % % % % trial_counter = 1;
% % % % % % for i = 1 : length(spikeTimes)
% % % % % %     temp_spike_vector = [];
% % % % % %     if spikeTimes(i) > start_times(trial_counter) && spikeTimes(i) < end_times(trial_counter)
% % % % % %
% % % % % %         temp_spike_vector = [temp_spike_vector (spikeTimes(i)-start_times(trial_counter))/20000];
% % % % % %
% % % % % %     elseif spikeTimes(i) >end_times(trial_counter) && spikeTimes(i) < start_times(trial_counter + 1)
% % % % % %         continue
% % % % % %     else
% % % % % %         all_spikes{fourD_in(trial_counter,:)} = temp_spike_vector;
% % % % % %         trial_counter = trial_counter +1;
% % % % % %
% % % % % %     end
% % % % % %
% % % % % % end
% % % % % %
% % % % % %
% % % % % %
% % % % % % for i = 1: length(all_spatial_periods)
% % % % % %     for j = 1: length(all_temporal_periods)
% % % % % %         for k = 1:length(all_orientations)
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %
% % % % % %         end
% % % % % %     end
% % % % % % end


% % % temp_obj = all_spikes{find(spatial(:)==all_spatial_periods(1)),find(temporal(:)==all_temporal_periods(1)),...
% % % find(orientation(:)==all_orientations(1))}(1,:);
% % % append(temp_obj, 120)
% % %
% % %
% % % unique_runs = [100, 200]
% % % append(unique_runs, 300)
% % % concat(unique_runs, 300)
% % % unique_runs(1,3) = 300

