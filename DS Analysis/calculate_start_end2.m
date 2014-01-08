function [start_times, end_times, fourD_in, num_trials_NEW] = calculate_start_end2(neuronFile, frames, all_spatial_periods, all_temporal_periods,...
    all_orientations, spatial, temporal, orientation, num_trials)

% create a matrix of start times and end times in sample #s corresponding
% to when stimuli were being shown based on TTL pulses
% add
temp=neuronFile.getTTLTimes();% create matrix of TTL pulse times
TTL=double(temp); % make TTL pulse matrix double

start_times = zeros(length(all_temporal_periods),1);
end_times = zeros(length(all_temporal_periods),1);

display_length = frames*20000/120; % frames * 1sec/120frames * 20,000samples/1sec

diff(1) = 0;
start_times(1) = TTL(1);
end_times(1) = start_times(1) + display_length;

%    keyboard

% calculate start and end times using beginnings
stim_counter = 1;
for i = 1 : length(TTL)
    
    testme = [];
    if i > 1 && TTL(i) < end_times(stim_counter)
        testme = TTL(i) - TTL(i-1);
    end
%     abs(testme - all_temporal_periods(stim_counter)*20000/120)
    if abs(testme - all_temporal_periods(stim_counter)*20000/120) > 300 % 0.015 sec off (300 samples)
        abs(testme - all_temporal_periods(stim_counter)*20000/120)
        fprintf('this temporal period doesn''t match the one in the textfile: TTL %d\n', i)
    end

    if TTL(i) > end_times(stim_counter)
        
        stim_counter = stim_counter + 1;
        
        start_times(stim_counter) = TTL(i); % next start time is this TTL
        end_times(stim_counter) = start_times(stim_counter) + display_length;% next end time is display-length after
        
        if TTL(i) - TTL(i-1) > 72667; % 3.63 seconds - longer than 2.1sec which is longest period,
            fprintf('there hasn''t been a ttl pulse in %d seconds\n', (TTL(i)-TTL(i-1))/20000)
        end
        
    end
    
end
% keyboard
% % % % % DIAGNOSTICS
% % % % diff_frames = zeros(length(TTL), 1);
% % % % for i = 2:length(TTL)
% % % %    diff_frames(i-1) = (TTL(i) - TTL(i-1))*120/20000;
% % % % end
% % % % 



% Create make 4D indices here

fourD_in = zeros(length(all_orientations), 4);

for i = 1 : length(all_orientations)    
    fourD_in(i, 1) = find(spatial(:) == all_spatial_periods(i));
    fourD_in(i, 2) = find(temporal(:) == all_temporal_periods(i));
    fourD_in(i, 3) = find(orientation(:) == all_orientations(i));
end

unique_runs = unique(fourD_in, 'rows');

for i = 1 : length(unique_runs)
    these_trials = fourD_in(:,1) == unique_runs(i,1) & fourD_in(:,2) == unique_runs(i,2) & fourD_in(:,3) == unique_runs(i,3);
    fourD_in(these_trials,4) = 1:num_trials;
end

% % % % keyboard

% Have to deal with the fact that if some trials were skipped, the stimuli
% and the end were not shown and thus will have fewer trials

% start_times = start_times(find(start_times,1,'first'): find(start_times, 1,'last')); %% this doesn't work if 1st start time = 0
start_times = start_times(1: find(start_times, 1,'last'));

end_times = end_times(find(end_times,1,'first'): find(end_times, 1,'last'));

fourD_in = fourD_in(1:length(start_times),:);

num_trials_NEW = zeros(length(spatial), length(temporal), length(orientation));


for i = 1:length(fourD_in)
    num_trials_NEW(fourD_in(i,1),fourD_in(i,2),fourD_in(i,3)) = num_trials_NEW(fourD_in(i,1),fourD_in(i,2),fourD_in(i,3)) +1;
end


end


% % % % % calcualte start times using differences
% % % % stim_counter = 1;
% % % % testme =[];
% % % % testme2 = [];
% % % % 
% % % % testall = [];
% % % % testall2 = [];
% % % % for i = 2: length(TTL)
% % % %     diff(i) = TTL(i) - TTL(i-1);
% % % %         
% % % %     if abs(diff(i) - diff(i - 1)) > 100 && abs(diff(i) - temporal_periods(stim_counter)*20000/120) > 10 ...
% % % %             && TTL(i) > start_times(stim_counter) + display_length
% % % % %          keyboard
% % % %         testme = [testme i];
% % % %         
% % % %         testme2 = [testme2 TTL(i) - (start_times(stim_counter) + display_length)];
% % % %         stim_counter = stim_counter + 1;
% % % %         start_times(stim_counter) = TTL(i);
% % % %         
% % % %         
% % % %         if TTL(i) - TTL(i-1) > 72667;
% % % %             i
% % % %             fprintf('there hasnt been a ttl pulse in %d seconds, trust nobody, ...
% % % %                         and calculate start and ends independently\n', (TTL(i)-TTL(i-1))/20000)
% % % %             break
% % % %         end
% % % %         
% % % %         
% % % %     end
% % % %     
% % % % end
% % % % 
% % % % 
% % % % % end_times = start_times + display_length;
% % % % 
% % % % %% calculate start and end times assuming that lengths of movies and grey screen times are always correct (last ditch effort)
% % % % 
% % % % for i = 2 : length(temporal_periods)
% % % %     
% % % %     start_times(i) = end_times(i-1) + 30000;% this is for 1.5 seconds grey screen - 
% % % %     end_times(i) = start_times(i)  + display_length;
% % % % 
% % % % 
% % % % end




% 
% for i = 2:length(diff)
%     test(i) = diff(i) - diff(i-1);
%     test2(i) = diff(i) - temporal_periods(stim_counter)*20000/120;
% end

%  figure
%  plot(diff)