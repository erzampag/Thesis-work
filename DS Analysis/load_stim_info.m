function [stim, spatial_periods, temporal_periods, orientation, spatial, temporal, num_trials, frames] = ...
    load_stim_info(date, stimfile2)

% Load stimuli information from a list of orientations displayed for a particular experiments
% in proper order from a .txt file in Matlab directory

% text file from LISP

stim_file = fopen(strcat('/Users/erinzampaglione/Documents/processed_data/', date, '/stimuli/',stimfile2)); % for runthru
% stim_file = fopen(strcat('/Users/erinzampaglione/processed_data/', date, '/stimuli/',stimfile2)); % for runthru
% keyboard
line1 = fgets(stim_file);

if length(line1) > 100 % Drifting squarwave stimulus
    frames =sscanf(line1, '(:TYPE :DRIFTING-SQUAREWAVE :RGB #(0.48 0.48 0.48) :X-START 0 :X-END 640 :Y-START 0 :Y-END 480 :FRAMES %d)');
%     stimFile = textscan(stim_file, '%s%d%s%d%s%f%s', 'headerlines', 1);
    stimFile = textscan(stim_file, '%s%f%s%f%s%f%s'); % for some reason 
else
    frames =sscanf(line1, '(:TYPE :DRIFTING-SINUSOID :X-START 0 :X-END 640 :Y-START 0 :Y-END 480 :FRAMES %d)');
    stimFile = textscan(stim_file, '%s%f%s%f%s%f%*[^\r\n]'); %%% needs to be changed to incorporate RGB info!!!
end

stim = stimFile{6};

spatial_periods = stimFile{2};
temporal_periods = stimFile{4};


orientation = unique(stimFile{6})';
spatial = unique(stimFile{2})';
temporal = unique(stimFile{4}');

num_trials = length(stimFile{6})./(length(orientation)*length(spatial)*length(temporal));

fclose(stim_file);

% text file from PTB

%%%%stim_file = fopen('/Users/erinzampaglione/Desktop/TEST/gratings_sets.txt')
%%%%stimFile = textscan(stim_file, '%d%d%d%d%d%d%d', 'headerlines', 1);
%%%%stim = stimFile{7}
%%%%spatial_periods =
%%%%temporal_periods =
end