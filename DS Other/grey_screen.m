function [idList, spike_rate_all_orientation, max_true_spike_rate, min_true_spike_rate, contrast_spike_rate, ...
    grey_spike_rate] = grey_screen(date, datafile, stimfile2)

%% Define path to raw data
% Edit this function to access Vision.jar and raw data in a different location
[pfile, neuronFile] = define_path_to_raw_data(date, datafile);

%% Determine how many figures to print:
% 0 - print nothing but SAVE VARIABLES!,
% 1 - print everything,
% 2 - print PSTH, FFT, tuningcurve, vector for neuron,
% 3 - print tuning curve and vector together,
% 4 - just final whole retina analysis
% 5 - von Mises fit and data
% 6 - tuning curve FFT
% 7 - OS von Mises fit
printfigures = 1;

%% Create directories to save figures for given cell types found

A= strcat('/Users/erinzampaglione/Documents/Lab_Work/SbC_cells/', date, '/', datafile);

% A= strcat('/Users/erinzampaglione/Documents/processed_data/datafromextHD/', date, '/', datafile); % for runthru

mkdir(A);
% B= strcat('/Users/erinzampaglione/Documents/Second term 1styear/DSCells', date, '/', datafile, '/OS');
% mkdir(B);
% C= strcat('/Users/erinzampaglione/Documents/Second term 1styear/DSCells',date, '/', datafile, '/DS');
% mkdir(C);
% D= strcat('/Users/erinzampaglione/Documents/Second term 1styear/DSCells', date, '/', datafile, '/NonOS');
% mkdir(D);


%% Get neuron IDs from Vision (before duplication removal)
%
idListPreD =  import_neuronIDs(neuronFile);

%% Get information from parameters file (classes, xOffDS, etc)
%
[data, EIx, EIy] = import_parameters(pfile, idListPreD);

classes = import_classes(pfile, idListPreD);
%% Remove Duplicates in IDlist and any other parameter
%
[idList, ind] = remove_duplicates(idListPreD, classes);

EIx = EIx(ind);
EIy = EIy(ind);
data = data(ind);

%% Define variables for use in analysis


% AllCells= zeros(length(idList), 2);   %these arrays will contain 2 numbers for each DS or OS cell: direction and magnitude.

% ratio = zeros(1,length(idList));
% noise = zeros(1,length(idList));
% noiseratio = zeros(1,length(idList));
% zerofreq = zeros(1,length(idList));
%
% DSI_error = zeros(1,length(idList));
%
% fwhh = zeros(1,length(idList));
% % size_confint = [];
% % fwhh_fractional_error = [];
% reduced_chi_square = zeros(1,length(idList));
% % reduced_chi_square_binom = [];
%
% MaxSpikes = zeros(1,length(idList));
%
% Fone_Fzero = zeros(1,length(idList));
% Ftwo_Fzero = zeros(1,length(idList));
% OS_x1 = zeros(1,length(idList)); OS_y1 = zeros(1,length(idList));
% OS_x2 = zeros(1,length(idList)); OS_y2 = zeros(1,length(idList));
% OS_counter = 0;
% Max_Spikes_absolute = zeros(1,length(idList));
%
spike_rate_all_orientation = {};

max_true_spike_rate = zeros(1, length(idList));
min_true_spike_rate = zeros(1, length(idList));
contrast_spike_rate = zeros(1, length(idList));
grey_spike_rate = zeros(1, length(idList));

%% %load stimuli information from a list of orientations displayed for a
% %particular experiment in proper order from a .txt file in Matlab directory

% text file from LISP

% % % % % % % % % % % % % % % stim_file = fopen(strcat('/Users/erinzampaglione/Documents/processed_data/', date, '/stimuli/',stimfile2)); % for runthru

stim_file = fopen(strcat('/Users/erinzampaglione/Documents/processed_data/datafromextHD/', date, '/stimuli/',stimfile2)); % for runthru


% stim_file = fopen(strcat('/Volumes/TEMP_MOUSE/', '2010-08-27-0', '/stimuli/', 's07.txt'));
stimFile = textscan(stim_file, '%s%d%s%d%s%f%s', 'headerlines', 1);
stim = stimFile{6};
% spatial_periods = stimFile{2};
% temporal_periods = stimFile{4};


% text file from PTB
%%%%stim_file = fopen('/Users/erinzampaglione/Desktop/TEST/gratings_sets.txt')
%%%%stimFile = textscan(stim_file, '%d%d%d%d%d%d%d', 'headerlines', 1);
%%%%stim = stimFile{7}
%%%%spatial_periods =
%%%%temporal_periods =
%% Use TTL pulses to determine the start and end times of the stimuli
%
[start_times, End_times] = calculate_start_end(neuronFile);

%% Start of the loop over neurons

%for given neuron, import the spikes times, just doing one right now so
%there's not too many plots popping up.

%2011-02-23-0
TROUBLE = [36; 146;168;182;186;220;234;247;252;274;306;398;471;475;488;514]; % grey_spike_rate > 20 & grey_spike_rate > contrast_spike_rate
TROUBLE = [36;46;50;146;168;182;186;220;234;247;252;274;306;391;398;421;471;475;488;507;514];%grey_spike_rate > 20 & grey_spike_rate > min_true_spike_rate
TROUBLE = [7;105;128;148;175;199;205;209;212;224;230;233;267;433;435;436;496;501;508]; %grey_spike_rate(:) > 5 & contrast_spike_rate(:)< 5

%2011-02-16-0
% TROUBLE = [1;18;45;90;114;192;193;397;450;518;524;545;580;607];%grey_spike_rate(:) > 5 & contrast_spike_rate(:)< 5
%
fprintf('Analyzing Retina from %s \n', date)
% for q=1:length(idList);
for q = 1: length(TROUBLE);
    
    if mod(q,50)==0,
        fprintf('processing neuron # %d \n',q);
        
    end
    
%     NeuronID=idList(q);
    
            NeuronID = idList(TROUBLE(q)); % if TROUBLE is the indices
    %
    %            NeuronID = TROUBLE(q); % if TROUBLE is the neuron IDs
    
    temp=neuronFile.getSpikeTimes(NeuronID);
    spikeTimes=double(temp);%converts the int matrix "spikeTimes" to a double matrix
    
    
    orientation=0:22.5:337.5;
    
    %% Creation of the raster plots, frequency (PSTH) plots for each orientation, in new figure window
    % Output PSTH for FFT and total spikes for finding rate / dir pref
    %
    [total_spikes, total_spikes_grey, combinedspikes, segments] = sort_spikes_by_stimulus_grey(orientation, stim, spikeTimes,...
        start_times, End_times, printfigures);
    
    %% Spike rate for grey screens
    
    all_spikes_grey = sum(sum(total_spikes_grey));
    grey_spike_rate(q) = all_spikes_grey / (3*80);
    
    
    
    %% Average spike rates
    % calculate average spike rate of cell during the display of each of 16 orientations of a moving bar
    trial_number=0;
    trial_period=zeros(16,5); %was (16,10)
    
    spike_rate=zeros(16,5); %was (16,10)
    unave_true_spike_rate = zeros(16,1);
    true_spike_rate=zeros(16,1); % this was here first!
    norm_true_spike_rate=zeros(16,1); %
    
    x=zeros(16,1);
    y=zeros(16,1);
    
    %     spikeSumx=0;  %to be used in adding up the vectors for different direction-selectivities
    %     spikeSumy=0;
    norm_spikeSumx = 0;
    norm_spikeSumy = 0;
    
    error = zeros(16,1);
    unave_error = zeros(16,1);
    norm_error = zeros(16,1);
    
    for i=1:16; %%% can i do this with allspikes?!
        
        for j=1:length(stim);
            if stim(j)==orientation(i);
                trial_number = trial_number+1;
                trial_period(i,trial_number)=(trial_period(i,trial_number)...
                    +End_times(j)-start_times(j))/20000;
                spike_rate(i,trial_number)=total_spikes(i,trial_number)/trial_period(i,trial_number);
                % spikes per second for a given trial for a given
                % orientation
                
            end
            
        end
        
        
        % Here we account for the real number of stimulus display epochs during
        % the recording. If there were some display epochs missed we will still
        % calculate a true average.
        for l = 1:trial_number;
            
            
            unave_true_spike_rate(i)=unave_true_spike_rate(i)+spike_rate(i,l);
            
            % probably total spikes could go here?
            
        end
        %%%%%% FOR NORMAL ANALYSIS
        true_spike_rate(i)=unave_true_spike_rate(i)/trial_number; % ave spike rate for a given orientation
        
        x(i)=x(i)+trial_number;
        
        %calculate  std error (on the mean!!)
        unave_error(i)=std(spike_rate(i,1:trial_number)');
        error(i)=unave_error(i)/sqrt(trial_number);
        
        trial_number=0;
    end
    
    
    spike_rate_all_orientation{q} = true_spike_rate;
    
    max_true_spike_rate(q) = max(true_spike_rate);
    min_true_spike_rate(q) = min(true_spike_rate);
    contrast_spike_rate(q) = sum(true_spike_rate)/16;
    
    % (At this point we could plot polar graph of neuron - this has been moved
    % to stand_alone_printed_figures.m)
    %% Getting the normalized spike
    %For each orientation, divide by the sum of the rates of all
    % orientations, and divide errorbars by the same(?)
    
    sum_of_rates = 0;
    for i = 1:16;
        sum_of_rates = sum_of_rates + unave_true_spike_rate(i);
    end
    
    
    for i = 1:16;
        norm_true_spike_rate(i) = unave_true_spike_rate(i) / sum_of_rates;
        norm_error(i) = unave_error(i) / sum_of_rates;
    end
    
    % Then I want normalized x and y vectors to get the coordinates for DS
    % vectors
    
    for j = 1: 16;
        norm_spikeSumx = norm_spikeSumx + norm_true_spike_rate(j)*(cos(orientation(j)*(2*pi/360)));
        norm_spikeSumy = norm_spikeSumy + norm_true_spike_rate(j)*(sin(orientation(j)*(2*pi/360)));
    end
    %         keyboard
    
    
    if printfigures == 1
        figure
        h= polar([0:0.1:2*pi], grey_spike_rate(q)*ones(1,length([0:0.1:2*pi])), '-k');
        set(h,'LineWidth',1.5)
        hold on
        h = polar(orientation(:)*(pi/180), true_spike_rate(:), '-b'); % data points
        set(h,'LineWidth',1.5)
        title(strcat('Tuning and grey spike rate for neuron-', num2str(NeuronID)));

        %     set(h, 'MarkerFaceColor', 'r');
    end
    
    
    
    
end

%     keyboard
figure
scatter(contrast_spike_rate, grey_spike_rate);
xlabel('Average spike rate for all orientations','FontSize', 25); ylabel('Average spike rate at grey screen','FontSize', 25);
hold on
plot((0:0.1:50),(0:0.1:50),'k-'); plot((0:0.1:50), 5*ones(1,length(0:0.1:50)), 'k-'); plot(5*ones(1,length(0:0.1:50)),(0:0.1:50), 'k-');
% % % % % % % % % %     figure
% % % % % % % % % %     scatter(max_true_spike_rate, grey_spike_rate);
% % % % % % % % % %     xlabel('Highest spike rate of all orientations','FontSize', 25); ylabel('Average spike rate at grey screen','FontSize', 25);
% % % % % % % % % %     figure
% % % % % % % % % %     scatter(min_true_spike_rate, grey_spike_rate);
% % % % % % % % % %     xlabel('Lowest spike rate of all orientations','FontSize', 25); ylabel('Average spike rate at grey screen','FontSize', 25);
% % % % % % % % % %



if printfigures == 0;
    %     filename = strcat(A, '_', datestr(now, 30)); % maybe you don't want this datestr for now
    filename = A;
    save(filename, 'idList', 'spike_rate_all_orientation', 'max_true_spike_rate', 'min_true_spike_rate', 'contrast_spike_rate', ...
        'grey_spike_rate');
    
    
end

%     keyboard