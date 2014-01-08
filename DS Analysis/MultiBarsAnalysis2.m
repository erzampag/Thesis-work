function [idList, AllCells, DSI, DSI_error, ratio, noiseratio,...
    noise, zerofreq, fwhh, reduced_chi_square, MaxSpikes, EIx, EIy,...
    Fone_Fzero, Ftwo_Fzero, OS_x1, OS_y1, OS_x2, OS_y2, Max_Spikes_absolute, spike_rate_all_orientation]=...
    MultiBarsAnalysis2(date, datafile, stimfile2)

% tic
% This function determines receptive field properties of neurons stimulated
% with a drifting square-wave stimulus presented in 16 different directions
% at five trials at each orientation.
% By determining the avererage spike rate at each orientation, a preferred orientation is found
%(direction closest to the vector sum of the average spike rates at each orientation).
% From this direction-selectivity index (DSI), width of tuning curve (fwhh),
% and EI position are determined.

%Erin Zampaglione, Summer 2013
%Adapted from drifting_squarewave_shorter_stimulus2, by Willie Tobin circa early 2010 and modified by Jason Triplett
%10-12-2010, Richard Smith in Fall 2011, and Erin Zampaglione in Winter
%2012

%% get neuron info
% 
[idListPreD, classes, ~, EIx, EIy, data, ~, neuronFile] = import_neuron_info(...
    '/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
    '/Users/erinzampaglione/Documents/processed_data/', date, datafile);

% % % [idListPreD, classes, ~, EIx, EIy, data, ~, neuronFile] = import_neuron_info(...
% % %     '/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
% % %     '/Users/erinzampaglione/processed_data/', date, datafile);

[idList, ind] = remove_duplicates(idListPreD, classes);
% keyboard
EIx = EIx(ind);
EIy = EIy(ind);
data = data(ind);

%% Determine how many figures to print:
% 0 - print nothing but SAVE VARIABLES!,
% 1 - print everything,
% 2 - print PSTH, FFT, tuningcurve, vector for neuron,
% 3 - print tuning curve and vector together,
% 4 - just final whole retina analysis
% 5 - von Mises fit and data
% 6 - tuning curve FFT
% 7 - OS von Mises fit
% 8 - test
printfigures = 0;

%% Create directories to save figures for given cell types found
B = strcat('/Users/erinzampaglione/Documents/Lab_Work/DS_param/', date);
A = strcat('/Users/erinzampaglione/Documents/Lab_Work/DS_param/', date, '/', datafile);
mkdir(B);

%% Define variables for use in analysis

AllCells = cell(length(idList), 2);   %these arrays will contain 2 numbers for each DS or OS cell: direction and magnitude.

spike_rate_all_orientation = cell(length(idList),1);

MaxSpikes = cell(length(idList), 1);
Max_Spikes_absolute = cell(length(idList),1);

DSI = cell(length(idList),1);
DSI_error = cell(length(idList),1);

ratio = cell(length(idList),1);
noise = cell(length(idList),1);
noiseratio = cell(length(idList),1);
zerofreq = cell(length(idList),1);

fwhh = cell(length(idList),1);
% size_confint = [];
% fwhh_fractional_error = [];
reduced_chi_square = cell(length(idList),1);
% reduced_chi_square_binom = [];

Fone_Fzero = cell(length(idList),1);
Ftwo_Fzero = cell(length(idList),1);
OS_x1 = cell(length(idList),1); OS_y1 = cell(length(idList),1);
OS_x2 = cell(length(idList),1); OS_y2 = cell(length(idList),1);

PSTHs = cell(length(idList), 1);
best_orientation = cell(length(idList), 1);


%% Load stimulus information

[all_orientations, all_spatial_periods, all_temporal_periods, orientation, spatial, temporal, num_trials, frames] =...
    load_stim_info(date, stimfile2);

for i  = 1: length(idList)
    AllCells{i,1} = zeros(length(spatial), length(temporal));
    AllCells{i,2} = zeros(length(spatial), length(temporal));
    
    MaxSpikes{i,1} = zeros(length(spatial), length(temporal));
    Max_Spikes_absolute{i,1} = zeros(length(spatial), length(temporal));
    
    DSI{i,1} = zeros(length(spatial), length(temporal));
    DSI_error{i,1} = zeros(length(spatial), length(temporal));
    
end
PSTH_bins = cell(1, length(temporal));
%% Use TTL pulses to determine the start and end times of the stimuli

% % % % % [start_times, end_times] = calculate_start_end(neuronFile);

[start_times, end_times, fourD_in, num_trials_NEW] = calculate_start_end2(neuronFile, frames, all_spatial_periods, all_temporal_periods,...
    all_orientations, spatial, temporal, orientation, num_trials);


%% Troubleshooting subsets of neurons using either their NeuronIDs or indices

% TROUBLE = [180;185;188;189;195;204;211;214;216;219;221;222;254;261;263;264;272;277;278;281;282;291;295]; %2012-02-23-PAPER FIGURE IN THIS GROUP!!

% TROUBLE = [20; 116; 128;138;144;153;155;175;178;191;194;195;203;214;227;228;239;244;251;260;267;271;285;287;377;443;444;448;...
%     482;516;528;530;532;]; % OS from FFT 2010-08-25
% TROUBLE = [6;16;19;33;42;47;48;56;57;61;64;65;68;69;71;75;86;87;89;90;94;99;101;110;111;112;129;130;142;151;157;162;163;170;171;...
%     180;185;188;189;195;197;204;210;211;214;216;219;221;229;254;261;262;263;264;272;277;278;281;282;291;295;304;305;308;313;317;318;...
%     319;324;334;335;347;348;349;353;355;357;358;360;368;372;373;376;380;395;404;414;418;420;458;484;494;516;528;536;]; % DS from FFT 2011-02-23
TROUBLE = [7;16;33;47;53;54;55;57;59;65;68;69;71;83;86;87;89;90;97;100;111;130;133;142;151;157;170;179;180;183;188;194;195;197;205;...
    211;213;223;237;254;263;264;272;280;281;282;294;295;304;305;307;317;319;324;334;347;348;349;353;354;355;357;358;368;370;372;373;377;...
    380;395;414;418;420;433;435;436;484;490;494;501;503;515;518;528;529;]; % OS from FFT 2011-02-23
%     TROUBLE = [25;33;66;74;85;94;98;105;109;118;122;124;130;140;147;165;170;177;187;217;219;235;259;273;323;325;331;332;334;338;363;372;378;...
%         386;394;402;412]; % OS from FFT 2011-06-28
%     TROUBLE = [32;45;52;61;62;68;84;88;92;112;168;200;223;224;249;254;260;...
%        277;307;308;324;330;351;390;392;407]; % OS from FFT 2011-03-25(DSCAM-/-)

% TROUBLE = [812;961;980;1277;3121;3393;3933;4187;4621;5088;6001;7667]; % 2011-08-24 neuron IDs of OS cells
% TROUBLE = [54 55];

TROUBLE = [2, 3];
% TROUBLE = 18;

%% Start of the loop over neurons

fprintf('Analyzing Retina from %s \n', date)
for q=1:length(idList);
% for q = 1: length(TROUBLE);
    
    if mod(q,50)==0,
        fprintf('processing neuron # %d \n',q); 
    end
    
    NeuronID=idList(q); % complete run-thru
%     
%     NeuronID = idList(TROUBLE(q)); % if TROUBLE is the indices
     
%     NeuronID = TROUBLE(q); % if TROUBLE is the neuron IDs
    
    temp=neuronFile.getSpikeTimes(NeuronID);
    spikeTimes=double(temp);%converts the int matrix "spikeTimes" to a double matrix
%     spikeTimes = spiketimes{NeuronID};

    if length(spikeTimes) < 100
        fprintf('Neuron ID %d had only %d spikes, and we''re gonna skip it\n', NeuronID, length(spikeTimes))
         continue
    end
    
    %% Creation of the raster plots, frequency (PSTH) plots for each orientation, in new figure window
    % Output PSTH for FFT and total spikes for finding rate / dir pref
    % If you want grey screen spike rates, use sort_spikes_by_stimulus_grey     

% %     [total_spikes, combined_spikes, segments] = sort_spikes_by_stimulus(orientation, all_orientations, spikeTimes,...
% %         start_times, end_times, printfigures);
    
    [total_spikes, combined_spikes, segments] = sort_spikes_by_stimulus2(orientation, temporal, spatial, frames, fourD_in, ...
        num_trials, num_trials_NEW, spikeTimes, start_times, end_times, printfigures);

    PSTHs{q} = combined_spikes;
    PSTH_bins = segments;
    
    %% Average and Normalize Spikes rates
%     keyboard
    
% %     [true_spike_rate, norm_true_spike_rate, Max_Spikes_absolute(q), max_index] =...
% %         ave_and_norm_spike_rates(orientation, all_orientations, num_trials, start_times, end_times, total_spikes);


    [true_spike_rate, norm_true_spike_rate, Max_Spikes_absolute{q}, max_index] =...
    ave_and_norm_spike_rates2(spatial, temporal, orientation, num_trials, num_trials_NEW, frames, total_spikes, printfigures);

    
    %%%neeed to make maxspikesabsolute to be a cell array
    spike_rate_all_orientation{q} = squeeze(true_spike_rate); % saving all true spike rates for later analysis
    


    
    %% Calculation of R and theta
    
    for i = 1 : length(spatial)
        for j = 1: length(temporal)
            [AllCells{q,2}(i,j), AllCells{q,1}(i,j)] = calculate_R_theta(squeeze(norm_true_spike_rate(i,j,:)), orientation); % R, THETA
        end
    end
    
    %% New Determination of preferred orientation! Also calculate error on DSI
    
    index_of_min = zeros(length(spatial), length(temporal));
    pref_O = zeros(length(spatial), length(temporal));
    
%     if q == 18
%         keyboard
%     end
    
    
    for i = 1 : length(spatial)
        for j = 1 : length(temporal)
            
            [DSI{q}(i,j), MaxSpikes{q}(i,j), index_of_min(i,j), ~, pref_O(i,j), DSI_error{q}(i,j)] =...
                determine_prefO_by_DSvector(AllCells{q,1}(i,j), orientation, squeeze(true_spike_rate(i,j,:)), q);
            
        end
    end
    
    if unique(pref_O) < 3
       fprintf('check out neuron %d', q); 
    end
   
%      keyboard % adapted thru here - do i want to keep these similar and loop over all spatial/temporals?

    %% Choosing the best stimulus to determine if it is ON/OFF or just ON or OFF
    
    for i  = 1 : length(spatial)
        for j = 1 : length(temporal)
                        
            if DSI{q}(i,j) < 0.5
                thebeststimulus = max_index(i,j); % anything that is not a DS cell
            else
                thebeststimulus = index_of_min(i,j); % for DS cells, we want to look at the direction closest to the vector sum
            end
            
            
            best_orientation{q}(i,j) = thebeststimulus;

            %             keyboard
            %     [ratio(q), noise(q), zerofreq(q), noiseratio(q)] = fourier_transform_PSTH(segments, thebeststimulus, combined_spikes, q);
            
            [ratio{q}(i,j), noise{q}(i,j), zerofreq{q}(i,j), noiseratio{q}(i,j)] = fourier_transform_PSTH2(segments{j}, frames,...
                thebeststimulus, orientation, squeeze(combined_spikes(i,j,:)), temporal(j), spatial(i), printfigures);
        end
    end
    
%     keyboard
    
%     noise{q}
%     noiseratio{q}
    
%     keyboard
        %% Fit tuning curve, print figure if ==5
        
        % I might only want to do this tuning curve fit if the cell is
        % identified as DS - either with DSI or AllCells(q,2)
    
% % % % % % % % %     [fwhh(q), reduced_chi_square(q)] = fit_tuning_curve(pref_Orate, pref_O, temp_theta_degrees, squeeze(true_spike_rate),...
% % % % % % % % %         orientation, NeuronID, q, printfigures, AllCells);
    
    % % %     keyboard
    
    %% calculation of DS and OS by FFT
    
% % % % % % % % %   
% % % % % % % % %     [Fone_Fzero(q), Ftwo_Fzero(q), OS_x1(q), OS_y1(q), OS_x2(q), OS_y2(q)] = ...
% % % % % % % % %         determine_prefO_by_FFT(squeeze(true_spike_rate),squeeze(norm_true_spike_rate),...
% % % % % % % % %         AllCells, pref_Orate, pref_O, orientation, printfigures, q, NeuronID, Max_Spikes_absolute);
% % % % % % % % %     
    
    
    %% Printing Figures!
    % FFT, single polar plot, etc - based on printed_figures switch
% % % %     if printfigures ~=0
% % % %         stand_alone_printed_figures(printfigures, combined_spikes, segments,...
% % % %             orientation, true_spike_rate, NeuronID, AllCells, index_of_min, q, Fone_Fzero, Ftwo_Fzero, power_tuning, thebeststimulus);
% % % %     end
%        keyboard
end
% toc % end of run through all neurons in retina
% keyboard
%% START OF WHOLE RETINA ANALYSIS

if printfigures == 4
    full_retina_printed_figures(Fone_Fzero, Ftwo_Fzero, AllCells, idList, DSI, ratio, noiseratio,...
    MaxSpikes, Max_Spikes_absolute, OS_x1, OS_y1, OS_x2, OS_y2)
end

% % % % % % % % 
% % % % % % % % DScellsindex = [];
% % % % % % % % DScellsindex = find(ratio(:) >1 & DSI(:) > 0.5 & noiseratio(:) >5);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % keyboard
if printfigures == 0;
    %     filename = strcat(A, '_', datestr(now, 30)); % maybe you don't want this datestr for now
    filename = A;
    
        save(filename, 'idList', 'AllCells', 'DSI','DSI_error', 'ratio',...
        'noiseratio', 'noise', 'zerofreq', ...
        'MaxSpikes', 'EIx', 'EIy',...
        'Max_Spikes_absolute', 'spike_rate_all_orientation',...
        'spatial', 'temporal', 'orientation', 'PSTHs', 'PSTH_bins', 'best_orientation', 'date');
    
    
% %     save(filename, 'idList', 'AllCells', 'DSI', 'OSI','DSI_error', 'ratio',...
% %         'noiseratio', 'noise', 'zerofreq', 'fwhh',...
% %         'reduced_chi_square', 'MaxSpikes', 'EIx', 'EIy',...
% %         'Fone_Fzero','Ftwo_Fzero', 'OS_x1', 'OS_y1', 'OS_x2', 'OS_y2', 'Max_Spikes_absolute', 'spike_rate_all_orientation',...
% %         'spatial', 'temporal', 'orientation');
end
%   keyboard
