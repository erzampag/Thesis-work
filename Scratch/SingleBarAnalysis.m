function [idList, AllCells, DSI, DSI_error, ratio, noiseratio,...
    DScellsindex, noise, zerofreq, fwhh, reduced_chi_square, MaxSpikes, EIx, EIy]=...
    SingleBarAnalysis(date, datafile, stimfile2)


% This function imports information from Vision and analyses neuronal
% responses to single moving bars.

% Written by Erin Zampaglione, Summer 2012, adapted from drifting_squarewave_shorter_stimulus2.



%% Define path to raw data
% Edit this function to access Vision.jar and raw data in a different location

[pfile, neuronFile] = define_path_to_raw_data(date, datafile);
%%  *****PRINTFIGURES:
% 0 - print nothing but SAVE VARIABLES!,
% 1 - print everything,
% 2 - print PSTH, FFT, tuningcurve, vector for neuron,
% 3 - print tuning curve and vector together,
% 4 - just final whole retina analysis
% 5 - 
printfigures = 0;

%% Define variables for use in analysis


AllCells=[];





%% Get neuron IDs from Vision (before duplication removal)
%
idListPreD =  import_neuronIDs(neuronFile);

%% Get information from parameters file (classes, xOffDS, etc)
% 
[classes, data, EIx, EIy] = import_parameters(pfile, idListPreD);


keyboard


%% Remove Duplicates in IDlist and any other parameter
%
[idList, ind] = remove_duplicates(idListPreD, classes);

EIx = EIx(ind);
EIy = EIy(ind);
data = data(ind);
% % keyboard

%% Use TTL pulses to determine the start and end times of the stimuli (will need to be different depending on what we decide on stim)
%
% how often is there ttl pulses? give iit similar timing to multibar?
% % % % % % [start_times, End_times] = calculate_start_end(neuronFile);

keyboard

%% Start of the loop over neurons

%for given neuron, import the spikes times, just doing one right now so
%there's not too many plots popping up.

%TROUBLE = ...

fprintf('Analyzing Retina from %s \n', date)
for q=1:length(idList);
    % % for q = 1: length(TROUBLE);
    
    if mod(q,50)==0,
        fprintf('processing neuron # %d \n',q);
    end
    
    NeuronID=idList(q);  %% i just chose a DS looking cell!
    
    % %     NeuronID = idList(TROUBLE(q)); % if TROUBLE is the indices
    
    %    NeuronID = TROUBLE(q); % if TROUBLE is the neuron IDs
    
    %% 
    temp=neuronFile.getSpikeTimes(NeuronID);
    spikeTimes=double(temp);%converts the int matrix "spikeTimes" to a double matrix
    
    Direction=0:22.5:337.5; %?
    orientation=0:22.5:337.5; %?
    
    
    
    
    
    
end


%% After all neurons have been analyzed


end