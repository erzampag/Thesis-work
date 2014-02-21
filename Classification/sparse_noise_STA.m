function [] = sparse_noise_STA(date, datafile)

date = '2014-01-29-0';
datafile = 'data003';



% import neuron info
[idList, classes, spike_times, EIx, EIy, ~, ~, neuronFile] = import_neuron_info(...
    '/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
    '/Users/erinzampaglione/Documents/processed_data/', date, datafile);

% import stimulus info

file_name = ['/Users/erinzampaglione/Documents/processed_data/' date  '/stimuli/' 'Sparse_noise.txt'];

stimFile = dlmread(file_name);




% import timing info (either use TTL pulses or VBL stamps if available)

temp2=neuronFile.getTTLTimes();% create matrix of TTL pulse times
TTL=double(temp2);%make TTL pulse matrix double

for i = 2:length(TTL)
TTL_diff(i-1) = TTL(i) - TTL(i-1);
end


% calculate STAs for each neuron
for NeuronID  = idList
        temp=neuronFile.getSpikeTimes(NeuronID);
    spikeTimes=double(temp);%converts the int matrix "spikeTimes" to a double matrix
    
end
