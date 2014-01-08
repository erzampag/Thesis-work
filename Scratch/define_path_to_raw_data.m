function [pfile, neuronFile] = define_path_to_raw_data(date, datafile)

% % javaaddpath('/Users/erinzampaglione/Desktop/workspace/vision8/Vision.jar');

javaaddpath('/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar');



% %define path to raw data file
% full_path=strcat('/Volumes/TEMP_MOUSE/', date, '/', datafile, '/' , datafile, '.bin');
% rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

% spikefile=strcat('/Volumes/TEMP_MOUSE/', date, '/', datafile, '/' , datafile, '.spikes');
% spikeFile=edu.ucsc.neurobiology.vision.io.SpikeFile(spikefile);

% % % % % % % neuronfile=strcat('/Volumes/TEMP_MOUSE/', date, '/', datafile, '/' , datafile, '.neurons');
% % % % % % % neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(neuronfile);
% % % % % % %
% % % % % % % param_path = strcat('/Volumes/TEMP_MOUSE/', date, '/', datafile, '/' , datafile,'.params');
% % % % % % % pfile = edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);

%%
neuronfile=strcat('/Users/erinzampaglione/Documents/processed_data/', date, '/', datafile, '/' , datafile, '.neurons');
neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(neuronfile);

param_path = strcat('/Users/erinzampaglione/Documents/processed_data/', date, '/', datafile, '/' , datafile,'.params');
pfile = edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);


%%
% % % % neuronfile=strcat('/Users/erinzampaglione/Documents/processed_data/datafromextHD/', date, '/', datafile, '/' , datafile, '.neurons');
% % % % neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(neuronfile);
% % % % 
% % % % param_path = strcat('/Users/erinzampaglione/Documents/processed_data/datafromextHD/', date, '/', datafile, '/' , datafile,'.params');
% % % % pfile = edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);






end