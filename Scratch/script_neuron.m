%%

date = '2013-10-17-0';
data_file = 'data001-map-data003';


javaaddpath('/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar');

neuronfile=strcat('/Users/erinzampaglione/Documents/processed_data/', date, '/', data_file, '/' , data_file, '.neurons');
neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(neuronfile);

for i = 1: length(idList)
    
    temp=neuronFile.getSpikeTimes(idList(i));
    spiketimes=double(temp);%converts the int matrix "spikeTimes" to a double matrix
    
    Neuron_objects(i) = neuron(idList(i), spiketimes, spike_rate_all_orientation{i});
    keyboard
end
