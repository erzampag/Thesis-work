function [neuronIDs, classes, spiketimes, EIx, EIy, data, pfile, neuronFile] = import_neuron_info(java_path, processed_data_path, date, datafile)

% This function pulls information from Vision into Matlab and outputs the
% neuron IDs, classes, spike times, and EI positions x and y.  If you want
% to pull other info, it can be placed in the"data" array.
% All inputs are taken as strings.
%
% ex: import_neuron_info('/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
%    '/Users/erinzampaglione/Documents/processed_data/', '2011-02-23-0', 'data000');
%
% ENZ, Spring 2013
keyboard
%% Define paths
javaaddpath(java_path);

neuronfile=strcat(processed_data_path, date, '/', datafile, '/' , datafile, '.neurons');
neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(neuronfile);

param_path = strcat(processed_data_path, date, '/', datafile, '/' , datafile,'.params');
pfile = edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);

%% Import neuron IDs
idListbefore = [];
idListbefore=neuronFile.getIDList();

%Remake the vector of neuronIDs so it doesn't have weird negative numbers
idListPreD = [];

for i = 1:length(idListbefore)
    if idListbefore(i) <=0,
        continue
    else
        idListPreD = [idListPreD, idListbefore(i)];
    end
end
neuronIDs = idListPreD';

%% Import classes

classInfo = pfile.getClassIDs(); % returns a JavaHashmap Object

classes = cell(length(neuronIDs), 2); % each row is [ID ClassInfo]
for i = 1:length(neuronIDs)
    classes{i,1} = neuronIDs(i); % because it's a cell it needs a curly bracket
    classes{i,2} = classInfo.get(int32(neuronIDs(i)));
end
%% Import spike times

for i = 1:length(neuronIDs)
    temp=neuronFile.getSpikeTimes(neuronIDs(i));
    spiketimes{i,1}=double(temp);%converts the int matrix "spikeTimes" to a double matrix
    
end

%% Import neuron locations

% Get other data

data = zeros(length(idListPreD), 5);
%
% for i = 1:length(idListPreD)
%     data(i,1) = idListPreD(i);
%     data(i,2) = pfile.getDoubleCell(idListPreD(i), 'xOffDS');
%     data(i,3) = pfile.getDoubleCell(idListPreD(i), 'yOffDS');
%     data(i,4) = pfile.getDoubleCell(idListPreD(i), 'magOffDS');
%     data(i,5) = pfile.getDoubleCell(idListPreD(i), 'angOffDS');
% end


EIx = zeros(length(neuronIDs),1);
EIy = zeros(length(neuronIDs),1);
% 
% keyboard
% for i = 1:length(neuronIDs)
%     if isempty(pfile.getDoubleCell(neuronIDs(i), 'EIx0'))
%         continue
%     else
%         
%         EIx(i,1) = pfile.getDoubleCell(neuronIDs(i), 'EIx0');
%         EIy(i,1) = pfile.getDoubleCell(neuronIDs(i), 'EIy0');
%     end
% end




end