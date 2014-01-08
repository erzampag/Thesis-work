 

% 
% [neuronIDs, classes, spiketimes, EIx, EIy, data, pfile, neuronFile] =...
%     import_neuron_info('/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
%      '/Users/erinzampaglione/Documents/processed_data/', '2012-04-30-0', 'data002');
 
 
%  
% [neuronIDs, classes, spiketimes, EIx, EIy, data, pfile, neuronFile] =...
%     import_neuron_info('/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
%      '/Users/erinzampaglione/Desktop/DS_parameterscan/', '2012-04-27-0', 'data002');
 
 
 [neuronIDs, classes, spiketimes, EIx, EIy, data, pfile, neuronFile] =...
    import_neuron_info('/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
     '/Users/erinzampaglione/Documents/processed_data/', '2012-05-07-0', 'data002');
 
 
 
 gratingResponse = {};
for i = 1:length(neuronIDs)
%     data(i,1) = idListPreD(i);
%     data(i,2) = pfile.getDoubleCell(idListPreD(i), 'xOffDS');
%     data(i,3) = pfile.getDoubleCell(idListPreD(i), 'yOffDS');
%     data(i,4) = pfile.getDoubleCell(idListPreD(i), 'magOffDS');
%     data(i,5) = pfile.getDoubleCell(idListPreD(i), 'angOffDS');
%     
   gratingResponse{i} = pfile.getArrayCell(neuronIDs(i), 'gratingResponse');
end

keyboard