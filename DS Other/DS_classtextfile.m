function [] = DS_classtextfile(idList, DScell_index, filename)    


formatted_sets = {};
    counter = 1;
    for j = 1 : length(idList)
        if counter <= length(DScell_index) % go through all the DS cells
            if j == DScell_index(counter) % if the neuron ID is one of the DS neurons
                formatted_sets{j,1} =  idList(j);
                formatted_sets{j,2} = 'All/DS';
                counter = counter+1;
            else % if the neuron ID is not a DS cell
                formatted_sets{j,1} =  idList(j);
                formatted_sets{j,2} = 'All/Other';
            end
        else % all cells after DS list is exhausted
            formatted_sets{j,1} =  idList(j);
            formatted_sets{j,2} = 'All/Other';
        end
    end
    
%     filename = ['/Users/erinzampaglione/Documents/processed_data/',textdata{i,1}, '/',...
%         'data000-map-', textdata{i,2}, '/DS_classification_matlab', '.txt'];
    
    fileID = fopen(filename, 'w'); %% datestr(now,30)
    
    for j = 1 : length(formatted_sets)
        fprintf(fileID, '%d %s\n', formatted_sets{j,:});
    end
    
    fclose(fileID);