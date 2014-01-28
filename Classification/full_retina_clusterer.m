function [] = full_retina_clusterer(id_List, PCA_var1, PCA_var2, scatter_var1, scatter_var2)

% not sure if i'll use this - idea is to have something that will take all
% neurons from a retina and cluster them by PCA, or just a scatter plot


    for i = 1: length(tuning_graph_index(1,:))
        for j = 1:length(tuning_graph_index(:,1))
            if tuning_graph_index(j,i) == 0
                continue
            else
                cluster_assignments(i) = j;
            end
        end
    end

keyboard


  keyboard
    formatted_sets = {};
    counter = 1;
    for i = 1 : length(idList)
%         if counter <= length(DSindex) % go through all the DS cells
%             if i == DSindex(counter) % if the neuron ID is one of the DS neurons
%                 for j = 1: length(cluster_assignments)
                    
                formatted_sets{i,1} =  idList(i);
                formatted_sets{i,2} = ['All/cluster' num2str(cluster_assignments(counter))];
                
%                 end
                counter = counter+1;
%             else % if the neuron ID is not a DS cell
%                 formatted_sets{i,1} =  idList(i);
%                 formatted_sets{i,2} = 'All/Other';
%             end
%         else % all cells after DS list is exhausted
%             formatted_sets{i,1} =  idList(i);
%             formatted_sets{i,2} = 'All/Other';
%         end
    end
    
        filename = ['/Users/erinzampaglione/Documents/processed_data/',date, '/',...
         data_file, '/CLUSTERED_matlab', datestr(now,30), '.txt'];


    fileID = fopen(filename, 'w'); %% datestr(now,30)
    
    for j = 1 : length(formatted_sets)
        fprintf(fileID, '%d %s\n', formatted_sets{j,:});
    end
    
    fclose(fileID);
    
end