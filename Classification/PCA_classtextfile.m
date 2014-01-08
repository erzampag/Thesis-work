function [] = PCA_classtextfile(date, data_file, PCA_indices, idList)

  
  cluster_num = zeros(length(idList),1);
  
    
    for i = 1: length(idList)
        for j = 1: length(PCA_indices(:,1))
            if PCA_indices(j,i) == 1
                cluster_num(i,1) = j;
            else
                continue
            end
        end
    end
  
    
formatted_sets = {};
    for j = 1 : length(idList)
                formatted_sets{j,1} =  idList(j);
                formatted_sets{j,2} = ['All/cluster_' num2str(cluster_num(j))];
    end
    


% % % % % % % formatted_sets = {};
% % % % % % %     counter = 1;
% % % % % % %     for j = 1 : length(idList)
% % % % % % %         if counter <= length(DScell_index) % go through all the DS cells
% % % % % % %             if j == DScell_index(counter) % if the neuron ID is one of the DS neurons
% % % % % % %                 formatted_sets{j,1} =  idList(j);
% % % % % % %                 formatted_sets{j,2} = ['All/DS' num2str(cluster_num)];
% % % % % % %                 counter = counter+1;
% % % % % % %             else % if the neuron ID is not a DS cell
% % % % % % %                 formatted_sets{j,1} =  idList(j);
% % % % % % %                 formatted_sets{j,2} = 'All/Other';
% % % % % % %             end
% % % % % % %         else % all cells after DS list is exhausted
% % % % % % %             formatted_sets{j,1} =  idList(j);
% % % % % % %             formatted_sets{j,2} = 'All/Other';
% % % % % % %         end
% % % % % % %     end
% % % % % % %     


     filename = ['/Users/erinzampaglione/Documents/processed_data/',date, '/',...
         data_file, '/PCA_class_matlab_', datestr(now,30) '.txt'];
    
    fileID = fopen(filename, 'w'); %% datestr(now,30)
    
    for j = 1 : length(formatted_sets)
        fprintf(fileID, '%d %s\n', formatted_sets{j,:});
    end
    
    fclose(fileID);
    
    
end