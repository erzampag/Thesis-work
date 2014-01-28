function [] = Directional_clusterer(date,datafile)


% file = strcat('/Users/erinzampaglione/Documents/Lab_Work/DSCells/DS_saveddata_no_OS/', date, '/', strcat(datafile, '.mat'));

% file = strcat('/Users/erinzampaglione/Documents/Lab_Work/DS_param/', date, '/', strcat(datafile, '.mat'));


file = strcat('/Users/erinzampaglione/Documents/Lab_Work/DSOSCells/', date, '/', strcat(datafile, '.mat'));

load(file)

keyboard
%% different for single paramter DS vs multi parameter DS

if iscell(AllCells)
    
    
    
    
    
else
    ON_OFF_DSindex = find(ratio(:) > 1 & DSI(:) > 0.5 & MaxSpikes(:) > 100 & noiseratio(:) > 5);
    DSindex = find(DSI(:) > 0.5 & MaxSpikes(:) > 100 & noiseratio(:) > 5);
end



%%
figure
p = polar(AllCells(DSindex,1), AllCells(DSindex,2), 'ob');
set(p,'LineWidth',1.5)
hold on
p = polar(AllCells(ON_OFF_DSindex,1), AllCells(ON_OFF_DSindex,2), 'ok');
set(p,'LineWidth',1.5)




pts = zeros(length(DSindex),2);
%         = zeros(length(ON_OFF_DSindex),1);
for  i = 1:length(DSindex)
    
    pts(i,1) = cos(AllCells(DSindex(i),1));
    pts(i,2) = sin(AllCells(DSindex(i),1));
    
end

keyboard
%% Instead of using the kmeans clusterer, I can just put these points in to
% the PCA clusterer and group them in that manner.  What is left to do then
% is make the output identical (ie "cluster assignments" rather than the
% rows of indices, and then use the same textfile generator.
    

    [~, tuning_graph_index] = show_classification_plots((1:length(DSindex)),...
        zeros(length(DSindex)), zeros(length(DSindex)),pts(:,1), pts(:,2));
    
    keyboard
    figure % Time Courses
    for i = 1 : length(tuning_graph_index(:,1))
        index = logical(tuning_graph_index(i,:));
        subplot(3,3,i);
        scatter(pts(index,1), pts(index,2), 'b');
        hold on

        title(strcat('nc',num2str(i)));
    end

    % re-write these as cluster assignments!!!!
    
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


%% This section uses the k-means clusterer, which automatically finds k clusters 
% based on how close each point is to the mean of the cluster

%     [cluster_assignments, error] = my_kmeans(pts, 4);

best_error = 1000;
for i  = 1: 100
    
    [cluster_assignments, error] = my_kmeans(pts, 4);
    
    if error < best_error
        
        cluster_assignments_KEEP = cluster_assignments;
        best_error = error
        
% % %         error_KEEPD = error;
    end
    
end
cluster_assignments = cluster_assignments_KEEP;

% Diagnostic from my_kmeans, just to look at final cluster assginments
    colors = [0 0 0; 0 0 1; 1 0 0; 0 1 0; 1 0 1; 0 1 1; 1 1 0]; %kbrgmcy
    figure; xlim([-1.2 1.2]); ylim([-1.2 1.2])
    hold on
    for i = 1: length(pts)
        plot(pts(i,1),pts(i,2), 'o', 'Color', colors(cluster_assignments(i), :))
    end




% [cluster_assignments, error] = my_kmeans(pts, 5)

keyboard


%% Make classification textfile
    keyboard
    formatted_sets = {};
    counter = 1;
    for i = 1 : length(idList)
        if counter <= length(DSindex) % go through all the DS cells
            if i == DSindex(counter) % if the neuron ID is one of the DS neurons
                for j = 1: length(cluster_assignments)
                    
                formatted_sets{i,1} =  idList(i);
                formatted_sets{i,2} = ['All/DS/cluster' num2str(cluster_assignments(counter))];
                
                end
                counter = counter+1;
            else % if the neuron ID is not a DS cell
                formatted_sets{i,1} =  idList(i);
                formatted_sets{i,2} = 'All/Other';
            end
        else % all cells after DS list is exhausted
            formatted_sets{i,1} =  idList(i);
            formatted_sets{i,2} = 'All/Other';
        end
    end
    
% % %     filename = ['/Users/erinzampaglione/Documents/processed_data/',date, '/',...
% % %         'data000-map-', datafile, '/DS_CLUSTERED_classification_matlab', '.txt'];
    
%         filename = ['/Users/erinzampaglione/Documents/processed_data/',date, '/',...
%          datafile, '/DS_CLUSTERED_classification_matlab', '.txt'];
     
     
         
    filename = ['/Users/erinzampaglione/Documents/processed_data/',date, '/',...
         datafile, '/DS_CLUSTERED_matlab', datestr(now,30), '.txt'];


    fileID = fopen(filename, 'w'); %% datestr(now,30)
    
    for j = 1 : length(formatted_sets)
        fprintf(fileID, '%d %s\n', formatted_sets{j,:});
    end
    
    fclose(fileID);

%%
figure
plot(pts(:,1),pts(:,2), 'o')
xlim([-1.2 1.2])
ylim([-1.2 1.2])



end