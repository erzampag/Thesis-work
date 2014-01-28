function [file, DS_indices, index, DS_index, mean_rates, max_rates, high_index...
    color, x_vector, y_vector, mean_rates_spatial, mean_rates_spatial_norm,mean_rates_temporal...
   mean_rates_temporal_norm, max_slice_spatial, max_slice_spatial_norm, max_slice_temporal...
   max_slice_temporal_norm, concat_no_orient, concat_no_orient_norm, color_long, mean_rate_for_each_dir...
   concat_orient, concat_orient_norm] = calculate_DS_param_props(date, data_file)

% separating out the property calculations so that it's not in the same place as
% graphing - i think this is better than keeping the parameters calculated
% by the graph they are used it if i end up wanting to use the same
% paramter in a different graph 

file = (['/Users/erinzampaglione/Documents/Lab_work/DS_param/' date '/' data_file '.mat']);
load(file)


%% DS_indices are made as just a first cut to the hand-picking of DS neurons

max_rate_absolute = cell(length(idList),1);
% build max spike rate array
for q  = 1: length(idList)
    for i = 1: length(spatial)
        for j = 1: length(temporal)
       max_rate_absolute{q}(i,j) = max(spike_rate_all_orientation{q}(i,j,:));     
%        mean_rate_absolute     
        end
    end
end

DS_indices = cell(length(idList),1); % for each neuron, we have a list of which spatial and temporal indices were automatically 
                                       % determined to be DS

for i  = 1 : length(idList)
    %    max_DSI(i) = max(max(DSI{i}));
    [spat, tem] = find(DSI{i} > 0.5 & max_rate_absolute{i} > 2); % Max rate for this SP/TEMP, any orientation is 2 sp/sec
%     [spat, tem] = find(DSI{i} > 0.5 &  max(max(max(spike_rate_all_orientation{i}))) > 2); % maximum spike rate ANYWHERE was greater than 2sp/sec
%     [spat, tem] = find(DSI{i} > 0.45 & MaxSpikes{i} > 100); % right now this is at 100, because even though there were
    % 25 sec of display, i mutiply by 50 so the threshold is also mutiplied by 2
    DS_indices{i} = [spat, tem];
    
end
DS_neuronIDs = zeros(length(idList), 1);
for i = 1 : length(idList)
%     DS_neuronIDs(i) = ~isempty(DS_indices{i}); % To make the first DS cut, there must be at least one SP/TEMP combination with DSI > 1 and 
                                                % Max rate ANYWHERE was greater than 2 spikes/sec
    
    DS_neuronIDs(i) = ~(length(DS_indices{i}(:,1)) < 2); % need at least 2 SP/TEMPS that fulfill requirement
end

index = (1: length(idList))';

DS_index = index(logical(DS_neuronIDs));

%% mean rates is the average response over all directions for a given spat/temp

mean_rates = cell(length(idList), 1); % 
max_rates = cell(length(idList), 1);
for q = 1:length(idList)
    for i = 1: length(spatial)
        for j = 1 : length(temporal)
            mean_rates{q}(i,j) = mean(spike_rate_all_orientation{q}(i,j,:));
            max_rates{q}(i,j) = max(spike_rate_all_orientation{q}(i,j,:));
        end
    end
end
 color = [0 0 0; 0 0 1; 1 0 0; 0 1 0; 1 0 1; 0 1 1; 1 1 0]; %kbrgmcy


%% remake X and Y's for compass plot
x_vector = cell(length(idList),1);
y_vector = cell(length(idList),1);

for q = 1:length(idList)
    for i = 1: length(spatial)
        for j = 1:length(temporal)
            x_vector{q}(i,j) = (AllCells{q,2}(i,j).*cos(AllCells{q,1}(i,j)));
            y_vector{q}(i,j) = (AllCells{q,2}(i,j).*sin(AllCells{q,1}(i,j)));
            
        end
    end
end


%% OKAY, here we want to condense rate(spatial) and rate(temporal) to a single line per neuron

mean_rates_spatial = zeros(length(idList), length(spatial));
mean_rates_spatial_norm = zeros(length(idList), length(spatial));
mean_rates_temporal = zeros(length(idList), length(temporal));
mean_rates_temporal_norm = zeros(length(idList), length(temporal));

max_slice_spatial = zeros(length(idList), length(spatial));
max_slice_spatial_norm = zeros(length(idList), length(spatial));
max_slice_temporal = zeros(length(idList), length(temporal));
max_slice_temporal_norm = zeros(length(idList), length(temporal));

for i = 1 : length(idList)
    mean_rates_spatial(i,:) = mean(mean_rates{i},2);
    mean_rates_temporal(i,:) = mean(mean_rates{i},1);
    
    mean_rates_spatial_norm(i,:) = mean_rates_spatial(i,:)/sum(mean_rates_spatial(i,:));
    mean_rates_temporal_norm(i,:) = mean_rates_temporal(i,:)/sum(mean_rates_temporal(i,:));
    
    high_point = max(max(mean_rates{i}));
    [spat_ind, temp_ind] = find(mean_rates{i} == high_point(1));
    
    max_slice_spatial(i,:) = mean_rates{i}(:, temp_ind(1));
    max_slice_temporal(i,:) = mean_rates{i}(spat_ind(1), :);
    
end

% additionally, what i'd like out of this is WHICH temp/spat period
% corresponds to the max and make ERROR bars to indicate if that cell is
% reliable

%% find the mean rate for just the few sp/temps that surround the highest point


for q = 1: length(idList)

    high_point = max(max(mean_rates{q}));
    [spat_ind, temp_ind] = find(mean_rates{q} == high_point(1),1);
    
    high_index(q,1) = spat_ind; 
    high_index(q,2)=  temp_ind;
%     
end

% % % % (high_index(q,1)-1)+(high_index(q,2)-1)*5 +1 % CONVERT TO VERT INDEX

%% make the concatinated vector like in Matt's thesis

concat_no_orient = zeros(length(idList), length(spatial)*length(temporal));
concat_no_orient_norm = zeros(length(idList), length(spatial)*length(temporal));


for i = 1 : length(idList)
  concat_no_orient(i, :) = cat(2, mean_rates{i}(1,:), mean_rates{i}(2,:),...
      mean_rates{i}(3,:), mean_rates{i}(4,:), mean_rates{i}(5,:));
  
  concat_no_orient_norm(i,:) = concat_no_orient(i,:)/sum(concat_no_orient(i,:));
        
end

%% Make the concatinated vector like in Matt's thesis but with the Directional information included

concat_orient = zeros(length(idList), length(spatial)*length(temporal)*length(orientation));
concat_orient_norm = zeros(length(idList), length(spatial)*length(temporal)*length(orientation));

for q = 1:length(idList)
    temp_vector = [];
    for i = 1:length(spatial)
        for j = 1:length(temporal)
            temp_vector = [temp_vector squeeze(spike_rate_all_orientation{q}(i,j,:))'];
        end 
        
    end
    concat_orient(q,:) = temp_vector;
    concat_orient_norm(q,:) = temp_vector/sum(temp_vector);
end
    



%% average over all spat/temps, to just have direction
% keyboard
mean_rate_for_each_dir = zeros(length(idList), length(spike_rate_all_orientation{q}(1,1,:)));
for q = 1: length(idList)
    
    for k = 1:length(spike_rate_all_orientation{q}(1,1,:))
       mean_rate_for_each_dir(q,k) =  mean(mean(spike_rate_all_orientation{q}(:,:,k),1),2);
    end
    
end


%% color vector used for plotting all the concatinated vectors on one graph so that we can
% distinguish individual 

color_long = zeros(27,3);
n = 0;
for i = [0,0.5,1]
    for j = [0,0.5,0.75]
        for k = [0,0.5,1]
            n = n+1;
        color_long(n,:) = [i, j, k];  
        end
    end
end

figure
scatter(1:length(color_long), 1:length(color_long),100, color_long, 'filled')
color_long = cat(1,color_long, color_long);

end