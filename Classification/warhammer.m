
% This script is the warhammer.
counter = 0;
 supervector = {};

for i = [1 2 4 5 10 15 17 19 21 23] % Loop over Retinas (Useable data: has OFFLBT, not nec. CRGs)
% i
    if i == 16 || i== 17 || i == 18 || i == 19
        [data, indices, celltypes, average_properties, norm_ave_properties] = properties_by_class(textdata{i,1}, 'data001');
    else
        [data, indices, celltypes, average_properties, norm_ave_properties] = properties_by_class(textdata{i,1}, 'data000');
    end
    counter = counter +1;
    
    
    for j = 1 : length(data)
        if ~strcmp(data{j,2},'All/Duplicates') && ~strcmp(data{j,2}(1:8), 'All/weak')
            
            supervector{counter}(j, :) = cat(2, i,  data{j, 7}, data{j, 6}, data{j,5});
            
        end
    end
  
%     keyboard
    supervector{counter}(all(supervector{counter}==0,2),:) = [];
    
  
end

all_neuron_vectors = [];

keyboard

for i = 1: length(supervector)
all_neuron_vectors = cat(1,all_neuron_vectors, supervector{i});

end


[coeff,score] = princomp(all_neuron_vectors(:, 2:end));

figure
scatter(score(:,1), score(:,2));

figure
scatter3(score(:,1), score(:,2), score(:,3));
xlabel('PC1', 'FontSize', 20); ylabel('PC2', 'FontSize', 20); zlabel('PC3', 'FontSize', 20);
title('PCA on Cat Vector for All Neurons', 'FontSize', 25);

All_TC = all_neuron_vectors(:, 2:26);
All_ACF = all_neuron_vectors(:, 27:226);
All_norm_RF = all_neuron_vectors(:, 227);

[nc] = show_classification_plots((1:length(all_neuron_vectors)), All_TC, All_ACF, All_norm_RF);

    graph_index = []; % Once clustering is done, use the "outside selection" array to make a logical array corresponding to each selection
for i = 1: length(nc)
    for j = 2:length(nc(:,i))
        if nc(j,i) ~= nc(j-1,i)
            graph_index(j-1,i) = true;
        else
            graph_index(j-1,i) = false;
        end
    end
end


ON_index = [];
OFF_index = [];
for i = 1: length(graph_index)
    if logical(graph_index(1,i))
        ON_index = [ON_index i];
    else
        OFF_index = [OFF_index i];
    end
end

keyboard






[nc] = show_classification_plots(ON_index, All_TC, All_ACF, All_norm_RF);

[nc] = show_classification_plots(OFF_index, All_TC, All_ACF, All_norm_RF);


keyboard


%%
% [coeff,score] = princomp(all_neuron_vectors(ON_index, 2:end));
[coeff,score] = princomp(all_neuron_vectors(ON_index, 2:end-1));


figure
scatter(score(:,1), score(:,2));

figure
scatter3(score(:,1), score(:,2), score(:,3));

% [coeff,score] = princomp(all_neuron_vectors(OFF_index, 2:end));
[coeff,score] = princomp(all_neuron_vectors(OFF_index, 2:end-1));



