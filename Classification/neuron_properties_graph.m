


[data, indices, celltypes, average_properties, norm_ave_properties] = properties_by_class(textdata{i,1}, 'data000');


for i = 1 : length(celltypes)
    
    
    average_properties{i,1} = celltypes{i};
    
    if isempty(indices{i})
        continue
    end
    
    degree_transiency = [];
    response_latency = [];
    receptive_field = [];
    %     autocorrelation = {};
    %     time_course = {};
    autocorrelation = [];
    time_course = [];
    
    for j = 1 : length(indices{i}) % temporary list of properties of each neuron of that celltype
        degree_transiency = [degree_transiency data{[indices{i}(j)],3}]; % DOT
        response_latency = [response_latency data{[indices{i}(j)],4}]; % RL
        receptive_field = [receptive_field data{[indices{i}(j)],5}]; % RF
        %         autocorrelation{j} = data{[indices{i}(j)],6};
        %         time_course{j} = data{[indices{i}(j)],7};
        autocorrelation = [autocorrelation; data{[indices{i}(j)],6}];
        time_course = [time_course; data{[indices{i}(j)],7}];
        
    end
    
    
    
end

