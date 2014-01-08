function [cluster_assignments, error] = my_kmeans(pts, num_means)

% initialize means
% keyboard
means = zeros(num_means, length(pts(1,:)));
for i  = 1: num_means
    means(i,:) = pts(round(1 + (length(pts)-1).*rand(1,1)), :);
end

dif_means = 1;

% keyboard
while dif_means > .00001
    
    % find all distances
    cluster_assignments = zeros(length(pts),1);
    for i = 1: length(pts)
        test_distances = zeros(length(means), 1);
        for j = 1: length(means)
            test_distances(j) = pdist([means(j,:); pts(i,:)]); % for each point, the distance from all means
        end
        
%         if length(unique(test_distances)) < length(test_distances)
%             fprintf('more than one test_distance\n')
% %             find(test_distances == min(test_distances))
%         end
        
        min_index = find(test_distances == min(test_distances));
        if length(min_index) > 1
            min_index = min_index(1); % if there's more than one minimum distance, choose the first one
        end
        cluster_assignments(i) = min_index;
    end
%     keyboard

    % take means
    for i = 1:length(means)
        new_means(i,:) = mean(pts(cluster_assignments(:) == i,:), 1);
    end

    
    
    
    % cluster plotter diagonstic
% % % % %     colors = [0 0 0; 0 0 1; 1 0 0; 0 1 0; 1 0 1; 0 1 1; 1 1 0]; %kbrgmcy
% % % % %     figure; xlim([-1.2 1.2]); ylim([-1.2 1.2])
% % % % %     hold on
% % % % %     for i = 1: length(pts)
% % % % %         plot(pts(i,1),pts(i,2), 'o', 'Color', colors(cluster_assignments(i), :))
% % % % %     end
% % % % %     % new means diagnostic
% % % % %     plot(means(:,1), means(:,2), 'xr', 'MarkerSize', 20)
% % % % %     plot(new_means(:,1), new_means(:,2), 'xk', 'MarkerSize', 20)
    
    
    
    
    % update step
    dif_means = new_means - means;
    
    dif_means = sum(sum(abs(dif_means)));
    means = new_means;
end

if unique(cluster_assignments) < num_means
%     cluster_assignments
    error = 1000;
%     fprintf('empty cluster - \n')
    return % does this work?
end

% error 
error = 0;
for i = 1: length(pts)
    error = error + (pdist([means(cluster_assignments(i),:); pts(i,:)]))^2;
end
% error
% keyboard
end