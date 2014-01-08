function [] = PCA_select(classes_only, date_only_cells, All_TC, All_ACF, All_norm_RF, All_FR, All_norm_DOT, All_norm_RL, color, unique_dates)
% This function should behave somewhat like Vision's Show Classification
% Plots. It obtains arrays of classes/properties from properties_by_class.m and
% class_properties_graph.m. Using rectangle select, it will perform PCA on
% the unselected data points (in show_classification_plots.m).  Afterward, 
% it splits each selection into it's own class.  It then graphs things by class.
%
% ENZ, Fall 2012

%% split ONs from OFFs
ON_index = [];
OFF_index = [];

% % % % % % %       [ACF_coeff,ACF_score] = princomp(All_ACF);  % All ON/OFF cells
% % % % % % %       [TC_coeff,TC_score] = princomp(All_TC);

for i = 1: length(classes_only)
    if strcmp(classes_only(i,2), 'N')
        ON_index = [ON_index i];
    else
        OFF_index = [OFF_index i];
    end
end


%% Supervector PCA
keyboard
% supervector = cat(2, All_TC, All_ACF, All_norm_RF);
supervector = cat(2, All_TC, All_ACF, All_FR);



[coeff,score] = princomp(supervector);

figure
scatter3(score(:,1), score(:,2), score(:,3))
xlabel('PC1', 'FontSize', 20); ylabel('PC2', 'FontSize', 20); zlabel('PC3', 'FontSize', 20);
title('PCA on Concatinated Vector', 'FontSize', 25);


% figure
% scatter3(score(:,2), score(:,3), score(:,4))


[coeff,score] = princomp(supervector(ON_index, :));
figure
scatter3(score(:,1), score(:,2), score(:,3))
xlabel('PC1', 'FontSize', 20); ylabel('PC2', 'FontSize', 20); zlabel('PC3', 'FontSize', 20);
title('PCA on Cat Vector for ON classes', 'FontSize', 25);



[coeff,score] = princomp(supervector(OFF_index, :));
figure
scatter3(score(:,1), score(:,2), score(:,3))
xlabel('PC1', 'FontSize', 20); ylabel('PC2', 'FontSize', 20); zlabel('PC3', 'FontSize', 20);
title('PCA on Cat Vector for OFF classes', 'FontSize', 25);


%% Classification plot maker

% ON PLOTS
[~, ON_graph_index] = show_classification_plots(ON_index, All_TC, All_ACF, All_norm_RF); %"outside selection" array
 
keyboard
% OFF PLOTS
[~, OFF_graph_index] = show_classification_plots(OFF_index,All_TC, All_ACF, All_norm_RF);
 
 savefile = strcat('PCA_indices_', datestr(now,30), '.mat');
 save(savefile, 'ON_graph_index', 'OFF_graph_index');
 
 keyboard

%% GRAPHING NEW CLASSES
%% Time Courses
figure 
for i = 1 : length(OFF_graph_index(:,1))
    index_temp = logical(OFF_graph_index(i,:));
    index = OFF_index(index_temp);
    subplot(4,4,i);
    hold on
    
    for j = 1: length(index)
        
        color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
        
        plot(All_TC(index(j),:)'+ (((color_index/2)-1)/5), 'Color',  color(color_index,:), 'LineWidth', 1.25);
        
%         plot(zeros(1,25),':');
%         plot((1:250),zeros(1,250)+(((color_index/2)-1)/5),':');
%         plot((1:241),zeros(1,241)+(((color_index/2)-1)/5),':');
        plot((1:length(All_TC)),zeros(1,length(All_TC))+(((color_index/2)-1)/5),':');

    end
    title(strcat('nc',num2str(i)));
    ylim([-0.7 2.3]);
end
hold on

end_of_OFF = i;
for i = 1: length(ON_graph_index(:,1))
    index_temp = logical(ON_graph_index(i,:));
    index = ON_index(index_temp);
    subplot(4,4,i+end_of_OFF);
    hold on
    
    for j = 1: length(index)
        
        color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
        plot(All_TC(index(j),:)'+ (((color_index/2)-1)/5), 'Color',  color(color_index,:), 'LineWidth', 1.25);

%         plot(All_TC(ON_index(index),:)', 'b');

    plot((1:length(All_TC)),zeros(1,length(All_TC))+(((color_index/2)-1)/5),':');
%     plot((1:241),zeros(1,241)+(((color_index/2)-1)/5),':');
    end
    title(strcat('nc',num2str(i)));
    
    ylim([-0.7 2.3]);
end

%% Autocorrelation Function
figure % Autocorrelation Function
for i = 1 : length(OFF_graph_index(:,1))
    index_temp = logical(OFF_graph_index(i,:));
    index = OFF_index(index_temp);
    subplot(4,4,i);
    hold on
    for j = 1: length(index)
        color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
        plot(All_ACF(index(j),:)'+ (((color_index/2)-1)/10), 'Color', color(color_index,:), 'LineWidth', 1.25);
        plot((0:length(All_ACF)),zeros(1,length(All_ACF)+1)+(((color_index/2)-1)/10),':');
    end
    ylim([0 1.6]);
    title(strcat('nc',num2str(i)));
    
end

hold on
end_of_OFF = i;

for i = 1 : length(ON_graph_index(:,1))
    index_temp = logical(ON_graph_index(i,:));
    index = ON_index(index_temp);
    subplot(4,4,i + end_of_OFF);
    hold on
    for j = 1 : length(index)
        color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
        plot(All_ACF(index(j),:)'+ (((color_index/2)-1)/10), 'Color', color(color_index,:), 'LineWidth', 1.25);
        plot((0:length(All_ACF)),zeros(1,length(All_ACF)+1)+(((color_index/2)-1)/10),':');    
    end
    title(strcat('nc',num2str(i)));
    ylim([0 1.6]);
end

%% All Flash responses
if ~isempty(All_FR)
    figure % Autocorrelation Function
    for i = 1 : length(OFF_graph_index(:,1))
        index_temp = logical(OFF_graph_index(i,:));
        index = OFF_index(index_temp);
        subplot(4,4,i);
        hold on
        for j = 1: length(index)
            color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
            plot(All_FR(index(j),:)'+ (((color_index/2)-1)/5), 'Color', color(color_index,:), 'LineWidth', 1.25);
            plot((0:length(All_FR)),zeros(1,length(All_FR)+1)+(((color_index/2)-1)/5),':');
        end
        ylim([0 1.6]);
        title(strcat('nc',num2str(i)));
        
    end
    
    hold on
    end_of_OFF = i;
    
    for i = 1 : length(ON_graph_index(:,1))
        index_temp = logical(ON_graph_index(i,:));
        index = ON_index(index_temp);
        subplot(4,4,i + end_of_OFF);
        hold on
        for j = 1 : length(index)
            color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
            plot(All_FR(index(j),:)'+ (((color_index/2)-1)/5), 'Color', color(color_index,:), 'LineWidth', 1.25);
            plot((0:length(All_FR)),zeros(1,length(All_FR)+1)+(((color_index/2)-1)/5),':');
            
        end
        title(strcat('nc',num2str(i)));
        ylim([0 1.6]);
    end
end


%% Plotting all on one graph
figure % Time Courses
for i = 1 : length(OFF_graph_index(:,1))
    index = logical(OFF_graph_index(i,:));
    subplot(4,4,i);
    plot(All_TC(OFF_index(index),:)', 'b');
    title(strcat('nc',num2str(i)));
    
    ylim([-0.7 0.7]);
    hold on
    plot(zeros(1,250),':');
    
end
hold on
end_of_OFF = i;
for i = 1: length(ON_graph_index(:,1))
    index = logical(ON_graph_index(i,:));
    subplot(4,4,i+end_of_OFF);
    plot(All_TC(ON_index(index),:)', 'b');
    %     title(strcat('nc',num2str(i+end_of_OFF)));
    
    title(strcat('nc',num2str(i)));
    
    ylim([-0.7 0.7]);
    hold on
    plot(zeros(1,250),':');
    
end

figure % Autocorrelation Function
for i = 1 : length(OFF_graph_index(:,1))
    index = logical(OFF_graph_index(i,:));
    subplot(4,4,i);
    plot(All_ACF(OFF_index(index),:)', 'b');
    title(strcat('nc',num2str(i)));
    
end
hold on
end_of_OFF = i;

for i = 1 : length(ON_graph_index(:,1))
    index = logical(ON_graph_index(i,:));
    subplot(4,4,i + end_of_OFF);
    plot(All_ACF(ON_index(index),:)', 'b');
    %     title(strcat('nc',num2str(i + end_of_OFF)));
    
    title(strcat('nc',num2str(i)));
    
    
    hold on
    plot(zeros(1,25),':');
end

% FR
figure % Autocorrelation Function
for i = 1 : length(OFF_graph_index(:,1))
    index = logical(OFF_graph_index(i,:));
    subplot(4,4,i);
    plot(All_FR(OFF_index(index),:)', 'b');
    title(strcat('nc',num2str(i)));
    
end
hold on
end_of_OFF = i;

for i = 1 : length(ON_graph_index(:,1))
    index = logical(ON_graph_index(i,:));
    subplot(4,4,i + end_of_OFF);
    plot(All_FR(ON_index(index),:)', 'b');
    %     title(strcat('nc',num2str(i + end_of_OFF)));
    
    title(strcat('nc',num2str(i)));
    
    
    hold on
    plot(zeros(1,25),':');
end


keyboard
%% 3D Scatter plots

color = zeros(8,3);
n = 0;
for i = 0:1 %0.1:0.4:1
    for j = 0:1
        for k = 0:1
            n = n+1;
        color(n,:) = [i, j, k];  
        end
    end
end


for i = 1: length(ON_graph_index(:,1)) + length(OFF_graph_index(:,1))
    class_list{i} = strcat('nc',num2str(i));
end


figure
for i = 1:length(ON_graph_index(:,1))
    
    index = logical(ON_graph_index(i,:));
%     color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
    
    scatter3(All_norm_DOT(ON_index(index)), All_norm_RL(ON_index(index)), All_norm_RF(ON_index(index)),150, color(i,:), 'filled', 'o');
    text(All_norm_DOT(ON_index(index)), All_norm_RL(ON_index(index)), All_norm_RF(ON_index(index)), date_only_cells(ON_index(index)));
    hold on
end
xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20);
title('Plot for Normalized ON Average Properties','FontSize', 25);
legend(class_list(1:length(ON_graph_index(:,1))));

figure
for i = 1:length(OFF_graph_index(:,1))
    
    index = logical(OFF_graph_index(i,:));
%     color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
    
    scatter3(All_norm_DOT(OFF_index(index)), All_norm_RL(OFF_index(index)), All_norm_RF(OFF_index(index)),...
        150, color(i,:), 'filled', 'o');
    text(All_norm_DOT(OFF_index(index)), All_norm_RL(OFF_index(index)), All_norm_RF(OFF_index(index)), date_only_cells(OFF_index(index)));
    hold on
end

xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20);
title('Plot for Normalized OFF Average Properties','FontSize', 25);
legend(class_list(1:length(OFF_graph_index(:,1))));
% legend('nc1', 'nc2', 'nc3', 'nc4', 'nc5', 'nc6', 'nc7', 'nc8', 'nc9', 'nc10', 'nc11');

keyboard
end
