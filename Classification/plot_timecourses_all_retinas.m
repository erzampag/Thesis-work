function [] = plot_timecourses_all_retinas(x_axes, y_axes, unique_classes, unique_dates, classes_only_cells, date_only_cells, color)



%% plotting all timecourses on one graph
figure
hold on
for i = 1:length(y_axes(:,1)) 
    plot(x_axes(i,:), y_axes(i,:))
end
    ylim([-1.1 1.1]);


%% plotting timecourses by class
figure
for i = 1 : length(unique_classes)
    
    index = find(strcmp(classes_only_cells(:), unique_classes(i)));
%     subplot(4,4,i); % suplot(5,4,i) if there's more than 16 cell types
    subplot(2,3,i); % suplot(5,4,i) if there's more than 16 cell types

    hold on
    
    for j = 1 : length(index)
                color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));

        plot(x_axes(index(j), :), y_axes(index(j), :), 'Color', color(color_index,:)); % time course from normalized fit
        hold on
        plot((-400:-1),zeros(1,400),':');
        
    end
    title(unique_classes(i));
    ylim([-1.1 1.1]);
    xlim([-400 0]);
    
end


%% plotting timecourses by classes and retinas

figure
for i = 1 : length(unique_classes)
    index = find(strcmp(classes_only_cells(:), unique_classes(i)));
    
%     subplot(4,4,i); % suplot(5,4,i) if there's more than 16 cell types
    subplot(2,3,i); % suplot(5,4,i) if there's more than 16 cell types
    hold on
    for j = 1 : length(index)
        color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
        
%         plot(x_axis_raw, All_TC(index(j),:)'+ (((color_index/2)-1)/5), 'Color', color(color_index,:), 'LineWidth', 1.25); %color by date

% %         plot(linspace(-400,0,length(All_TC_fit_norm(index(j),:))), All_TC_fit_norm(index(j),:)'+ (((color_index/2)-1)/5),...
% %             'Color', color(color_index,:), 'LineWidth', 1.25); % SCALED color by date FIT
        
                plot(x_axes(index(j),:), y_axes(index(j),:)'+ (((color_index/2)-1)/5),...
            'Color', color(color_index,:), 'LineWidth', 1.25); % SCALED X AND Y color by date FIT
        
        
%         plot(linspace(-330,0,length(All_TC_spline_norm(index(j),:))), All_TC_spline_norm(index(j),:)'+ (((color_index/2)-1)/5),...
%             'Color', color(color_index,:), 'LineWidth', 1.25); % SCALED color by date SPLINE

        plot((-400:-1),zeros(1,400)+(((color_index/2)-1)/5),':');
        
%         text(-400,(((color_index/2)-1)/5), unique_dates{color_index/2});All_DOT        
    end
    title(unique_classes(i), 'FontSize', 20);
    ylim([-0.7 (length(unique_dates)-1)/5+0.7]);
    
    xlim([min(min(x_axes)) 0]);

end


end
