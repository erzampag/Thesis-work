function [] = scatter3_by_classes(on_handle, off_handle, on_title, off_title, x, x_label, y, y_label, z, z_label,...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color)

% this function uses the matlab scatter function to color code and label
% classes by type and retina.

figure(on_handle);
figure(off_handle);
for j = 1:length(unique_classes)
    
    index = find(strcmp(classes_only_cells(:), unique_classes(j)));
    color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
    
    if strcmp(unique_classes{j}(2), 'N')
        figure(on_handle)
        scatter3(x(index), y(index), z(index),150, color(color_index,:), 'filled', 'o');
        
%         text(x(index), y(index), z(index), date_only_cells(index));
        
        set(gca,'FontSize', 15);
        hold on
    else
        figure(off_handle)
        scatter3(x(index), y(index), z(index),150, color(color_index,:), 'filled', 'o');
        
%         text(x(index), y(index), z(index), date_only_cells(index));
        
        set(gca,'FontSize', 15);
        hold on
    end
end

OFF_index = [];
for i = 1 : length(unique_classes)
    OFF_index = [OFF_index strcmp(unique_classes{i}(2), 'F')];
end

figure(on_handle)

legend(unique_classes(logical(~OFF_index)));
xlabel(x_label, 'FontSize', 20); ylabel(y_label, 'FontSize', 20); zlabel(z_label, 'FontSize', 20)
title(on_title,'FontSize', 25);

figure(off_handle)

legend(unique_classes(logical(OFF_index)));
xlabel(x_label, 'FontSize', 20); ylabel(y_label, 'FontSize', 20); zlabel(z_label, 'FontSize', 20)
title(off_title,'FontSize', 25);


end