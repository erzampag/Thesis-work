function All_TC_xynorm = TC_xy_scaling(All_TC_fit, x_axis, All_RL, date_only_cells, unique_classes, unique_dates,...
    classes_only_cells, full_description_cells, color)

keyboard

plot_timecourses_all_retinas(x_axis, All_TC_fit, unique_classes, unique_dates, classes_only_cells, date_only_cells, color)


%%%%%%% Take TC_fit and scale all to the same zero crossing

x_axis_norm_one = linspace(-400,0,length(All_TC_fit));
x_axis_norm_one = x_axis_norm_one([ones(length(classes_only_cells),1)],:); % repeat that vector for all time courses (y values)

master_retina = find(strcmp(full_description_cells(:), '2011-02-17-0-OFFLBT'));
for i = 1 : length(unique_dates)
    for j = 1 : length(date_only_cells)
        if strcmp(date_only_cells{j}, unique_dates{i})
            
            % Scaling to RL
            x_axis_norm_one(j,:) = x_axis(j,:).*(All_RL(master_retina) / All_RL(j));
            
        end
    end
end


plot_timecourses_all_retinas(x_axis_norm_one, All_TC_fit, unique_classes, unique_dates, classes_only_cells, date_only_cells, color)



%%%%%%%% Re-fit the time courses
All_TC_fit_norm_one = zeros(length(date_only_cells), length(All_TC_fit));
for i = 1 : length(date_only_cells)

[All_TC_fit_norm_one(i,:),~, ~, ~] = ...
    fit_time_course((x_axis_norm_one(i,:)'+400*25/24)*24/400, All_TC_fit(i,:)');

end

x_axis = linspace(-400,0,length(All_TC_fit_norm_one));
x_axis = x_axis([ones(length(classes_only_cells),1)],:); % repeat that vector for all time courses (y values)

plot_timecourses_all_retinas(x_axis, All_TC_fit_norm_one, unique_classes, unique_dates, classes_only_cells, date_only_cells, color)


%%%%%%% Normalize OFFLBT Y 
All_TC_xynorm = zeros(length(date_only_cells), length(All_TC_fit));

scaling_factor = zeros(length(date_only_cells), length(All_TC_fit));

for i = 1 : length(unique_dates)
    for j = 1 : length(date_only_cells)
        if strcmp(date_only_cells{j}, unique_dates{i})
            that_retina = strcmp(full_description_cells(:), strcat(unique_dates{i}, '-OFFLBT'));         
            
            scaling_factor(j,:) = (All_TC_fit_norm_one(master_retina,:) ./ All_TC_fit_norm_one(that_retina,:));
            
            All_TC_xynorm(j,:) = All_TC_fit_norm_one(j,:) .* (All_TC_fit_norm_one(master_retina,:) ./ All_TC_fit_norm_one(that_retina,:));
        end
    end
end


plot_timecourses_all_retinas(x_axis, All_TC_xynorm, unique_classes, unique_dates, classes_only_cells, date_only_cells, color)


%%%%%% Undo the X Scaling 
x_axis_undone = linspace(-400,0,length(All_TC_fit_norm_one));
x_axis_undone = x_axis_undone([ones(length(classes_only_cells),1)],:); % repeat that vector for all time courses (y values)

for i = 1 : length(unique_dates)
    for j = 1 : length(date_only_cells)
        if strcmp(date_only_cells{j}, unique_dates{i})
            
% % % % % %             % Scaling to RL
            x_axis_undone(j,:) = x_axis(j,:).*(All_RL(j) / All_RL(master_retina));      
            
        end
    end
end


plot_timecourses_all_retinas(x_axis_undone, All_TC_xynorm, unique_classes, unique_dates, classes_only_cells, date_only_cells, color)


%%%%%% RE-do the X Scaling
x_axis_correctscale = linspace(-400,0,length(All_TC_fit_norm_one));
x_axis_correctscale = x_axis_correctscale([ones(length(classes_only_cells),1)],:); % repeat that vector for all time courses (y values)

for i = 1 : length(unique_dates)
    for j = 1 : length(date_only_cells)
        if strcmp(date_only_cells{j}, unique_dates{i})
            that_retina = strcmp(full_description_cells(:), strcat(unique_dates{i}, '-OFFLBT'));
                       
            x_axis_correctscale(j,:) = x_axis_undone(j,:).*(All_RL(master_retina) / All_RL(that_retina));
            
        end
    end
end

plot_timecourses_all_retinas(x_axis_correctscale, All_TC_xynorm, unique_classes, unique_dates, classes_only_cells, date_only_cells, color)




%%%%%% re-sample using interpolation

spline_x_axis = linspace(-330,0,241); % chosen from min of x_axis_norm
spline_x_axis = spline_x_axis([ones(length(classes_only_cells),1)],:); % repeat that vector for all time courses (y values)

for i = 1 : length(date_only_cells);
    All_TC_spline_norm(i,:) = interp1(x_axis_correctscale(i,:), All_TC_xynorm(i,:), spline_x_axis(i,:));
end



plot_timecourses_all_retinas(spline_x_axis, All_TC_spline_norm, unique_classes, unique_dates, classes_only_cells, date_only_cells, color)



%%% Look at PCA space - either with clustering or just plotting by class?








end