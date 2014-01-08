function [nc, master_graph_index] = show_classification_plots(ON_index, All_TC, All_ACF, All_norm_RF,NEW_VAR)
% This is just the classification plots portion of PCA select (split off
% into it's own function so I can do the ON/OFF division first).  ON_index
% is called that way because it was first just used as a script looking at
% the ON cells only.
%
% ENZ, Fall 2012

counter = 0; % how many times over the while loop
nc = true(1,length(ON_index)); % array of everything outside the rectangle
reply = 'n';
% % % chosen_plot = 'string';
    figure;

while lower(reply) ~= 'y'
% % % % while ~strcmp(chosen_plot, 'done')
%% Calculuate PCA
    
    % %     if counter == 0 %% | in == false(1,length(ON_index))  % no longer
    % using this to index into anything!
    % %         in = true(1,length(ON_index)); %
    % %     end
    
    correct_index = nc(counter+1,:); % use to index into the things you want to PCA (TC, ACF) or plot (RF)
    
    [ACF_ON_coeff,ACF_ON_score] = princomp(All_ACF(ON_index(correct_index),:));
    [TC_ON_coeff,TC_ON_score] = princomp(All_TC(ON_index(correct_index),:));
    
%         TC_ON_score
    %% Vision Plots

    %     clf(handle);
    %     figure(handle);
    %     classification_plots = figure; % ON classes
    figure
    plot1 = subplot(2,3,1);
    % gscatter(TC_ON_score(:,1),TC_ON_score(:,2),classes_only_cells(ON_index))
    scatter(TC_ON_score(:,1),TC_ON_score(:,2));
    xlabel('TF1', 'FontSize', 20); ylabel('TF2', 'FontSize', 20)
    
    plot2 = subplot(2,3,2); scatter(All_norm_RF(ON_index(correct_index)), TC_ON_score(:,1));
    xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF1', 'FontSize', 20);
    title('Vision PCA plots for ON classes')
    
    plot3 = subplot(2,3,4); scatter(ACF_ON_score(:,1),ACF_ON_score(:,2));
    xlabel('Autocorrelation 1', 'FontSize', 20); ylabel('Autocorrelation 2', 'FontSize', 20);
    
    plot4 = subplot(2,3,5); scatter(All_norm_RF(ON_index(correct_index)),ACF_ON_score(:,1));
    xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
    
    plot5 = subplot(2,3,6); scatter(TC_ON_score(:,1),ACF_ON_score(:,1))
    xlabel('TF1', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
    
    
    plot6 = subplot(2,3,3); scatter(All_norm_RF(ON_index(correct_index)),NEW_VAR(ON_index(correct_index)))
    xlabel('VAR', 'FontSize', 20); ylabel('VAR', 'FontSize', 20);
    
%     linkdata(plot1); linkdata(plot2); linkdata(plot3); linkdata(plot4); linkdata(plot5);
% linkdata('on');
    
    % gname(full_description(ON_index,:))
    % gname(classes_only(ON_index,:))
    
    %% Cluster
    chosen_plot = [];
    dummyvar = 1234;
    while isempty(chosen_plot)
        chosen_plot = input('which plot do you want? "plot_":');

        if isempty(chosen_plot)
            chosen_plot = dummyvar;
        end
        
        if chosen_plot == plot1 || chosen_plot == plot2 || chosen_plot == plot3...
                || chosen_plot == plot4 || chosen_plot == plot5 || chosen_plot == plot6
            break
        else
            fprintf('not a correct graph handle\n')
            chosen_plot = [];
        end
    end
    
   
    
    % User selects rectangle, and needs to double click on rectangle when done positioning (from drp_arash.m)
    h=imrect(chosen_plot);
    position = wait(h);
    hold on
    % Get position of rectangle in [x_min, y_min, width, height]
    pos=getPosition(h);
    
    x_min = pos(1,1);
    y_min = pos(1,2);
    width = pos(1,3);
    height = pos(1,4);
    x_max = x_min+width;
    y_max = y_min+height;
    
    % Draws the ROI
    rectangle('position', [x_min, y_min, width, height]);
    % Specifies the vertices of the ROI
    xv = [x_min  x_min  x_max  x_max  x_min];
    yv = [y_min  y_max  y_max  y_min  y_min];
    
    % % %     keyboard
    
    % Detect neurons within the ROI.
    if chosen_plot == plot1
        in = inpolygon(TC_ON_score(:,1),TC_ON_score(:,2),xv,yv);
    elseif chosen_plot == plot2
        in = inpolygon(All_norm_RF(ON_index(correct_index)),TC_ON_score(:,1),xv,yv);
    elseif chosen_plot == plot3
        in = inpolygon(ACF_ON_score(:,1),ACF_ON_score(:,2),xv,yv);
    elseif chosen_plot == plot4
        in = inpolygon(All_norm_RF(ON_index(correct_index)),ACF_ON_score(:,1),xv,yv);
    elseif chosen_plot == plot5
        in = inpolygon(TC_ON_score(:,1),ACF_ON_score(:,1),xv,yv);
    elseif chosen_plot == plot6
        in = inpolygon(All_norm_RF(ON_index(correct_index)),NEW_VAR(ON_index(correct_index)),xv,yv);
    end
    
    % % %     graph_index = in'
    
    in = ~in'; % Convert logical from inside rectangle to outside rectangle (bad naming scheme)
    
    counter = counter + 1; % how many times over the while loop
    
    nc_counter = 0; % should be how many trues in new logical
    
    for i = 1: length(nc) % Expand new logical to correspond to the trues in the previous array "logical array"
        if nc(counter,i)
            nc_counter = nc_counter + 1;
            nc(counter+1,i) = in(nc_counter);
        end
    end
    
    
    
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
    

    figure % Time Courses
for i = 1 : length(graph_index(:,1))
    index = logical(graph_index(i,:));
    subplot(5,3,i);
    plot(All_TC(ON_index(index),:)', 'b');
    title(strcat('nc',num2str(i)));
    
        ylim([-0.7 0.7]);
    hold on
    plot(zeros(1,25),':');
    
end

figure % Autocorrelation Function
for i = 1 : length(graph_index(:,1))
    index = logical(graph_index(i,:));
    subplot(5,3,i);
    plot(All_ACF(ON_index(index),:)', 'b');
    title(strcat('nc',num2str(i)));

end

% % % % 
% % % % subclass = input('would you like to subclassify? Y/N:', 's');
% % % % if subclass == 'y'
% % % %     
% % % %     [~, sub_graph_index] = show_classification_plots(ON_index(index), All_TC, All_ACF, All_norm_RF)
% % % %     keyboard
% % % %     
% % % % %         sub_graph_index = []; %  a logical array corresponding to each selection, (as for graphing) for ease of putting back
% % % % % for i = 1: length(sub_nc)
% % % % %     for j = 2:length(sub_nc(:,i))
% % % % %         if sub_nc(j,i) ~= sub_nc(j-1,i)
% % % % %             sub_graph_index(j-1,i) = true;
% % % % %         else
% % % % %             sub_graph_index(j-1,i) = false;
% % % % %         end
% % % % %     end
% % % % % end
% % % %     
% % % % put_back = input('would you like to put neurons back in the higher cluster? 0, 1, 2, etc:');    
% % % % 
% % % % if put_back == 0
% % % %     
% % % % else
% % % % 
% % % %     sub_graph_counter = 0;
% % % %     for i = 1 : length(nc)
% % % %        if ~nc(counter+1, i) % editing most recent nc, looking at spaces (neurons that were in a cluster)
% % % %         sub_graph_counter = sub_graph_counter +1;
% % % %         nc(counter+1,i) = sub_graph_index(put_back,sub_graph_counter); % if it was in the graph, put it back in the classification pool
% % % %        end
% % % %     end
% % % %     
% % % % end
% % % % 
% % % %     keyboard
% % % %      
% % % %     
% % % % end
    
    
%     nc
    reply = input('Are you finished clustering? Y/N:', 's');
    if isempty(reply)
        reply = 'n';
    end
end


keyboard
% linkdata('off');

master_graph_index = []; % Once clustering is done, use the "outside selection" array to make a logical array corresponding to each selection
for i = 1: length(nc)
    for j = 2:length(nc(:,i))
        if nc(j,i) ~= nc(j-1,i)
            master_graph_index(j-1,i) = true;
        else
            master_graph_index(j-1,i) = false;
        end
    end
end
 master_graph_index

    
end