% function [] = DSparam_runthru(textdata)
% This script takes a list of experimental runs (organized like
% DSnotes.xls), and either uses drifting_squarewave_shorter_stimulus2.m or
% saved data from that function (or DS_analysis_trial1.m which is all
% commented out at the bottom) to do analysis on all retinas together.

% Written by Erin Zampaglione in Summer 2013.


date = '2012-04-30-0'; data_file = 'data002';
date = '2012-05-03-0'; data_file = 'data002';
date = '2012-05-07-0'; data_file = 'data002';
date = '2013-08-14-0'; data_file = 'data001';
% % % % % % % % % % % % date = '2013-08-20-0'; data_file = 'data002';
date = '2013-09-06-0'; data_file = 'data001';
date = '2013-09-13-0'; data_file = 'data001';
% % % % % % % % % % % % date = '2013-09-17-0'; data_file = 'data001';
date = '2013-09-19-0'; data_file = 'data001';
date = '2013-10-17-0'; data_file = 'data001';
date = '2013-10-17-0'; data_file = 'data001-map-data003';
date = '2013-11-06-0'; data_file = 'data001';
date = '2013-11-06-0'; data_file = 'data002';
date = '2014-01-09-0'; data_file = 'data002';
date = '2014-01-09-0'; data_file = 'data002-map-data001';


[file, DS_indices, index, DS_index, mean_rates, max_rates, high_index...
    color, x_vector, y_vector, mean_rates_spatial, mean_rates_spatial_norm,mean_rates_temporal...
   mean_rates_temporal_norm, max_slice_spatial, max_slice_spatial_norm, max_slice_temporal...
   max_slice_temporal_norm, concat_no_orient, concat_no_orient_norm, color_long, mean_rate_for_each_dir...
   concat_orient, concat_orient_norm] = calculate_DS_param_props(date, data_file);

load(file);
%% Create a list of indices that are the neurons chose as DS

file_id = ['/Users/erinzampaglione/Documents/Lab_Work/DS_param_output/' date '/' data_file];
mkdir(file_id);

if exist([file_id '/good_neurons_index.mat'], 'file')
    load([file_id '/good_neurons_index.mat'])
else
    [good_neurons_index] =  DSparam_manual_examination(1, date, data_file, spatial, temporal,...
    DS_index, idList, orientation, spike_rate_all_orientation, DS_indices);
end


DS_classtextfile(idList, good_neurons_index,...
    ['/Users/erinzampaglione/Documents/processed_data/' date '/' data_file '/DS_HAND_' datestr(now,30) '.txt']);
% [good_neurons_index] =  DSparam_manual_examination(1, date, data_file, spatial, temporal,...
%     DS_index, idList, orientation, spike_rate_all_orientation, DS_indices);


%% Only using single param sp/temp, or using group around that sp/temp
% looking at just original sp/temp (sp64, temp128)

% % % % [R, theta] = calculate_R_theta(spikerate, orientation);
% % % % [DSI, MaxSpikes, index_of_min, index_of_opp, pref_O, DSI_error]  = determine_prefO_by_DSvector(theta, orientation, spikerate, q);
for i = 1:length(idList)
    
    correct_square = zeros(9, length(orientation));
    counter = 0;
    for j = [2 3 4]
        for k = [3 4 5]
            counter = counter+1;
    correct_square(counter,:) = squeeze(spike_rate_all_orientation{i}(j,k,:));
        end
    end
    
        spikerate_64_128(i,:) = spike_rate_all_orientation{i}(3,4,:);
        spikerate_64_128_surround(i,:) = sum(correct_square,1);
% keyboard
end

% looking at the 9 sp/temps surrounding single param sp/temp
for i = 1:length(idList)

[R, theta64_128(i)] = calculate_R_theta(spikerate_64_128(i,:)', orientation);
% % % % % [DSI_64_128, MaxSpikes, index_of_min, index_of_opp, pref_O, DSI_error]  
[DSI_64_128(i), ~, index_of_min, index_of_opp, pref_O, DSI_error_64_128(i)]  =...
    determine_prefO_by_DSvector(theta64_128(i), orientation, spikerate_64_128(i,:)', i);

end

% % % DS_indices_64_128 = cell(length(idList),1);

good_neurons_index_64_128 = [];
for i  = 1 : length(idList)
%     find_output = find(DSI_64_128(i) > 0.5 & max_rate_absolute{i}(3,4) > 2); % Max rate for this SP/TEMP, any orientation is 2 sp/sec
%     good_neurons_index_64_128(i) = find_output;
    
%     if DSI_64_128(i) > 0.5 %&& max_rate_absolute{i}(3,4) > 2; % Max rate for this SP/TEMP, any orientation is 2 sp/sec
    if DSI{i}(3,4) > 0.5
    good_neurons_index_64_128 =[ good_neurons_index_64_128 i];    
    end
%     DS_indices_64_128{i} = [spat, tem];
end
% % % DS_neuronIDs = zeros(length(idList), 1);
% % % for i = 1 : length(idList)
% % %     DS_neuronIDs(i) = ~(length(DS_indices{i}(:,1)) < 1); % need at least 2 SP/TEMPS that fulfill requirement
% % % end
% % % index = (1: length(idList))';
% % % DS_index = index(logical(DS_neuronIDs));
% % % 
% % % [good_neurons_index_] =  DSparam_manual_examination(1, date, data_file, spatial, temporal,...
% % %     DS_index, idList, orientation, spike_rate_all_orientation, DS_indices);





% looking at all sp/temps including and surrounding that

for i = 1:length(idList)

[R, theta_64_128_sur(i)] = calculate_R_theta(spikerate_64_128_surround(i,:)', orientation);

[DSI_64_128_sur(i), MaxSpikes, index_of_min, index_of_opp, pref_O, DSI_error_sur(i)]  =...
    determine_prefO_by_DSvector(theta_64_128_sur(i), orientation, spikerate_64_128_surround(i,:)', i);

end
good_neurons_index_64_128_sur = [];
for i  = 1 : length(idList)
    if DSI_64_128_sur(i) > 0.5
    good_neurons_index_64_128_sur =[ good_neurons_index_64_128_sur i];    
    end
%     DS_indices_64_128{i} = [spat, tem];
end



keyboard 



%% From this point on, use "good_neurons_index"

%% 5x5 polar plot neuron display that shows every SP/TEMP/ORIENTATION

% index= 1;
saveme = 1; % change to 1 if you want to save this neuron
% for q = 1:length(index)
% for q = find(idList(:) == 3917)
% for q = DS_20130919_hand' 
% for q = DS_index'
for q = good_neurons_index
% for q = 15
% for q = good_neurons_index_64_128_sur
% for q = this_cluster_index'


%     if isempty(DS_indices{q}) || length(DS_indices{q}(:,1)) < 2 % now there has to be at least 2 instances where DSI > 0.5
%          continue
%     end
    
    figure
    counter = 1;
    for i = 1: length(spatial)
        for j = 1:length(temporal)
            subplot(5,5,counter)
%             subplot(6,5,counter)
            
            p = polar2([orientation(:)*pi/180'; orientation(1)],...
                [squeeze(spike_rate_all_orientation{index(q)}(i,j,:)); spike_rate_all_orientation{index(q)}(i,j,1)],...
                [0, ceil(max(max(max(spike_rate_all_orientation{index(q)}))))]);
            % all plots go from [0:highest spike rate for that neuron rounded up]
%             set(p, 'LineWidth', 3); % Maek all lines fat
            
            counter = counter +1;
            
            if i ~= 1 || j ~= 3
                title(['Spat: ', num2str(spatial(i)), '  Temp: ' num2str(temporal(j))])
            end
            
            for k = 1:length(DS_indices{q}(:,1)) % loop over all the sp/temps found to be DS
                if DS_indices{q}(k,1) == i  && DS_indices{q}(k,2) == j
                    set(p, 'LineWidth', 3); % if this sp/temp matches one that made first cut DS, thicken the line
%                     set(p, 'LineColor', 'b');
                  
                else
%                     set(p, 'LineColor', 'r')
                end
            end
        end
    end
    
% % %     subplot(6,5,26)
% % %      hold on
% % %      for i  =  1: length(temporal)
% % %          plot(spat_conv(spatial), mean_rates{index(q)}(:,i), '-o', 'Color', color(i,:));
% % %      end
% % %      xlabel('spatial freq (cyc/deg)'); ylabel('rate (spikes/sec)')

        mtit(['Index ' num2str(q) ' Neuron' num2str(idList(q))])
       
if saveme ==1
    set(gcf, 'PaperPosition', [0,0,15,15.5]);
    print('-depsc', [file_id '/' 'Neuron' num2str(idList(q))])
            close(gcf);

end

%         close(gcf);
end


keyboard


%% Try representing a retina in a similar manner to the neuron display, by only graphing vectors
%%% that came from finding DSI > 0.5 and MaxSpikes > 50


hist_counter = zeros(length(spatial), length(temporal));

figure
counter = 1;
for i = 1: length(spatial)
    for j = 1:length(temporal)
        subplot(5,5,counter) % make one plot at a time, for each sp/temp
        
        x_fake=[0 1 0 -1];
        y_fake=[1 0 -1 0];
        h_fake=compass(x_fake,y_fake, '--');
        hold on;
        
%         for q = 1 : length(idList) % go thru every neuron for each plot
%         for q = DS_20130919_hand' % go thru every neuron for each plot
        for q = good_neurons_index    

%             if isempty(DS_indices{q})|| length(DS_indices{q}(:,1)) < 2
%                 continue % if there were no DS responses found in that neuron, skip to the next neuron
%             end
            
            for k = 1:length(DS_indices{q}(:,1)) % go thru every sp/temp found to be DS for this neuron
                
                if DS_indices{q}(k,1) ~= i  || DS_indices{q}(k,2) ~= j
                    continue % if the DS cells in that neuron don't correspond to this sp/temp, skip to next DS sp/temp
                end
                
                hist_counter(i,j) = hist_counter(i,j) +1;
%                 testme=  testme+1;
                ch = compass(x_vector{q}(i,j), y_vector{q}(i,j),'-r'); % Plot DS vector
                
                if ratio{q}(i,j) < 1
                    set(ch, 'Color', 'b'); % hopefully this only colors them if they are ON or OFF
%                     fprintf('count\n')
                else
                    set(ch, 'Color', 'k'); % hopefully this only colors them if they are ON or OFF
                end
                
                xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
                set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
                set(ch,'LineWidth',1.25) % was 1.5
                %Removing the label
                % set(findall(gcf, 'String', '30','String','60', 'String', '180') ,'String', ' ')
                
            end
            
            ch = title(['Spatial: ', num2str(spatial(i)), '  Temporal: ' num2str(temporal(j))]);
            %             ch = title(strcat('Pref Dir/Mag for DS Cells in Retina-',...
            %                 textdata{i,1}),'FontSize', 10);
            %         delete(findall(gcf, 'type', 'text'));
            
            set(h_fake,'Visible','off')
                        
        end
        counter = counter + 1;
        
    end
end

        set(gcf, 'PaperPosition', [0,0,15,15.5]);
        print('-depsc', [file_id '/' 'All_DS_vectors']);
        close(gcf);

keyboard

%% How bout for each neuron, show it's preferred direction / max rate PSTH
% % index = 1;
% for q = 1:length(index)
% for q = 444
% for q = find(idList(:) == 3917)
for q = good_neurons_index
    
%     if isempty(DS_indices{q}) || length(DS_indices{q}(:,1)) < 2 % now there has to be at least 3 instances where DSI > 0.5
%         continue
%     end
    
    figure
    counter = 1;
    for i = 1: length(spatial)
        for j = 1:length(temporal)
            subplot(5,5,counter)
            
            ori_index = best_orientation{q}(i,j);
           
            h = plot(PSTH_bins{j}, PSTHs{q}{i,j,ori_index});
            
            if ratio{q}(i,j) < 1
                set(h, 'Color', 'r'); % hopefully this only colors them if they are ON or OFF
            end
            
            for k = 1:length(DS_indices{q}(:,1)) % loop over all the sp/temps found to be DS
                if DS_indices{q}(k,1) == i  && DS_indices{q}(k,2) == j
                    set(h, 'LineWidth', 2); % if this sp/temp matches one that was DS, thicken the line
                end
            end
            
            title(['Spat: ', num2str(spatial(i)), '  Temp: ' num2str(temporal(j)) 'Ori: ' num2str(orientation(ori_index))])

            counter = counter +1;
        end
    end
%     mtit(['Index ' num2str(q) ' Neuron' num2str(idList(q))])
        set(gcf, 'PaperPosition', [0,0,15,10]); 
        print('-depsc', [file_id '/' 'PSTH' num2str(idList(q))])
        
        close(gcf)

end

%% Represent neurons as 4D matrix?

for q = 1:length(index)
    
    if isempty(DS_indices{q}) || length(DS_indices{q}(:,1)) < 2 % now there has to be at least 3 instances where DSI > 0.5
        continue
    end
    figure
    hold on
    
    for i = 1:length(spatial)
        for j = 1:length(temporal)
            for k = 1: length(orientation)
                scatter3(spatial(i), temporal(j), orientation(k), 150, spike_rate_all_orientation{index(q)}(i,j,k),'filled', 'o')
                % scatter3(log(spatial(i)), log(temporal(j)), orientation(k), 150, spike_rate_all_orientation{index}(i,j,k),'filled', 'o')
            end
        end
    end
    set(gca, 'CLim', [0, max(max(max(spike_rate_all_orientation{index(q)})))])
    %    set(gca, 'XScale', 'log');   set(gca, 'YScale', 'log')
    ylim([10 max(temporal)]); xlim([10 max(spatial)]);
    set(gca, 'XScale', 'log');   set(gca, 'YScale', 'log')
    xlabel('spatial'); ylabel('temporal'); zlabel('orientation');
    title(['Index ' num2str(q) ' Neuron' num2str(idList(q))])
    
end

%% Let's put together a section where orientation isn't taken into account - Just functions of sp/temp/velocity
%%% Heat-map maker:
figure
imagesc([spat_conv(spatial(1)) spat_conv(spatial(5))],...
    [temp_conv(temporal(1)) temp_conv(temporal(5))], max_rates{index})

imagesc([spatial(1) spatial(5)], [temporal(1) temporal(5)], max_rates{index})

ylabel('Spatial periods'); xlabel('Temporal Periods'); title('Average Spike Rate over all Orientations');
colormap(gray); colorbar


%%
%%%%%% let's try plotting spike rate as a function of SPATIAL period,
%%%%%% colored by TEMPORAL frequency for one neuron
%  for q = 1: length(index)

% for q = DS_20130919_hand'
% for q = find(idList(:) == 3917)
% for q = 444
for q = good_neurons_index

%      if isempty(DS_indices{q}) || length(DS_indices{q}(:,1)) < 2 % now there has to be at least 3 instances where DSI > 0.5 %
%          continue
%      end
%          keyboard
     figure
     hold on
     for i  =  1: length(temporal)
         plot(spat_conv(spatial), mean_rates{index(q)}(:,i), '-o', 'Color', color(i,:));
         
     end
     
     plot(spat_conv(spatial), mean_rates_spatial(index(q),:), '-x', 'Color', 'k');

     xlabel('spatial freq (cyc/deg)'); ylabel('rate (spikes/sec)')
     % % temp_conv(temporal)
     legend('7.5 Hz', '3.75 Hz', '1.875 Hz', '0.9375 Hz', '0.4688 Hz')
     title(['Index ' num2str(q) ' Neuron' num2str(idList(q))])
     
     set(gcf, 'PaperPosition', [0,0,6,4]); 
     print('-depsc', [file_id '/' 'Spatial' num2str(idList(q))])
     
     close(gcf)
 end


%%%%%% let's try plotting spike rate as a function of TEMPORAL period,
%%%%%% colored by SPATIAL frequency for one neuron

%  for q = 1: length(index)
% for q = DS_20130919_hand'  
for q = good_neurons_index
     
%      if isempty(DS_indices{q}) || length(DS_indices{q}(:,1)) < 2 % now there has to be at least 3 instances where DSI > 0.5 %
%          continue
%      end
     %     q
     figure
     hold on
     for i  =  1: length(spatial)
         plot(temp_conv(temporal), mean_rates{index(q)}(i,:), '-o', 'Color', color(i,:));
%          plot(temp_conv(temporal), mean_rates{q}(i,:), '-o', 'Color', color(i,:));

     end
     
     plot(temp_conv(temporal), mean_rates_temporal(index(q),:), '-x', 'Color', 'k');
     
     xlabel('temp freq cyc/sec (Hz)'); ylabel('rate (spikes/sec)')
     % % temp_conv(temporal)
     legend('0.2153 cyc/deg', '0.1076 cyc/deg', '0.0538 cyc/deg', '0.0269 cyc/deg', '0.0135 cyc/deg')
     title(['Index ' num2str(q) ' Neuron' num2str(idList(q))])
%      
          set(gcf, 'PaperPosition', [0,0,6,4]); 
     print('-depsc', [file_id '/' 'Temporal' num2str(idList(q))])

     close(gcf)
 end

 
% % % the last this to try is plotting rate as a function of SPEED: some speeds
% % % were tried more than others, ie 1 pixel/frame

% for q = 1: length(index)
% for q = DS_20130919_hand'
% for q = find(idList(:) == 3917)   
for q = good_neurons_index   

%     if isempty(DS_indices{q}) || length(DS_indices{q}(:,1)) < 2 % now there has to be at least 3 instances where DSI > 0.5 %
%         continue
%     end
    
    figure
    hold on
    for i = 1 : length(spatial)
        for j = 1: length(temporal)
            
            scatter(spatial(i)/temporal(j), mean_rates{index(q)}(i,j), 50, color(j,:), 'fill')
        end
    end
    set(gca, 'XScale', 'log');
    xlabel('pixels/frame'); ylabel('rate (spikes/sec)')
    
    title(['Index ' num2str(q) ' Neuron' num2str(idList(q))])
    
         set(gcf, 'PaperPosition', [0,0,6,4]); 
     print('-depsc', [file_id '/' 'Speed' num2str(idList(q))])
     
     close(gcf)
end

%%
% THREE-D HISTOGRAM WHERE EACH NEURON IS PLACED IN THE HIGHEST MEAN SP/TEMP
% BOX
% hist_DS = zeros(length(good_neurons_index),2);
hist_DS = zeros(length(index),2);

hist_counter = 0;
% for i  = (good_neurons_index)
 for i  = (index') % all neurons?

    

    high_point = max(max(mean_rates{i}));
    [spat_ind, temp_ind] = find(mean_rates{i} == high_point(1));
    
    hist_counter = hist_counter +1;
    hist_DS(hist_counter,1) = spat_ind; 
    hist_DS(hist_counter,2)=  temp_ind;
%     
end

figure
hist3(hist_DS, 'Edges', {1:5, 1:5}, 'FaceAlpha', .65);
set(gcf, 'renderer', 'opengl')
xlabel('spatial periods'), ylabel('temporal periods'); zlabel('number of neurons'), 
title(['Retina ' date]);
view([70 30]);

set(gcf, 'PaperPosition', [0,0,6,5]);
print('-depsc', [file_id '/' 'All_DS_S-T_hist'])


% 
figure
hold on
for q = good_neurons_index
    q
    %         plot(temp_conv(temporal), mean_rates_temporal(index(q),:), '-x', 'Color', 'k');
    plot(temp_conv(temporal), mean_rates_temporal_norm(index(q),:), '-x', 'Color', 'k');
    xlabel('temp freq cyc/sec (Hz)'); ylabel('rate (spikes/sec)')
    title(['All DS cells average rate as function of temp for' date])
end

figure
hold on
for q = good_neurons_index   
    
    %         plot(spat_conv(spatial), mean_rates_spatial(index(q),:), '-x', 'Color', 'k');
    plot(spat_conv(spatial), mean_rates_spatial_norm(index(q),:), '-x', 'Color', 'k');
    xlabel('spatial freq (cyc/deg)'); ylabel('rate (spikes/sec)')
    title(['All DS cells average rate as function of spat for' date])
end

%% graphing concatinated vectors

for q = good_neurons_index
figure
% bar(concat_no_orient(q,:))

bar(concat_no_orient_norm(q,:))

hold on
% plot([5.5  5.5],[0 10], '--');plot([10.5  10.5],[0 10], '--');plot([15.5  15.5],[0 10], '--');
% plot([20.5  20.5],[0 10], '--');

plot([5.5  5.5],[0 .15], '--');plot([10.5  10.5],[0 .15], '--');plot([15.5  15.5],[0 .15], '--');
plot([20.5  20.5],[0 .15], '--');
        title(['Index ' num2str(q) ' Neuron' num2str(idList(q))])
        
                 set(gcf, 'PaperPosition', [0,0,6,4]); 
     print('-depsc', [file_id '/' 'ConCatVect' num2str(idList(q))])
     
     close(gcf)

end

%% Concatinated vectors all on one graph, max slices on one graph

figure
color_counter = 1;
for q = good_neurons_index
    plot(concat_no_orient_norm(q,:), 'Color', color_long(color_counter, :));
    color_counter = color_counter +1;
hold on
end
plot([5.5  5.5],[0 .15], '--');plot([10.5  10.5],[0 .15], '--');plot([15.5  15.5],[0 .15], '--');
plot([20.5  20.5],[0 .15], '--');
title(['Retina ' date]);

set(gcf, 'PaperPosition', [0,0,6,4]);
print('-depsc', [file_id '/' 'All_DS_ConCatVect'])
close(gcf)


% % % % % % % matrix = [];
% % % % % % % for q = good_neurons_index
% % % % % % %     matrix = cat(1, matrix, concat_no_orient_norm(q,:));
% % % % % % % end
% % % % % % % figure
% % % % % % % plot(matrix')



figure
color_counter = 1;
for q = good_neurons_index
    plot(spat_conv(spatial), max_slice_spatial(q,:), '-o', 'Color', color_long(color_counter, :))
    color_counter = color_counter +1;
hold on
end
title(['Retina ' date ' Spatial Slice']);
 xlabel('spatial freq (cyc/deg)'); ylabel('rate (spikes/sec)')
set(gcf, 'PaperPosition', [0,0,6,4]);
print('-depsc', [file_id '/' 'All_DS_SpatSlice2'])
close(gcf)


figure
color_counter = 1;
for q = good_neurons_index
    plot(temp_conv(temporal), max_slice_temporal(q,:),'-o',  'Color', color_long(color_counter, :))
    color_counter = color_counter +1;
hold on
end
title(['Retina ' date ' Temporal Slice']);
xlabel('temp freq cyc/sec (Hz)'); ylabel('rate (spikes/sec)')
set(gcf, 'PaperPosition', [0,0,6,4]);
print('-depsc', [file_id '/' 'All_DS_TempSlice2'])
close(gcf)


%% Concatinated vectors WITH DIRECTIONAL information
% for q = 1:10
for q = good_neurons_index

figure
bar(concat_orient_norm(q,:))
hold on
% separators by spatial period
plot([40.5  40.5],[0 .07], '--');plot([80.5  80.5],[0 .07], '--');plot([120.5  120.5],[0 .07], '--');
plot([160.5  160.5],[0 .07], '--');

% separators by temporal period
for i = 8.5 : 8 : 192.5
    plot([i i], [0 .05], '--')
end
        title(['Index ' num2str(q) ' Neuron' num2str(idList(q))])
        xlim([0 201])
% % %                  set(gcf, 'PaperPosition', [0,0,6,4]); 
% % %      print('-depsc', [file_id '/' 'ConCatVect' num2str(idList(q))])
% % %      
%      close(gcf)

end

% PCA with the whole turning curve - maybe we can get out DS cells?
        new_tuning_matrix = [];
    for j = 1: length(index)
        % % % % for j = 1: length(AllCells)
        new_tuning_matrix(j, :) = concat_orient_norm(index(j),:);
        % % %         new_tuning_matrix(j, :) = norm_spike_rate{j};
    end
    
    [~, tuning_graph_index] = show_classification_plots((1:length(new_tuning_matrix)),...
        new_tuning_matrix, zeros(length(new_tuning_matrix)), zeros(1,length(new_tuning_matrix))',...
        zeros(1,length(new_tuning_matrix))');

    
    figure % Concatinated vectors clustered by PCA
    for i = 1 : length(tuning_graph_index(:,1))
        index = logical(tuning_graph_index(i,:));
        subplot(3,3,i);
        plot(new_tuning_matrix(index,:)', 'b');
        hold on
        plot([40.5  40.5],[0 .07], '--');plot([80.5  80.5],[0 .07], '--');plot([120.5  120.5],[0 .07], '--');
        plot([160.5  160.5],[0 .07], '--');
        title(strcat('nc',num2str(i)));
    end
    
%% I should really try the FFT I bet with so many trials I might see
% the peak for DS cells
for q = 1:32
% for q = good_neurons_index
    Y=fft(concat_orient_norm(q,:));
    n=length(Y);
    % % if n>1;
    % %     Y(1)=[];
    % % end
    % % power = abs(Y(1:floor(n/2))).^2
    power_tuning = abs(Y(1:ceil(n/2)+1)).^2;
    
    amp_tuning = abs(Y(1:ceil(n/2)+1));
    
    power = abs(Y(1:floor(n/2))).^2;
nyquist = 1/2;
    
    figure
    subplot(2,1,1)
    plot(power_tuning(2:end))
    title(['Index ' num2str(q) ' Neuron' num2str(idList(q))])
    
subplot(2,1,2)
plot(concat_orient(q,:));
end


%% Here we can graph just the whole neuron's response at the orientation

for q = good_neurons_index
    figure
%     plot(mean_rate_for_each_dir(q,:))
polar([orientation orientation(1)]*pi/180, [mean_rate_for_each_dir(q,:) mean_rate_for_each_dir(q,1)])


% hold on
end
title(['Retina ' date ' Temporal Slice']);
xlabel('temp freq cyc/sec (Hz)'); ylabel('rate (spikes/sec)')
set(gcf, 'PaperPosition', [0,0,6,4]);
print('-depsc', [file_id '/' 'All_DS_TempSlice2'])
close(gcf)





%% maybe here we can put in the PCA?


    new_tuning_matrix = [];
    for j = 1: length(good_neurons_index)
        % % % % for j = 1: length(AllCells)
        new_tuning_matrix(j, :) = concat_no_orient_norm(good_neurons_index(j),:);
        % % %         new_tuning_matrix(j, :) = norm_spike_rate{j};
    end
    
    
        new_tuning_matrix = [];
    for j = 1: length(index)
        % % % % for j = 1: length(AllCells)
        new_tuning_matrix(j, :) = concat_no_orient_norm(index(j),:);
        % % %         new_tuning_matrix(j, :) = norm_spike_rate{j};
    end
    
    
    [~, tuning_graph_index] = show_classification_plots((1:length(new_tuning_matrix)),...
        new_tuning_matrix, zeros(length(new_tuning_matrix)), zeros(1,length(new_tuning_matrix))',...
        zeros(1,length(new_tuning_matrix))');
    
    figure % Concatinated vectors clustered by PCA
    for i = 1 : length(tuning_graph_index(:,1))
        index = logical(tuning_graph_index(i,:));
        subplot(3,3,i);
        plot(new_tuning_matrix(index,:)', 'b');
        hold on
        plot([5.5  5.5],[0 .15], '--');plot([10.5  10.5],[0 .15], '--');plot([15.5  15.5],[0 .15], '--');
        plot([20.5  20.5],[0 .15], '--');
        title(strcat('nc',num2str(i)));
    end
    
    % Something I really need is to be able to make subclasses based on the
    % PCA clustering
    
    %% I want to plot ALL DS neurons from one retina on one cartesian graph,
    % using matlab vertical indexing to color them by their MAX spike rate!
    
    figure 
    hold on
    for q = good_neurons_index
      vert_index = (high_index(q,1)-1)+(high_index(q,2)-1)*5 +1;

        plot(orientation, mean_rate_for_each_dir(q,:), 'Color', color_long(vert_index,:))
          
    end
    
    %% Work on things that can be used on all retinas
    
    %% comparing two DSIs with different spatial and temporal periods
    DSI_comparison = zeros(length(idList), 2);
    vector_dist_2_4 = zeros(length(idList),1);
    vector_dist_4_2 = zeros(length(idList),1);
    
    figure
    counter = 1;
    for i = 1:length(idList)
        DSI_comparison(i,:) = [DSI{i}(2,4) DSI{i}(4,2)];
        
        vector_dist_2_4(i) = sqrt(x_vector{i}(2,4)^2 + y_vector{i}(2,4)^2);
        vector_dist_4_2(i) = sqrt(x_vector{i}(4,2)^2 + y_vector{i}(4,2)^2);

%         scatter(vector_dist_2_4,vector_dist_4_2, 50, 'k')
        
%         scatter(DSI{i}(2,4), DSI{i}(4,2), 40, 'k')
        
        if counter > length(good_neurons_index)
            continue
        end
        if i == good_neurons_index(counter)
%             scatter(vector_dist_2_4(i),vector_dist_4_2(i), 50, 'b', 'filled')
           
%             scatter(DSI{i}(2,4), DSI{i}(4,2), 30, 'b')
            
            counter = counter + 1;
        end
        hold on
    end
    
    
            scatter(vector_dist_2_4,vector_dist_4_2, 50, [0.5 0.5 0.5])
hold on
for q = good_neurons_index
                scatter(vector_dist_2_4(q),vector_dist_4_2(q), 50, 'r', 'filled')

end
    
    xlabel(['Spatial: ' num2str(spatial(2)) ' Temporal: ' num2str(temporal(4))],'FontSize', 25);
    ylabel(['Spatial: ' num2str(spatial(4)) ' Temporal: ' num2str(temporal(2))],'FontSize', 25);
    title(['Comparing vectorlength in retina' date])
    
    set(gcf, 'PaperPosition', [0,0,12,8]);
    print('-depsc', [file_id '/' 'Vector_Length_Comp'])
%     close(gcf)
    
    
        [~, tuning_graph_index] = show_classification_plots((1:length(idList)),...
        zeros(length(idList)), zeros(length(idList)), vector_dist_2_4,...
        vector_dist_4_2);
    
    this_cluster_index = index(logical(tuning_graph_index(1,:)));
    
 %% how about this: for each neuron find the max spatial frequency (color code these)
 % and plot spike rate as a function of temporal frequency, then do the
 % same for max temporal frequency (color code) and plot rate as a function
 % of spatial frequency.
 
 
 %% use the average spike rate along spatial periods, the average spike rate along
 % temporal periods, and use those as weights for the spatial / temporal
 % periods (or log of that) then plot as a scatter plot, bold if they were
 % identified as DS - a way to see the whole retina's spatial/temporal
 % tuning
 
 weighted_spatial = zeros(1,length(idList));
 weighted_temporal = zeros(1,length(idList));
 for q = 1:length(idList)
% for q = good_neurons_index

%      weighted_spatial(q) = mean_rates_spatial_norm(q,:)*spatial';
%      weighted_temporal(q) = mean_rates_temporal_norm(q,:)*temporal';

          weighted_spatial(q) = mean_rates_spatial_norm(q,:)*log10(spat_conv(spatial))';
     weighted_temporal(q) = mean_rates_temporal_norm(q,:)*log10(temp_conv(temporal))';
     
 end
    
 figure
scatter(weighted_temporal, weighted_spatial, 50, [0.5 0.5 0.5])
hold on
for q = good_neurons_index
    
    scatter(weighted_temporal(q), weighted_spatial(q), 50, 'r', 'filled')
end
xlabel('log_1_0 temporal freq(Hz)', 'FontSize', 25), ylabel('log_1_0 spatial freq(cyc/deg)', 'FontSize', 25)
title(['Weighted spat/temp ' date])

    set(gcf, 'PaperPosition', [0,0,12,8]);
    print('-depsc', [file_id '/' 'Weighted_Spat_Temp'])

    
    [~, tuning_graph_index] = show_classification_plots((1:length(idList)),...
        zeros(length(idList)), zeros(length(idList)), weighted_spatial,...
        weighted_temporal);
    
[] = PCA_classtextfile(date, data_file, tuning_graph_index, idList);

