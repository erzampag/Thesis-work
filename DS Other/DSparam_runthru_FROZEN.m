% Frozen 11/11/13 to attempt better organization 

% function [] = DSparam_runthru(textdata)
% This script takes a list of experimental runs (organized like
% DSnotes.xls), and either uses drifting_squarewave_shorter_stimulus2.m or
% saved data from that function (or DS_analysis_trial1.m which is all
% commented out at the bottom) to do analysis on all retinas together.

% Written by Erin Zampaglione in Summer 2013.



% for i  = 1 : length(textdata)
% for i = [18 19];
% for i = [18];
%     
%     % This is when you have saved data!
%     %     file = strcat('/Users/erinzampaglione/Documents/Lab_Work/DSCells/DS_saveddata_no_OS/'...
%                                 , textdata{i,1}, '/', strcat(textdata{i,2}, '.mat'));
%     file = strcat('/Users/erinzampaglione/Documents/Lab_Work/DSOSCells/', textdata{i,1}, '/', 'data002', '.mat');
%     
%     load(file)
%     
%     % % % % %     % This is when you have to use raw data:
%     % % % % %     [idList, AllCells, DSI, DSI_error, ratio, noiseratio,...
%     % % % % %         noise, zerofreq, fwhh, reduced_chi_square, MaxSpikes, EIx, EIy,...
%     % % % % %         Fone_Fzero, Ftwo_Fzero, OS_x1, OS_y1, OS_x2, OS_y2, Max_Spikes_absolute, spike_rate_all_orientation]=...
%     % % % % %         MultiBarsAnalysis2(textdata{i,1}, 'data002', 's02.txt');
%     
% end

% file = (['/Users/erinzampaglione/Documents/Lab_work/DS_param/' '2013-09-19-0/' 'data001.mat']);
% 
% load(file)

% keyboard

%%% try finding any cell with any DSI > 0.5

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


DS_indices = cell(length(idList),1);

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




DS_20130919_hand = ...
[33;74;85;104;105;107;123;127;147;159;168;198;234;236;254;...
    282;286;290;307;314;343;351;373;383;391;411;422;423;437;460;471;485;488];

mean_rates = cell(length(idList), 1);
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



% % DS_classtextfile(idList, DS_index,...
% %     ['/Users/erinzampaglione/Documents/processed_data/' '2013-08-14-0/' 'data000-map-data001/' 'DS_class.txt']);

% % DS_classtextfile(idList, DS_index,...
% %     ['/Users/erinzampaglione/Documents/processed_data/' '2012-05-03-0/' 'data001-map-data002/' 'DS_class.txt']);

% DS_classtextfile(idList, DS_index,...
%     ['/Users/erinzampaglione/Documents/processed_data/' '2013-09-19-0/' 'data000-map-data001/' 'DS_class.txt']);
% 
% DS_classtextfile(idList, DS_index,...
%     ['/Users/erinzampaglione/Documents/processed_data/' '2013-09-13-0/' 'data000-map-data001/' 'DS_class.txt']);
% 
% DS_classtextfile(idList, DS_20130919_hand,...
%     ['/Users/erinzampaglione/Documents/processed_data/' '2013-09-19-0/' 'data000-map-data001/' 'DS_HAND.txt']);
% 
% DS_classtextfile(idList, good_neurons_index,...
%     ['/Users/erinzampaglione/Documents/processed_data/' '2013-10-17-0/' 'data000-map-data001/' 'DS_HAND.txt']);
% 
% DS_classtextfile(idList, good_neurons_index,...
%     ['/Users/erinzampaglione/Documents/processed_data/' '2012-05-03-0/' 'data001-map-data002/' 'DS_HAND.txt']);
% 
% DS_classtextfile(idList, good_neurons_index,...
%     ['/Users/erinzampaglione/Documents/processed_data/' '2012-04-30-0/' 'data002/' 'DS_HAND.txt']);
% 
% DS_classtextfile(idList, good_neurons_index,...
%     ['/Users/erinzampaglione/Documents/processed_data/' '2013-10-17-0/' 'data001-map-data003/' 'DS_HAND.txt']);


keyboard 

%% let's create a neuron display that shows every SP/TEMP/ORIENTATION
good_neurons = [];
good_neurons_index = [];
mkdir(date);
% index= 1;
% for q = 1:length(index)
% for q = find(idList(:) == 3917)
% for q = DS_20130919_hand' 
for q = DS_index'
% for q = good_neurons_index

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
        
%         keyboard

    response = input('Is this a good cell?Y/N:', 's');
    
    if lower(response) == 'y'
        
        good_neurons = [good_neurons idList(q)];
            
        set(gcf, 'PaperPosition', [0,0,15,15.5]);
        print('-depsc', [date '/' 'Neuron' num2str(idList(q))])
        
    end
        close(gcf);
end


    good_neurons_index = zeros(1,length(good_neurons));
    for i = 1: length(good_neurons)
        good_neurons_index(i) = find(idList(:) == good_neurons(i));
    end
    save([date '/' 'good_neurons_index'], 'good_neurons_index');
    


keyboard




%% Try representing a retina in a similar manner to the neuron display, by only graphing vectors
%%% that came from finding DSI > 0.5 and MaxSpikes > 50

% remake X and Y's for compass plot
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
%                 
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
        print('-depsc', [date '/' 'All_DS_vectors']);
        close(gcf);

%% I want to make a 3D histogram with the count of how many neurons were DS in each SP/TEMP!
figure
h = bar3(hist_counter) %, 'FaceAlpha', 0.65);
colorbar

for k = 1:length(h)
    zdata = get(h(k), 'ZData');
    set(h(k), 'CData', zdata, 'FaceColor', 'interp')
end
set(gcf,'renderer', 'opengl')


% Actually, instead let's put each neuron in it's favorite SP/TEMP ONCE!!!
hist_all = zeros(length(idList),2);

for i  = 1 : length(idList)
    [spat_all, tem_all] = find(mean_rates{i} == max(max(mean_rates{i})));
    if length(tem_all) > 1
        fprintf('neuron with more than one max mean rate\n')
        
        continue % ignore all neurons that have more than one max spike rate
    end

    hist_all(i,1) = spat_all; hist_all(i,2)=  tem_all;
%     hist_all
end

figure
hist3(hist_all, [5,5])
figure
hist(hist_all(:,1))


for i  = good_neurons_index
    [spat_all, tem_all] = find(mean_rates{i} == max(max(mean_rates{i})));
    if length(tem_all) > 1
        fprintf('neuron with more than one max mean rate\n')
        
        continue % ignore all neurons that have more than one max spike rate
    end

    hist_all(i,1) = spat_all; hist_all(i,2)=  tem_all;
%     hist_all
end





keyboard

%% How bout for each neuron, show it's preferred direction / max rate PSTH
% % index = 1;
% for q = 1:length(index)
% for q = 444
% for q = find(idList(:) == 3917)
for q = good_neurons_index

% for q = DS_20130919_hand'
%     DS_indices{q}
%     keyboard
    
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
        print('-depsc', [date '/' 'PSTH' num2str(idList(q))])
        
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


%%
% attempt at surface plot

% % % % %     figure
% % % % %     %     subplot(1,3,1)
% % % % %     [X,Y] = meshgrid(log(spatial), log(temporal));
% % % % %     surf(Y, X, max_rates{index});
% % % % %     colormap(gray); colorbar
% % % % %     ylabel('Spatial periods'); xlabel('Temporal Periods'); zlabel('Spike Rate');

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
     print('-depsc', [date '/' 'Spatial' num2str(idList(q))])
     
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
     print('-depsc', [date '/' 'Temporal' num2str(idList(q))])

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
     print('-depsc', [date '/' 'Speed' num2str(idList(q))])
     
     close(gcf)
end


%%

% THREE-D HISTOGRAM WHERE EACH NEURON IS PLACED IN THE HIGHEST MEAN SP/TEMP
% BOX
hist_DS = zeros(length(good_neurons_index),2);
hist_counter = 0;
for i  = (good_neurons_index)

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
print('-depsc', [date '/' 'All_DS_S-T_hist'])


% % % % % % % % % % % % % % % 3D bar graph of the mean rates
% % % % % % % % % % % % % % X = 262;
% % % % % % % % % % % % % % figure
% % % % % % % % % % % % % % h = bar3(mean_rates{X}); %, 'FaceAlpha', 0.65);
% % % % % % % % % % % % % % colorbar
% % % % % % % % % % % % % % set(gca, 'XTickLabel', (temporal)); xlabel('Temporal Periods')
% % % % % % % % % % % % % % set(gca, 'YTickLabel', spatial); ylabel('Spatial Periods')
% % % % % % % % % % % % % % title(['Index ' num2str(X) ' Neuron' num2str(idList(X))])
% % % % % % % % % % % % % % for k = 1:length(h)
% % % % % % % % % % % % % %     zdata = get(h(k), 'ZData');
% % % % % % % % % % % % % %     set(h(k), 'CData', zdata, 'FaceColor', 'interp')
% % % % % % % % % % % % % % end
% % % % % % % % % % % % % % set(gcf,'renderer', 'opengl')
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % figure
% % % % % % % % % % % % % % surf(mean_rates{X})

% 
figure
hold on
for q = good_neurons_index
    
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

%% make the concatinated vector like in Matt's thesis

concat_no_orient = zeros(length(idList), length(spatial)*length(temporal));
concat_no_orient_norm = zeros(length(idList), length(spatial)*length(temporal));


for i = 1 : length(idList)
  concat_no_orient(i, :) = cat(2, mean_rates{i}(1,:), mean_rates{i}(2,:),...
      mean_rates{i}(3,:), mean_rates{i}(4,:), mean_rates{i}(5,:));
  
  concat_no_orient_norm(i,:) = concat_no_orient(i,:)/sum(concat_no_orient(i,:));
        
end

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
     print('-depsc', [date '/' 'ConCatVect' num2str(idList(q))])
     
     close(gcf)

end

% Concatinated vectors all on one graph, max slices on one graph
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
color_long = cat(1,color_long, color_long);

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
print('-depsc', [date '/' 'All_DS_ConCatVect'])
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
print('-depsc', [date '/' 'All_DS_SpatSlice2'])
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
print('-depsc', [date '/' 'All_DS_TempSlice2'])
close(gcf)

%% average over all spat/temps, to just have direction






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
    % % % % % % % % % % % %
    % % % % % % % % % % % %  [coeff, score] = princomp(new_tuning_matrix);
    % % % % % % % % % % % %
    % % % % % % % % % % % % %  figure
    % % % % % % % % % % % % subplot(3,1,3)
    % % % % % % % % % % % %  scatter(score(:,1),score(:,2));
    % % % % % % % % % % % %  title(strcat('PCA for OS Cells in Retina-',...
    % % % % % % % % % % % %         textdata{i,1}),'FontSize', 25);
    
    
    [~, tuning_graph_index] = show_classification_plots((1:length(good_neurons_index)),...
        new_tuning_matrix, zeros(length(good_neurons_index)), zeros(1,length(good_neurons_index))');
    
    figure % Time Courses
    for i = 1 : length(tuning_graph_index(:,1))
        index = logical(tuning_graph_index(i,:));
        subplot(3,2,i);
        plot(new_tuning_matrix(index,:)', 'b');
        title(strcat('nc',num2str(i)));
    end
    
    % Something I really need is to be able to make subclasses based on the
    % PCA clustering
    
    

