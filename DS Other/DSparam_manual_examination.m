function [good_neurons_index] =  DSparam_manual_examination(manual_input, date, data_file, spatial, temporal,...
    DS_index, idList, orientation, spike_rate_all_orientation, DS_indices)

% if input = 1, do manual examination
% if input  = 0, just print all 

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
                [squeeze(spike_rate_all_orientation{(q)}(i,j,:)); spike_rate_all_orientation{(q)}(i,j,1)],...
                [0, ceil(max(max(max(spike_rate_all_orientation{(q)}))))]);
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
    
    
    
    if manual_input == 1
        response = input('Is this a good cell?Y/N:', 's');
        
        if lower(response) == 'y'
            
            good_neurons = [good_neurons idList(q)];
            
            set(gcf, 'PaperPosition', [0,0,15,15.5]);
            print('-depsc', [date '/' 'Neuron' num2str(idList(q))])
            
        end
        
        close(gcf);
    end
end

if manual_input == 1
    good_neurons_index = zeros(1,length(good_neurons));
    for i = 1: length(good_neurons)
        good_neurons_index(i) = find(idList(:) == good_neurons(i));
    end
    save([date '/' data_file 'good_neurons_index'], 'good_neurons_index');
end


end
