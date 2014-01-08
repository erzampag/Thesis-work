% for i = 1 : 23 
% for i = 17:23
% for i = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23];
for i = 5;
    
      %     % This is when you have saved data!
%     
        file = strcat('/Users/erinzampaglione/Documents/Lab_Work/SbC_cells/', textdata{i,1}, '/', strcat(textdata{i,2}, '.mat'));
        load(file)
    
%     keyboard
    %% This is when you have to do actual analysis on raw data - spikes, etc
%     
% [idList, spike_rate_all_orientation, max_true_spike_rate, min_true_spike_rate, contrast_spike_rate, ...
%         grey_spike_rate] = grey_screen(textdata{i,1},textdata{i,2}, textdata{i,6});
    


    figure
    if strncmp(textdata{i,4}, 'W', 1)
        scatter(contrast_spike_rate, grey_spike_rate);
        hold on;
        
    elseif strncmp(textdata{i,4}, 'H', 1)
        scatter(contrast_spike_rate, grey_spike_rate, 'go');
        hold on;
        
    elseif strncmp(textdata{i,4}, 'D', 1)
        scatter(contrast_spike_rate, grey_spike_rate, 'rx');
        hold on;

    end
    plot((0:0.1:50),(0:0.1:50),'k-'); plot((0:0.1:50), 5*ones(1,length(0:0.1:50)), 'k-'); plot(5*ones(1,length(0:0.1:50)),(0:0.1:50), 'k-');

    xlabel('Average spike rate for all orientations','FontSize', 25); ylabel('Average spike rate at grey screen','FontSize', 25);
    title(strcat('Cells from Retina-', textdata{i,1}), 'FontSize', 25);
    
    cellsofinterest = find(grey_spike_rate(:) > 5 & contrast_spike_rate(:)< 5);
    
end
keyboard