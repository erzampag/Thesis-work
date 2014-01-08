% This script takes a list of experimental runs (organized like
% DSnotes.xls), and either uses drifting_squarewave_shorter_stimulus2.m or
% saved data from that function (or DS_analysis_trial1.m which is all
% commented out at the bottom) to do analysis on all retinas together.

% Written by Erin Zampaglione in Fall 2012.

num_healthy_WT = [];
num_OS_WT = [];
num_DS_WT = [];

num_healthy_DSCAM = [];
num_OS_DSCAM = [];
num_DS_DSCAM = [];

%turn to 1 if you want to do PCA
PCA = 1;

% for i = 1 : length(textdata)
% for i = 4 : length(textdata)
% for i = 9;
% for i = 7;
% for i = [1 5 17 18 19 20 21 22 23]; % Runs that have EIx0 and EIy0 data
% for i = [1 2 4 5 6 7 8 9 10 11 12 13 14 15]; % Runs with OS variablesO
% for i = [1 2 4 5 6 7 8 9 10 11 12 13 14 15 16 18 19] % New runs with EI data and OS
% for i = [1 2 4 5 6 7 8 9 10 11 12 13 14 15 16 18 19 32 33] % New runs with EI data and OS
% for i = [2 9];
% for i = [4 5];
% for i = 9 % 2011-06- 28
    % for i = 15
% for i = 5; %2011-02-23
% for i = 4; %2011-02-17
% for i = 2; %2010-08-27
% for i = 11
% for i = 15; %2012-01-31
% for i = 10; %2011-06-30
for i = 1; %2010-08-25
        
    %% This is when you have saved data!
    
    file = strcat('/Users/erinzampaglione/Documents/Lab_Work/DSOSCells/', textdata{i,1}, '/', strcat(textdata{i,2}, '.mat'));
    load(file)
    
    %% This is when you have to do actual analysis on raw data - spikes, etc
    
% % %         [idList, AllCells, DSI, DSI_error, ratio, noiseratio,...
% % %             DScellsindex, noise, zerofreq, fwhh, reduced_chi_square, MaxSpikes, EIx, EIy,...
% % %             Fone_Fzero, Ftwo_Fzero, OS_x1, OS_y1, OS_x2, OS_y2, Max_Spikes_absolute, spike_rate_all_orientation]=...
% % %             MultiBarsAnalysis(textdata{i,1},textdata{i,2}, textdata{i,6});
    
    %%
    Healthy_cells = find(noiseratio(:) > 5);
    
    OScellsindex = find(Ftwo_Fzero(:) > sqrt(0.03) & Ftwo_Fzero(:) > Fone_Fzero(:) & Max_Spikes_absolute(:) > 100); % & noiseratio(:) > 5);
    DScellsindex = find(ratio(:) > 1 & DSI(:) > 0.5 & MaxSpikes(:) > 100 & noiseratio(:) > 5);
    
    if strncmp(textdata{i,4}, 'W', 1) %|| strncmp(textdata{i,4}, 'H', 1)
        num_healthy_WT = [num_healthy_WT length(Healthy_cells)]; % # of strongly responding cells
        num_OS_WT = [num_OS_WT length(OScellsindex)];
        num_DS_WT = [num_DS_WT length(DScellsindex)];
    end
    if strncmp(textdata{i,4}, 'D', 1)
        num_healthy_DSCAM = [num_healthy_DSCAM length(Healthy_cells)]; % # of strongly responding cells
        num_OS_DSCAM = [num_OS_DSCAM length(OScellsindex)];
        num_DS_DSCAM = [num_DS_DSCAM length(DScellsindex)];
    end
    
    norm_spike_rate = {};
    for j = 1: length(idList)
    norm_spike_rate{j} = spike_rate_all_orientation{j}./sum(spike_rate_all_orientation{j});
    end

%     keyboard
    %% Retina Graphs
% % % % % %     figure
% % % % % %     
% % % % % %    %     subplot(2,1,1)
% % % % % % %     subplot(2,1,1)
% % % % % %     x_fake=[0 1 0 -1];
% % % % % %     y_fake=[1 0 -1 0];
% % % % % %     
% % % % % %     h_fake=compass(x_fake,y_fake, '--');
% % % % % %     hold on;
% % % % % %     set(h_fake,'Visible','off');
% % % % % %     
% % % % % %     for j = 1 : length(OScellsindex)
% % % % % %         ch = compass(OS_x1(OScellsindex(j)), OS_y1(OScellsindex(j)), '-k');
% % % % % %         xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
% % % % % %         set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
% % % % % %         set(ch,'linewidth',1.5);
% % % % % %         
% % % % % %         if strncmp(textdata{i,4}, 'D', 1),
% % % % % %             set(ch, 'Color', 'r')
% % % % % %         end
% % % % % %         
% % % % % %         ch = compass(OS_x2(OScellsindex(j)), OS_y2(OScellsindex(j)), '-k');
% % % % % %         xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
% % % % % %         set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
% % % % % %         set(ch,'linewidth',1.5);
% % % % % %         
% % % % % %         
% % % % % %         
% % % % % %         if strncmp(textdata{i,4}, 'D', 1),
% % % % % %             set(ch, 'Color', 'r')
% % % % % %         end
% % % % % %     end
% % % % % %     
% % % % % %     ch = title(strcat('Vectors for OS Cells in Retina-',...
% % % % % %         textdata{i,1}),'FontSize', 25);
    
    
    %% Tuning curves
   
    
    figure

%     subplot(2,1,2)
%     subplot(2,1,2)
    hold on
    for j = 1 : length(OScellsindex)
        
%         h = plot(spike_rate_all_orientation{OScellsindex(j)});
        h = plot((0:22.5:337.5), norm_spike_rate{OScellsindex(j)});
        
        if strncmp(textdata{i,4}, 'D', 1),
            set(h, 'Color', 'r')
        end
    end
    title(strcat('Tunings for OS Cells in Retina-',...
        textdata{i,1}),'FontSize', 25);
    
%     keyboard
    %% EI Plotting (not many have this information)        

% Offsets taken from Georges's plotPositionEI.m

% % % % posElectrodes = dlmread('/Users/erinzampaglione/Documents/Lab_Work/matlabscriptsforactivationmaps/512coords.txt');
% % % % mapXOffset = 1042;
% % % % mapYOffset = 845;
% % % % figure
% % % % scatter((posElectrodes(:,1) + mapXOffset)*2,(posElectrodes(:,2) + mapYOffset)*2,'.k');
% % % % axis equal
% % % % hold on
% % % % plot((EIx+mapXOffset)*2, (EIy+ mapYOffset)*2, 'og')
% % % %         for j = 1: length(OScellsindex)
% % % %             plot((EIx(OScellsindex(j))+mapXOffset)*2, (EIy(OScellsindex(j))+ mapYOffset)*2, 'or', 'MarkerFaceColor', 'r')
% % % %         end
% % % %         
% % % % %         for j = 1 : length(DScellsindex_ONorOFF)
% % % % %             plot((EIx(DScellsindex_ONorOFF(j))+mapXOffset)*2, (EIy(DScellsindex_ONorOFF(j))+ mapYOffset)*2, 'ob', 'MarkerFaceColor', 'b')
% % % % %         end        
% % % %         
% % % %         title(strcat('EI postition for OS Cells in Retina-', textdata{i,1}),'FontSize', 25);
    
% keyboard

%% Ratio
figure
hist(log10(ratio(OScellsindex)), (-2.5:0.2:1.5))
    title(strcat('Log(F2/F1) Hist for OS Cells in ',...
        textdata{i,1}),'FontSize', 25);
    xlim([-1.5,1.5])

%  keyboard
idList(OScellsindex)

%% PCA for tuning curve (taken from Vision-like PCA script)

if PCA == 1;
    new_tuning_matrix = [];
    for j = 1: length(OScellsindex)
        % % % % for j = 1: length(AllCells)
        new_tuning_matrix(j, :) = norm_spike_rate{OScellsindex(j)};
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
    
    
    [~, tuning_graph_index] = show_classification_plots((1: length(OScellsindex)),...
        new_tuning_matrix, zeros(length(OScellsindex)), zeros(1,length(OScellsindex))');
    
    tuning_graph_index
    
    figure % Time Courses
    for i = 1 : length(tuning_graph_index(:,1))
        index = logical(tuning_graph_index(i,:));
        subplot(3,2,i);
        plot((0:22.5:337.5), new_tuning_matrix(index,:)', 'b');
        title(strcat('nc',num2str(i)));
    end
    
     
    
    idList(OScellsindex(logical(tuning_graph_index(1,:))))
    keyboard

end
end
keyboard
%% Boxplot OS cells and DS cells
percent_OS_WT = (num_OS_WT./num_healthy_WT).*100;
percent_DS_WT = (num_DS_WT./num_healthy_WT).*100;

percent_OS_DSCAM = (num_OS_DSCAM./num_healthy_DSCAM).*100;
percent_DS_DSCAM = (num_DS_DSCAM./num_healthy_DSCAM).*100;

all_things = cat(2, percent_OS_WT, percent_DS_WT, percent_OS_DSCAM, percent_DS_DSCAM);

box_string = {};
% for i = 1 : length(percent_OS_WT)
%     box_string = strcat(box_string, 'OS WT');
% end
% for i = 1 : length(percent_DS_WT)
%     box_string = strcat(box_string, 'DS WT');
% end
% for i = 1 : length(percent_OS_DSCAM)
%     box_string = strcat(box_string, 'OS DSCAM');
% end
% for i = 1 : length(percent_DS_DSCAM)
%     box_string = strcat(box_string, 'DS DSCAM');
% end
counter = 0;
for i = 1 : length(percent_OS_WT)
    box_string(i) = {'OS WT'};
    counter = counter+1;
end
counter = counter +1;
for i = counter : counter + length(percent_DS_WT) -1
    box_string(i) = {'DS WT'};
    counter = counter +1;
end
% counter = counter +1;
for i = counter : counter + length(percent_OS_DSCAM) -1
    box_string(i) = {'OS DSCAM'};
    counter = counter +1;
end
% counter = counter +1;
for i = counter : counter + length(percent_DS_DSCAM) -1
    box_string(i) = {'DS DSCAM'};
end

keyboard

figure
bp = boxplot(all_things, box_string, 'whisker', 6, 'widths', 0.4, 'Color', 'k');
set(bp, 'linewidth', 1.2)

hold on
box off

keyboard