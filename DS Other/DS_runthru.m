% This script takes a list of experimental runs (organized like
% DSnotes.xls), and either uses drifting_squarewave_shorter_stimulus2.m or
% saved data from that function (or DS_analysis_trial1.m which is all
% commented out at the bottom) to do analysis on all retinas together.

% Written by Erin Zampaglione in Winter 2012.
%%
% [num txt textdata] = xlsread('DSnotes.xls', 'Sheet1');

wildindex = [];
wildindexDS = [];
wildindexHealth = [];

dscamindex = [];
dscamindexDS = [];
dscamindexHealth = [];

hetindex = [];
hetindexDS = [];
hetindexHealth = [];

DScellsindex = 0;
DScellsindex_new = [];
DScellsindex_ONorOFF = [];


wild_ave_DS_DSI = [];
wild_ave_DS_Mag = [];

dscam_ave_DS_DSI = [];
dscam_ave_DS_Mag = [];

het_ave_DS_DSI = [];
het_ave_DS_Mag = [];

wild_allneurons_DS_DSI = [];
dscam_allneurons_DS_DSI = [];
het_allneurons_DS_DSI = [];

wild_allneurons_DS_fwhh = [];
dscam_allneurons_DS_fwhh = [];
het_allneurons_DS_fwhh = [];

wild_allneurons_DS_maxspikes = [];
dscam_allneurons_DS_maxspikes = [];
het_allneurons_DS_maxspikes = [];

wild_allneurons_DS_chisquare = [];
dscam_allneurons_DS_chisquare = [];
het_allneurons_DS_chisquare = [];

ave_DSI = [];
ave_Mag = [];

%% change to 1 if you're doing analysis on already saved data
saved_data = 1;
%%
for i = 45
    % for i = (36:45)
    
    % for i = 1 : length(textdata) %textdata should be DSnotes.xls imported
    % for i = 1 : 23 % Used in DSCAM paper
    % for i = 19; %2012-05-03
    % for i = 5  % 2011-02-23 WT data002
    
    % for i = 11  %2011-08-24 DSCAM data002
    % for i = 7 % 2011-03-25 DSCAM data003
    % for i = 3 % 2011-02-08 DSCAM data002
    % for i = 12 % 2011-08-25 DSCAM data002
    % for i = 8 % 2011-03-28 DSCAM data 003
%     for i = 9 % 2011-06-28 WT  (weird chi-square)
    % for i = 17 % 2012-04-26 WT (weird chi-square)
    % for i = 18 % 2012-04-30 WT (weird chi-square)
    % for i = [20, 21, 22, 23]; % Het runs
    % for i = 20; % really big het run
    % for i = [1 2 4 5 6 9 10 15 16 17 18 19 20 21 22 23]  % WT and HET runs only
    % for i = [1 5 17 18 19 20 21 22 23]; % Runs that have EIx0 and EIy0 data
    % for i = 1
    % for i = [4 5 9 10 15 21 23 34]; % retinas WITH CRGs and FFFs (9 has no CRG))
    % for i = [4 5 9 10 15 34]; % retinas IN POSTER WITH CRGs and FFFs (9 has no CRG))  ***CHANGE AT COLOR SPECIFICATION
    
    % for i = [1 2 4 5 9 10 15]; % Retinas from the OS mapping
    % for i = [17 21 22 23]
    
    %     % This is when you have saved data!
    
    % % % % %         file = strcat('/Users/erinzampaglione/Documents/Lab_Work/DSCells/DS_saveddata_no_OS/',...
    % % % % %             textdata{i,1}, '/', strcat(textdata{i,2}, '.mat'));
    
    file = strcat('/Users/erinzampaglione/Documents/Lab_Work/DSOSCells/',...
        textdata{i,1}, '/', strcat(textdata{i,2}, '.mat'));
    
    
    load(file)
    
    output_dir = ['/Users/erinzampaglione/Documents/Lab_Work/DSOSCells//',...
        textdata{i,1}, '/', strcat(textdata{i,2})];
    %     keyboard
    %% This is when you have to do actual analysis on raw data - spikes, etc
    
    
    %     [idList, AllCells, DSI, DSI_error, ratio, noiseratio,...
    %         DScellsindex, noise, zerofreq, fwhh, reduced_chi_square, MaxSpikes, EIx, EIy,...
    %         Fone_Fzero, Ftwo_Fzero, OS_x1, OS_y1, OS_x2, OS_y2, Max_Spikes_absolute, spike_rate_all_orientation]=...
    %         MultiBarsAnalysis(textdata{i,1},textdata{i,2}, textdata{i,6});
    
    %
    % %     [idList, AllCells, DSI, DSI_error, ratio, noiseratio, DScellsindex, noise, zerofreq, fwhh...
    % %         reduced_chi_square, MaxSpikes, EIx, EIy] =...
    % %         drifting_squarewave_shorter_stimulus2(textdata{i,1},textdata{i,2}, textdata{i,6});% date, datafile, stimfile
    % % %
    
    %%
    %     DScellsindex_new = find(ratio(:) >1 & DSI(:) > 0.5 & noiseratio(:) > 5);
    Healthy_cells = find(noiseratio(:) > 5);
    
    
    ONOFF_cells = find(noiseratio(:) >5 & ratio(:) > 1);
    
    
    %     DScellsindex_new = find(ratio(:) >1 & DSI(:) > 0.5 & MaxSpikes(:) > 100);
    %     Healthy_cells = find(MaxSpikes(:) > 100);
    
    DScellsindex_all = find(DSI(:) > 0.5 & MaxSpikes(:) > 100 & noiseratio(:) > 5);
    
    DScellsindex_ONorOFF = find(ratio(:) < 1 & DSI(:) > 0.5 & MaxSpikes(:) > 100 & noiseratio(:) > 5); % ON or OFF cells
    
    DScellsindex_new = find(ratio(:) > 1 & DSI(:) > 0.5 & MaxSpikes(:) > 100 & noiseratio(:) > 5);
    % Healthy_cells = find(MaxSpikes(:) > 100 & noiseratio(:) > 5);
    
    % DScellsindex_new = find(ratio(:) >1 & DSI(:) > 0.5 & noiseratio(:) > 5 & reduced_chi_square(:) < 50);
    %
    %         figure
    %         hist(log(ratio(DScellsindex_all)))
    %         figure
    %         hist(log(ratio(DScellsindex_ONorOFF)))
    %         figure
    %         hist(log(ratio(DScellsindex_new)))
    %         keyboard
    
    %% Angle Histogram
    % % % %         figure
    % % % %         hist(AllCells(DScellsindex_new(:),1))
    % % % % title(strcat('Hist of angles of DS Cells in Retina-', textdata{i,1}),'FontSize', 25);
    
    %% EI plotting?
    
    % % %         figure
    % % %         plot(EIx, EIy, 'og')
    % % %
    % % %         hold on;
    % % %         for j = 1: length(DScellsindex_new)
    % % %             plot(EIx(DScellsindex_new(j)), EIy(DScellsindex_new(j)), 'ok');
    % % %         end
    % % %
    % % %
    % % %         for j = 1 : length(DScellsindex_ONorOFF)
    % % %                     plot(EIx(DScellsindex_ONorOFF(j)), EIy(DScellsindex_ONorOFF(j)), 'or');
    % % %         end
    
    % % % % % % % % % % % % % % Offsets taken from Georges's plotPositionEI.m
    % % % % % % % % % % % % %
    % % % % % % % % % % % % % posElectrodes = dlmread('/Users/erinzampaglione/Documents/Lab_Work/matlabscriptsforactivationmaps/512coords.txt');
    % % % % % % % % % % % % % mapXOffset = 1042;
    % % % % % % % % % % % % % mapYOffset = 845;
    % % % % % % % % % % % % % figure
    % % % % % % % % % % % % % scatter((posElectrodes(:,1) + mapXOffset)*2,(posElectrodes(:,2) + mapYOffset)*2,'.k');
    % % % % % % % % % % % % % axis equal
    % % % % % % % % % % % % % hold on
    % % % % % % % % % % % % % plot((EIx+mapXOffset)*2, (EIy+ mapYOffset)*2, 'og')
    % % % % % % % % % % % % %         for j = 1: length(DScellsindex_new)
    % % % % % % % % % % % % %             plot((EIx(DScellsindex_new(j))+mapXOffset)*2, ...
    % % % % % % % % % % % % %               (EIy(DScellsindex_new(j))+ mapYOffset)*2, 'or', 'MarkerFaceColor', 'r')
    % % % % % % % % % % % % %         end
    % % % % % % % % % % % % %
    % % % % % % % % % % % % %         for j = 1 : length(DScellsindex_ONorOFF)
    % % % % % % % % % % % % %             plot((EIx(DScellsindex_ONorOFF(j))+mapXOffset)*2, ...
    % % % % % % % % % % % % %               (EIy(DScellsindex_ONorOFF(j))+ mapYOffset)*2, 'ob', 'MarkerFaceColor', 'b')
    % % % % % % % % % % % % %         end
    % % % % % % % % % % % % %
    % % % % % % % % % % % % %         title(strcat('EI postition for DS Cells in Retina-', textdata{i,1}),'FontSize', 25);
    
    
    %%  DSI error stuff
    %     figure
    %     subplot(3,1,1)
    %     hist(DSI_error)
    %     title(strcat('DSI errors in Retina-', textdata{i,1}));
    %         xlim([0 0.4])
    %     hold on
    %
    %     if strncmp(textdata{i,4}, 'D', 1),
    %         h = findobj(gca,'Type','patch');
    %         set(h,'FaceColor','r')
    %     end
    %
    %     x = [0.1, 0.1];
    %     y = [0, 400];
    %     plot(x,y, 'k-');
    %
    %     no_noise_DS = find(ratio(:) >1 & DSI(:) > 0.5);
    %
    %     subplot(3,1,2)
    %     hist(DSI_error(no_noise_DS))
    %         xlim([0 0.4])
    %     hold on
    %
    %
    %     if strncmp(textdata{i,4}, 'D', 1),
    %         h = findobj(gca,'Type','patch');
    %         set(h,'FaceColor','r')
    %     end
    %     x = [0.1, 0.1];
    %     y = [0, 30];
    %     plot(x,y, 'k-');
    %
    %     subplot(3,1,3)
    %
    %     hist(DSI_error(DScellsindex_new))
    %         xlim([0 0.4])
    %     hold on
    %
    %     if strncmp(textdata{i,4}, 'D', 1),
    %         h = findobj(gca,'Type','patch');
    %         set(h,'FaceColor','r')
    %     end
    %     x = [0.1, 0.1];
    %     y = [0, 30];
    %     plot(x,y, 'k-');
    %
    %
    %
    %     DSI_ratio = [];
    %     for j = 1:length(idList)
    %         DSI_ratio(j) = DSI(j)/DSI_error(j);
    %     end
    %
    
    
    %% Make Vectors of Number of DS and total Cells
     fprintf('test\n')

    if strncmp(textdata{i,4}, 'W', 1), % wildtype runs
        

        wildindex = [wildindex length(idList)]; % number of wt neurons in each run
        wildindexDS = [wildindexDS length(DScellsindex_new)]; % number of wt DS neurons
        wildindexHealth = [wildindexHealth length(Healthy_cells)]; % # of strongly responding cells
        
        % % % % % % %         wildindex
        
        if DSI(DScellsindex_new) >0
            wild_ave_DS_DSI = [wild_ave_DS_DSI mean(DSI(DScellsindex_new))]; % average DSI for each retina
            wild_allneurons_DS_DSI = [wild_allneurons_DS_DSI DSI(DScellsindex_new)]; % All DSIs from all retinas
            
            wild_allneurons_DS_fwhh = [wild_allneurons_DS_fwhh fwhh(DScellsindex_new)]; % all fwhhs from all retinas
            wild_allneurons_DS_maxspikes = [wild_allneurons_DS_maxspikes MaxSpikes(DScellsindex_new)];% all MaxSpikes from all retinas
            wild_allneurons_DS_chisquare = [wild_allneurons_DS_chisquare reduced_chi_square(DScellsindex_new)]; % all X^2 from all retinas
            
        else
            %             wild_ave_DS_DSI = [wild_ave_DS_DSI 0];
%             continue
        end
        if AllCells(DScellsindex_new, 2) > 0
            wild_ave_DS_Mag = [wild_ave_DS_Mag mean(AllCells(DScellsindex_new, 2))];%ave mag for each retina
        else
            %             wild_ave_DS_Mag = [wild_ave_DS_Mag 0];
%             continue
            
        end
        
    elseif strncmp(textdata{i,4}, 'D', 1), % dscam runs
        
        dscamindex = [dscamindex length(idList)]; % number of dscam neurons in each run
        dscamindexDS = [dscamindexDS length(DScellsindex_new)]; % number of dscam DS neurons
        
        dscamindexHealth = [dscamindexHealth length(Healthy_cells)]; % # of strongly responding cells
        
        if DSI(DScellsindex_new) >0
            dscam_ave_DS_DSI = [dscam_ave_DS_DSI mean(DSI(DScellsindex_new))];% average DSI for each retina
            dscam_allneurons_DS_DSI = [dscam_allneurons_DS_DSI DSI(DScellsindex_new)]; % All DSIs from all retinas
            
            dscam_allneurons_DS_fwhh = [dscam_allneurons_DS_fwhh fwhh(DScellsindex_new)]; % all fwhhs from all retinas
            dscam_allneurons_DS_maxspikes = [dscam_allneurons_DS_maxspikes MaxSpikes(DScellsindex_new)]; % all MaxSpikes from all retinas
            dscam_allneurons_DS_chisquare = [dscam_allneurons_DS_chisquare reduced_chi_square(DScellsindex_new)]; % all X^2 from all retinas
            
        else
            %dscam_ave_DS_DSI = [dscam_ave_DS_DSI 0]
%             continue
        end
        if AllCells(DScellsindex_new, 2) > 0
            dscam_ave_DS_Mag = [dscam_ave_DS_Mag mean(AllCells(DScellsindex_new, 2))];%ave mag for each ret
        else
            %dscam_ave_DS_Mag = [dscam_ave_DS_Mag 0];%ave mag for each retina
%             continue
        end
        
        
    elseif strncmp(textdata{i,4}, 'H', 1), % heterozygous runs
        
        hetindex = [hetindex length(idList)];
        hetindexDS = [hetindexDS length(DScellsindex_new)]; % number of dscam DS neurons
        hetindexHealth = [hetindexHealth length(Healthy_cells)]; % # of strongly responding cells
        
        
        
        if DSI(DScellsindex_new) >0
            het_allneurons_DS_DSI = [het_allneurons_DS_DSI DSI(DScellsindex_new)]; % All DSIs from all retinas
            
            het_allneurons_DS_fwhh = [het_allneurons_DS_fwhh fwhh(DScellsindex_new)]; % all fwhhs from all retinas
            het_allneurons_DS_maxspikes = [het_allneurons_DS_maxspikes MaxSpikes(DScellsindex_new)]; % all MaxSpikes from all retinas
            het_allneurons_DS_chisquare = [het_allneurons_DS_chisquare reduced_chi_square(DScellsindex_new)]; % all X^2 from all retinas
        end
    end
    


    %% reduced chi squared histograms
    %     figure
    %     if strncmp(textdata{i,4}, 'W', 1),
    %
    %         %         hist(fwhh(DScellsindex_new), (0.25:0.5:4.75));
    %         hist(reduced_chi_square(DScellsindex_new), (2.5:5:52.5))
    %         hold on
    %
    %
    %     elseif strncmp(textdata{i,4}, 'D', 1), % dscam runs
    %         %         hist(fwhh(DScellsindex_new),(0.25:0.5:4.75));
    %         hist(reduced_chi_square(DScellsindex_new), (2.5:5:52.5))
    %         hold on
    %         h = findobj(gcf,'Type','patch');
    %         set(h,'FaceColor','r')
    %
    %     end
    %     hold on
    %
    %     title(strcat('Distribution of reduced chi-squared in Retina-',...
    %         textdata{i,1}));
    %     xlabel('Reduced Chi-Square', 'FontSize', 15);
    %     ylabel('Number of DS Neurons', 'FontSize', 15);
    %
    %
    % %      fprintf('%s %f \n',textdata{i,4}, length(find(reduced_chi_square(DScellsindex_new)>30)))
    %
    %
    %% Full Widths at Half Height
    %
    %     if strncmp(textdata{i,4}, 'W', 1)
    %
    %         figure
    %         hist(fwhh(DScellsindex_new), (1.25:.5:3.75));
    %         hold on;
    %
    %     elseif strncmp(textdata{i,4}, 'D', 1),
    %
    %         figure
    %         hist(fwhh(DScellsindex_new), (1.25:.5:3.75));
    %         hold on;
    %         h = findobj(gcf,'Type','patch');
    %         set(h,'FaceColor','r')
    %         %
    %     end
    %     hold on
    %
    %     title(strcat('Distribution of full width at half height in Retina-',...
    %         textdata{i,1}));
    %     xlabel('Full Width at Half Height', 'FontSize', 15);
    %     ylabel('Number of DS Neurons', 'FontSize', 15);
    
    %% Noise Ratio vs MaxSpikes
    
    %     figure
    %     plot(
    
    %% FULL RETINA GRAPHS

    if saved_data == 1

        color = zeros(27,3);
        n = 0;
        for l = 0:0.5:1
            for j = 0:0.5:1
                for k = 0:0.5:1
                    n = n+1;
                    color(n,:) = [l, j, k];
                end
            end
        end
        
        unique_dates = [];
        k=1;
        for j = 1:23; %%% Change to match i
            unique_dates{k} = textdata{j,1};
            k=k+1;
        end
        unique_dates = unique(unique_dates);
        
        % All DS neurons for each retina (probably only use this when you have
        % saved data!
        figure
        x_fake=[0 1 0 -1];
        y_fake=[1 0 -1 0];
        
        h_fake=compass(x_fake,y_fake, '--');
        hold on;
        
        for j = 1 : length(DScellsindex_new)
            
            
            color_index = 2*find(strcmp(textdata{i,1}, unique_dates));
            
            
            %     for i = 1 : length(AllCells)
            %     if AllCells(i,2) >0.4 %normalized magnitude vector
            %         if DSI(i) > 0.5 && ratio(i) > 1 % Direction selectivity index and f2/f1
            %  if strncmp(textdata{i,4}, 'W', 1)
            if strncmp(textdata{i,4}, 'W', 1) || i == 11
                
                ch = compass((AllCells(DScellsindex_new(j),2).*cos(AllCells(DScellsindex_new(j),1))),...
                    (AllCells(DScellsindex_new(j),2).*sin(AllCells(DScellsindex_new(j),1))),'-k');
                
            elseif strncmp(textdata{i,4}, 'D', 1),
                
                ch = compass((AllCells(DScellsindex_new(j),2).*cos(AllCells(DScellsindex_new(j),1))),...
                    (AllCells(DScellsindex_new(j),2).*sin(AllCells(DScellsindex_new(j),1))), '-r');
                
            elseif strncmp(textdata{i,4}, 'H', 1), % heterozygous runs
                
                ch = compass((AllCells(DScellsindex_new(j),2).*cos(AllCells(DScellsindex_new(j),1))),...
                    (AllCells(DScellsindex_new(j),2).*sin(AllCells(DScellsindex_new(j),1))), '-g');
                
            end
            
            xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
            set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
            
            %Removing the label
            % set(findall(gcf, 'String', '30','String','60', 'String', '180') ,'String', ' ')
            
            %polar2(AllCells(i,1), AllCells(i,2),[0 1], '-k');% this is in the Matlab folder
            set(ch,'LineWidth',1.5)
            %             set(ch, 'Color', color(color_index,:)); % comment out to not color for poster
            hold on
            
            
        end
        
        % %         % PLOTTING ALL CELLS
        % %         figure
        % %         polar2(AllCells(:,1), AllCells(:,2),[0 1], 'o')
        
        
        for j = 1 : length(DScellsindex_ONorOFF)
            
            ch = compass((AllCells(DScellsindex_ONorOFF(j),2).*cos(AllCells(DScellsindex_ONorOFF(j),1))),...
                (AllCells(DScellsindex_ONorOFF(j),2).*sin(AllCells(DScellsindex_ONorOFF(j),1))), '-b'); % ON or OFF cells
            
            xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
            set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
            set(ch,'LineWidth',1.5)
            
            
        end
        
        
        ch = title(strcat('Pref Dir/Mag for DS Cells in Retina-',...
            textdata{i,1}),'FontSize', 25);
        
        %         delete(findall(gcf, 'type', 'text'));
        
        set(h_fake,'Visible','off')
        
        set(gcf, 'PaperPosition', [0,0,12,12]);
        print('-depsc', [output_dir '/' 'DS_PolarPlot'])
        %         keyboard
        
        
        %% AND print each neuron (I'm gonna stick this in MultibarAnalysis too i guess?
        orientation = [0:22.5:360-22.5];
        for q = DScellsindex_all'
            figure
            
            p = polar([orientation(:)*pi/180'; orientation(1)],...
                [spike_rate_all_orientation{q}; spike_rate_all_orientation{q}(1)]);
            
            title(['Spike Rates for Neuron ' num2str(idList(q)) ', Index ' num2str(q)],'FontSize', 25)
            
            if ratio(q) <1
                set(p, 'Color', 'r')
            end
            
            set(gcf, 'PaperPosition', [0,0,12,12]);
            print('-depsc', [output_dir '/' 'Tuning' num2str(idList(q))])
            
            close(gcf)
        end
        %     keyboard
        
    end
    
    %%
    
    %         %% For choosing paper figures
    %         if i == 5 || i == 11
    %
    %             % (magnitude distibution for one retina)
    %             %             figure
    % % % % % % % % % % % % % % % % % % %             rgreat1ind = find(ratio(:) >1 & noiseratio(:) > 5 & DSI(:) > 0 & MaxSpikes(:) > 100);
    %             rgreat1ind = find(noiseratio(:) > 5 & DSI(:) > 0);
    %
    % %                         rgreat1ind = find(ratio(:) >1 & noiseratio(:) > 5 & DSI(:) > 0);
    %
    %             %             hist(AllCells((rgreat1ind),2), (0.0125:0.025:0.9875))
    %             %             xlim([0 1])
    %             %             ylim([0 41])
    %             %             xlabel('Magnitude of DS Vector', 'FontSize', 15);
    %             %             ylabel('# of Healthy ON/OFF Neurons', 'FontSize', 15);
    %             %             title(strcat...
    %             %                 ('Distribution of the Magnitudes for all Healthy ON/OFF Cells in Retina-',...
    %             %                 textdata{i,1}));
    %             %             hold on
    %             %             x = [0.5, 0.5];
    %             %             y = [0, 50];
    %             %             plot(x,y, 'k-');
    %             figure
    %             xlim([0 1]);
    % % %             ylim([0 25]);
    %             ylim([0 100]);
    %             hist(DSI(rgreat1ind), (0.025:.05:.975));
    %
    %             xlabel('Direction Selectivity Index', 'FontSize', 25);
    % % % % % %             ylabel('Healthy ON/OFF RGCs', 'FontSize', 25);
    %                         ylabel('Healthy RGCs', 'FontSize', 25);
    %
    %             title(strcat('Distribution of DSIs for all Healthy ON/OFF Cells in Retina-', textdata{i,1}));
    %             hold on
    %             box off
    %             h = findobj(gcf,'Type','patch');
    %             set(h,'FaceColor',[0.5 0.5 0.5])
    %             set(gca, 'fontsize', 20);
    %             x = [0.5, 0.5];
    % %             y = [0, 25];
    %             y = [0,100];
    %             plot(x,y, 'k-');
    %
    %
    %         end
    
    %     end
    
    %% Make classification textfile!
    % % % % %     keyboard
    % % % % %     formatted_sets = {};
    % % % % %     counter = 1;
    % % % % %     for j = 1 : length(idList)
    % % % % %         if counter <= length(DScellsindex_new) % go through all the DS cells
    % % % % %             if j == DScellsindex_new(counter) % if the neuron ID is one of the DS neurons
    % % % % %                 formatted_sets{j,1} =  idList(j);
    % % % % %                 formatted_sets{j,2} = 'All/DS';
    % % % % %                 counter = counter+1;
    % % % % %             else % if the neuron ID is not a DS cell
    % % % % %                 formatted_sets{j,1} =  idList(j);
    % % % % %                 formatted_sets{j,2} = 'All/Other';
    % % % % %             end
    % % % % %         else % all cells after DS list is exhausted
    % % % % %             formatted_sets{j,1} =  idList(j);
    % % % % %             formatted_sets{j,2} = 'All/Other';
    % % % % %         end
    % % % % %     end
    % % % % %
    % % % % %     filename = ['/Users/erinzampaglione/Documents/processed_data/',textdata{i,1}, '/',...
    % % % % %         'data000-map-', textdata{i,2}, '/DS_classification_matlab', '.txt'];
    % % % % %
    % % % % %     fileID = fopen(filename, 'w'); %% datestr(now,30)
    % % % % %
    % % % % %     for j = 1 : length(formatted_sets)
    % % % % %         fprintf(fileID, '%d %s\n', formatted_sets{j,:});
    % % % % %     end
    % % % % %
    % % % % %     fclose(fileID);
    
    
end
keyboard
%% Distributions of Maximum spike rates for all neurons
%     figure
%     hist(wild_allneurons_DS_maxspikes, (25:50:1300))
%     figure
%     hist(dscam_allneurons_DS_maxspikes, (25:50:300))
%     figure
%     plot(wild_allneurons_DS_maxspikes, wild_allneurons_DS_chisquare, 'ok');
%     hold on
%     plot(dscam_allneurons_DS_maxspikes, dscam_allneurons_DS_chisquare, 'rx');
%     xlabel('Max Spikes'); ylabel('Reduced Chi Square')
%
% figure
% hist(wild_allneurons_DS_fwhh)
% figure
% hist(het_allneurons_DS_fwhh)
% figure
% hist(dscam_allneurons_DS_fwhh)

%% Analysis on All Retinas

mean(wildindexHealth);
std(wildindexHealth)/sqrt(length(wildindexHealth));

mean(dscamindexHealth);
std(dscamindexHealth)/sqrt(length(dscamindexHealth));

mean(hetindexHealth);
std(hetindexHealth)/sqrt(length(hetindexHealth));

% percentwild = (wildindexDS./wildindex)*100;
% percentdscam = (dscamindexDS./dscamindex)*100;

percentwildHealth = (wildindexDS./wildindexHealth)*100;
percentdscamHealth = (dscamindexDS./dscamindexHealth)*100;
percenthetHealth = (hetindexDS./hetindexHealth)*100;

percentwildhetHealth = cat(2, percentwildHealth, percenthetHealth);


totalDSwild = sum(wildindexDS);
% totalcellwild = sum(wildindex);
%
totalDSdscam = sum(dscamindexDS);
% totalcelldscam = sum(dscamindex);

totalDS_het = sum(hetindexDS);
% totalpercentwild = (totalDSwild/totalcellwild)*100;
% totalpercentdscam = (totalDSdscam/totalcelldscam)*100;

% n_percentwildHealth = percentwildHealth';
% n_percentdscamHealth = percentdscamHealth';
% percentwildHealth_false = [2.8571    1.5707    5.0467   14.2566    3.4000    1.5660    8.3521];



%% Boxplot

% healthMatrix = cat(2, percentwildHealth, percentdscamHealth);
% healthMatrix = cat(2,healthMatrix, percenthetHealth);
healthMatrix = cat(2, percentwildhetHealth, percentdscamHealth);


mean(percentwildhetHealth);
mean(percentdscamHealth);


%
G1 = {'Wildtype' 'Wildtype', 'Wildtype', 'Wildtype','Wildtype' 'Wildtype', 'Wildtype', 'Wildtype',...
    'Wildtype', 'Wildtype', 'Wildtype', 'Wildtype','Wildtype', 'Wildtype', 'Wildtype', 'Wildtype'...
    'DSCAM', 'DSCAM', 'DSCAM', 'DSCAM', 'DSCAM', 'DSCAM', 'DSCAM'};

% G1 = {'Wildtype', 'Wildtype', 'Wildtype', 'Wildtype', 'Wildtype', 'Wildtype', 'Wildtype',...
%     'DSCAM', 'DSCAM', 'DSCAM', 'DSCAM', 'DSCAM', 'DSCAM', 'DSCAM'};
% G1 = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2];
% G2 = [9:15];

figure
bp = boxplot(healthMatrix, G1, 'whisker', 6, 'widths', 0.4, 'Color', 'k');
set(bp, 'linewidth', 1.2)

hold on
box off
% plot([0.96:0.01:1.03], percentwildHealth, 'ok'); % old data
plot((0.88:0.02:1.10), percentwildHealth, 'ok', 'MarkerSize', 8); % with new data
plot((1.94:0.02:2.06),percentdscamHealth, 'ok', 'MarkerSize', 8);
plot((0.98:0.02:1.04),percenthetHealth, 'dk', 'MarkerSize', 8 , 'MarkerFaceColor', [0.5 0.5 0.5]);

% plot((0.94:0.01:1.05), percentwildHealth, 'ok', 'MarkerSize', 8); % with new data
% plot((1.97:0.01:2.03),percentdscamHealth, 'ok', 'MarkerSize', 8);

% ylabel('Percent DS for Healthy Cells', 'FontSize', 20)
ylabel('ON-OFF DSRGCs/Total RGCs (%)', 'FontSize', 20);

set(findobj(gca,'Type','text'),'FontSize',10)
set(gca, 'ytick', (0:2:16), 'FontSize', 20);


plot([0.8,1.2], [mean(percentwildhetHealth), mean(percentwildhetHealth)], '--k', 'LineWidth', 2);
plot([1.8,2.2], [mean(percentdscamHealth), mean(percentdscamHealth)], '--k', 'LineWidth', 2);



% keyboard

%% historgram percentage of DS neurons for Healthy Cells

figure
hist(percentdscamHealth, (0.75:1.5:15));
% hist(percentdscamHealth, (2:4:14));
% hist(percentdscamHealth, (1:2:15));

hold on
box off
h = findobj(gcf,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w', 'facealpha', 0.75)

hist(percentwildHealth, (0.75:1.5:15));
% hist(percentwildHealth, (2:4:14));
% hist(percentwildHealth, (1:2:15));

i = findobj(gcf, 'Type', 'patch');
set(i, 'facealpha', 0.75)
xlim([-0 16]);
xlabel('Percent DS for Healthy Cells', 'FontSize', 25)
ylabel('Retinas', 'FontSize', 25)
% set(gca, 'XTickLabel', x, 'FontSize', 20);

%% bar graph Average DSI (now without averaging over retinas)
figure

wildhet_allneurons_DS_DSI = cat(2, wild_allneurons_DS_DSI, het_allneurons_DS_DSI);

std(wildhet_allneurons_DS_DSI)/sqrt(length(wildhet_allneurons_DS_DSI));
% ave_DSI = [mean(wild_ave_DS_DSI) mean(dscam_ave_DS_DSI)];
% error_DSI = [std(wild_ave_DS_DSI)/sqrt(length(wild_ave_DS_DSI))...
%     std(dscam_ave_DS_DSI)/sqrt(length(dscam_ave_DS_DSI))];

% ave_DSI = [mean(wild_allneurons_DS_DSI) mean(dscam_allneurons_DS_DSI) mean(het_allneurons_DS_DSI)];
% error_DSI = [std(wild_allneurons_DS_DSI)/sqrt(length(wild_allneurons_DS_DSI))...
%     std(dscam_allneurons_DS_DSI)/sqrt(length(dscam_allneurons_DS_DSI))...
%     std(het_allneurons_DS_DSI)/sqrt(length(het_allneurons_DS_DSI))];

ave_DSI = [mean(wildhet_allneurons_DS_DSI) mean(dscam_allneurons_DS_DSI)];
error_DSI = [std(wildhet_allneurons_DS_DSI)/sqrt(length(wildhet_allneurons_DS_DSI))...
    std(dscam_allneurons_DS_DSI)/sqrt(length(dscam_allneurons_DS_DSI))];


%  ave_DSI = [mean(wild_allneurons_DS_DSI) mean(dscam_allneurons_DS_DSI)];
% error_DSI = [std(wild_allneurons_DS_DSI)/sqrt(length(wild_allneurons_DS_DSI))...
%     std(dscam_allneurons_DS_DSI)/sqrt(length(dscam_allneurons_DS_DSI))];


bar(ave_DSI);
hold on
box off
h = findobj(gcf,'Type','patch');
set(h,'FaceColor',[0.5 0.5 0.5])
errorbar(ave_DSI, error_DSI, '.k')
ylim([0 1]);
x = {'Wildtype DS' 'DSCAM DS^-/-' 'Het DS'};

ylabel('average DSI', 'FontSize', 20);
set(gca, 'XTickLabel', x, 'FontSize', 20);

%% bar graph average fwhh
figure

wildhet_allneurons_DS_fwhh = cat(2, wild_allneurons_DS_fwhh, het_allneurons_DS_fwhh);

% USE THIS FOR THREE BARS IN A ROW
ave_fwhh = [mean(wild_allneurons_DS_fwhh) mean(dscam_allneurons_DS_fwhh) mean(het_allneurons_DS_fwhh)];
error_fwhh = [std(wild_allneurons_DS_fwhh)/sqrt(length(wild_allneurons_DS_fwhh))...
    std(dscam_allneurons_DS_fwhh)/sqrt(length(dscam_allneurons_DS_fwhh))...
    std(het_allneurons_DS_fwhh)/sqrt(length(het_allneurons_DS_fwhh))];



% error_fwhh = [std(wild_allneurons_DS_fwhh)...
%     std(dscam_allneurons_DS_fwhh)];


ave_fwhh = [mean(wildhet_allneurons_DS_fwhh) mean(dscam_allneurons_DS_fwhh)];
error_fwhh = [std(wildhet_allneurons_DS_fwhh)/sqrt(length(wildhet_allneurons_DS_fwhh))...
    std(dscam_allneurons_DS_fwhh)/sqrt(length(dscam_allneurons_DS_fwhh))];

% error_fwhh = [std(wildhet_allneurons_DS_fwhh)...
%     std(dscam_allneurons_DS_fwhh)];


bar(ave_fwhh);
hold on
box off
h = findobj(gcf,'Type','patch');
set(h,'FaceColor',[0.5 0.5 0.5])
errorbar(ave_fwhh, error_fwhh, '.k')

x = {'Wildtype DS' 'DSCAM DS' 'Het DS'};

ylabel('average FWHM (radians)', 'FontSize', 20);
set(gca, 'XTickLabel', x, 'FontSize', 20);


%% Statistics
% wild vs dscam
[hks_h, hks_p] = kstest2(percentwildHealth, percentdscamHealth);
[ht_h, ht_p] = ttest2(percentwildHealth, percentdscamHealth);
[tdsi_h, tdsi_p] = ttest2(wild_allneurons_DS_DSI, dscam_allneurons_DS_DSI);
[tfwhh_h, tfwhh_p] = ttest2(wild_allneurons_DS_fwhh, dscam_allneurons_DS_fwhh);

% het vs dscam
[iks_h, iks_p] = kstest2(percenthetHealth, percentdscamHealth);
[it_h, it_p] = ttest2(percenthetHealth, percentdscamHealth);
[udsi_h, udsi_p] = ttest2(het_allneurons_DS_DSI, dscam_allneurons_DS_DSI);
[ufwhh_h, ufwhh_p] = ttest2(het_allneurons_DS_fwhh, dscam_allneurons_DS_fwhh);

% wild vs het
[jks_h, jks_p] = kstest2(percentwildHealth, percenthetHealth);
[jt_h, jt_p] = ttest2(percentwildHealth, percenthetHealth);
[vdsi_h, vdsi_p] = ttest2(wild_allneurons_DS_DSI, het_allneurons_DS_DSI);
[vfwhh_h, vfwhh_p] = ttest2(wild_allneurons_DS_fwhh, het_allneurons_DS_fwhh);

% wild including hets vs DSCAM
[kks_h, kks_p] = kstest2(percentwildhetHealth, percentdscamHealth);
[kt_h, kt_p] = ttest2(percentwildhetHealth, percentdscamHealth);
[wdsi_h, wdsi_p] = ttest2(wildhet_allneurons_DS_DSI, dscam_allneurons_DS_DSI);
[wfwhh_h, wfwhh_p] = ttest2(wildhet_allneurons_DS_fwhh, dscam_allneurons_DS_fwhh);


fprintf('t-test p value for  percent DS cells (DSCAM vs WT): %d \n', ht_p);
fprintf('ks-test p value for  percent DS cells (DSCAM vs WT): %d \n', hks_p);
fprintf('t-test p value for DSI of DS (DSCAM vs WT): %d \n', tdsi_p);
fprintf('t-test p value for fwhh of DS (DCAM vs WT): %d \n', tfwhh_p);

fprintf('\nt-test p value for  percent DS cells (DSCAM vs Het): %d \n', it_p);
fprintf('ks-test p value for  percent DS cells (DSCAM vs Het): %d \n', iks_p);
fprintf('t-test p value for DSI of DS (DSCAM vs Het): %d \n', udsi_p);
fprintf('t-test p value for fwhh of DS (DCAM vs Het): %d \n', ufwhh_p);

fprintf('\nt-test p value for  percent DS cells (Het vs WT): %d \n', jt_p);
fprintf('ks-test p value for  percent DS cells (Het vs WT): %d \n', jks_p);
fprintf('t-test p value for DSI of DS (Het vs WT): %d \n', vdsi_p);
fprintf('t-test p value for fwhh of DS (Het vs WT): %d \n', vfwhh_p);

fprintf('\nt-test p value for  percent DS cells (DSCAM vs Het and WT): %d \n', kt_p);
fprintf('ks-test p value for percentDS cells (DSCAM vs Het and WT): %d \n', kks_p);
fprintf('t-test p value for DSI of DS (DSCAM vs Het and WT): %d \n', wdsi_p);
fprintf('t-test p value for fwhh of DS (DCAM vs Het and WT): %d \n', wfwhh_p);

keyboard

save(strcat('/Users/erinzampaglione/Documents/Lab_Work/DSCells/',...
    strcat('runthroughvariables_', datestr(now,30))))

fprintf('done \n')