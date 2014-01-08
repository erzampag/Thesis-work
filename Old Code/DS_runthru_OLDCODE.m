
%% bar graph percent DS neurons over Healthy cells

% figure
% ave_per_DS_Health = [mean(percentwildHealth) mean(percentdscamHealth)];
% 
% error_per_Health = [std(percentwildHealth)/sqrt(length(percentwildHealth))...
%     std(percentdscamHealth)/sqrt(length(percentdscamHealth))];
% 
% % error_per_Health = [std(percentwildHealth)...
% %     std(percentdscamHealth)];
% 
% bar(ave_per_DS_Health);
% hold on
% box off
% errorbar(ave_per_DS_Health, error_per_Health, '.k')
% x = {'Wildtype' 'DSCAM'};
% % ylabel('Percent DS cells of Healthy Cells', 'FontSize', 20);
% ylabel('Percent of cells identified as DS', 'FontSize', 20);
% 
% set(gca, 'XTickLabel', x, 'FontSize', 20);

%% number of DS vs healthy cells scatterplot
% figure
% plot(wildindexHealth, wildindexDS, 'ko');
% hold on
% plot(dscamindexHealth, dscamindexDS, 'rx');
% xlabel('Total Healthy neurons'); ylabel('Total DS neurons')

%% Old or unused graphs


%             DScellsindex_new = find(ratio(:) >1 & DSI(:) > 0.5 & DSI_ratio(:) > 5);
%         Healthy_cells = find(DSI_ratio(:) > 5);

%     
% %             DScellsindex_new = find(ratio(:) >1 & DSI(:) > 0.5 & noiseratio(:) > 10);
% %         Healthy_cells = find(noiseratio(:) > 5);
%     
%     %     DScellsindex_new = find(ratio(:) >1 & AllCells(:,2) > 0.5); %


% (magnitude distibution for one retina)
%             figure
%             rgreat1ind = find(ratio(:) >1 & noiseratio(:) > 5 & DSI(:) > 0);
%             hist(AllCells((rgreat1ind),2), (0.0125:0.025:0.9875))
%             xlim([0 1])
%             ylim([0 41])
%             xlabel('Magnitude of DS Vector', 'FontSize', 15);
%             ylabel('# of Healthy ON/OFF Neurons', 'FontSize', 15);
%             title(strcat...
%                 ('Distribution of the Magnitudes for all Healthy ON/OFF Cells in Retina-',...
%                 textdata{i,1}));
%             hold on
%             x = [0.5, 0.5];
%             y = [0, 50];
%             plot(x,y, 'k-');

% % historgram number of DS neurons

% figure
% hist(dscamindexDS, (2.5:5:20));
% hold on
% h = findobj(gcf,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w', 'facealpha', 0.75)
% hist(wildindexDS, (2.5:5:77.5));
% i = findobj(gcf, 'Type', 'patch');
% set(i, 'facealpha', 0.75)
% %  xlim([-0 20]);
% xlabel('Number DS', 'FontSize', 15)
% ylabel('Retinas', 'FontSize', 15)


% % historgram percentage of DS neurons

% figure
% hist(percentdscam, (1:2:15));
% % hist(percentdscam, (0.5:1:14.5));
% hold on
% h = findobj(gcf,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w', 'facealpha', 0.75)
% hist(percentwild, (1:2:15));
% % hist(percentwild, (0.5:1:14.5));
% i = findobj(gcf, 'Type', 'patch');
% set(i, 'facealpha', 0.75)
% xlim([-0 20]);
% xlabel('Percent DS', 'FontSize', 15)
% ylabel('Retinas', 'FontSize', 15)


% % Bar graph total DS cells over total cells

% figure
% lumped_percents = [totalpercentwild totalpercentdscam];
% lumped_error_per = [sqrt(totalpercentwild) sqrt(totalpercentdscam)];
% bar(lumped_percents);
% hold on
% errorbar(lumped_percents, lumped_error_per, '.k')
% x = {'Wildtype' 'DSCAM'};
% ylabel('Percent DS cells for all Retinas Lumped', 'FontSize', 15);
% set(gca, 'XTickLabel', x, 'FontSize', 15);


% % bar graph number of DS neurons

% figure
% ave_num_DS = [mean(wildindexDS) mean(dscamindexDS)];
% error = [std(wildindexDS)/sqrt(length(wildindexDS)) std(dscamindexDS)/sqrt(length(dscamindexDS))];
% bar(ave_num_DS);
% hold on
% errorbar(ave_num_DS, error, '.k')
% x = {'Wildtype' 'DSCAM'};
% ylabel('Number of DS cells', 'FontSize', 15);
% set(gca, 'XTickLabel', x, 'FontSize', 15);

% % bar graph percent  DS neurons

% figure
% ave_per_DS = [mean(percentwild) mean(percentdscam)];
% error_per = [std(percentwild)/sqrt(length(percentwild))...
%     std(percentdscam)/sqrt(length(percentdscam))];
% bar(ave_per_DS);
% hold on
% errorbar(ave_per_DS, error_per, '.k')
% x = {'Wildtype' 'DSCAM'};
% ylabel('Percent DS cells', 'FontSize', 15);
% set(gca, 'XTickLabel', x, 'FontSize', 15);

% figure % hbar?

% hbar(ave_per_DS);
% hold on
% errorbar(ave_per_DS, error_per, '.k')
% y = {'Wildtype' 'DSCAM'};
% xlabel('Percent DS cells', 'FontSize', 15);
% set(gca, 'XTickLabel', x, 'FontSize', 15);

% % bar graph Average Mag

% figure
% ave_Mag = [mean(wild_ave_DS_Mag) mean(dscam_ave_DS_Mag)];
% error_Mag = [std(wild_ave_DS_Mag)/sqrt(length(wild_ave_DS_Mag))...
%     std(dscam_ave_DS_Mag)/sqrt(length(dscam_ave_DS_Mag))];
% bar(ave_Mag);
% hold on
% errorbar(ave_Mag, error_Mag, '.k')
% x = {'Wildtype' 'DSCAM'};
% ylabel('Average DS vector Magnitude for ON/OFF DS cells', 'FontSize', 15);
% set(gca, 'XTickLabel', x, 'FontSize', 15);

% % number of DS vs number cells
% figure
% plot(wildindex, wildindexDS, 'ko');
% hold on
% plot(dscamindex, dscamindexDS, 'rx');
% xlabel('Total neurons'); ylabel('Total DS neurons')

% STATISTICS
% [h, p] = ttest2(totalpercentwild, totalpercentdscam)
% [ks_h, ks_p] = kstest2(percentwild, percentdscam);
% [t_h, t_p] = ttest2(percentwild, percentdscam);
% [tmag_h, tmag_p] = ttest2(wild_ave_DS_Mag, dscam_ave_DS_Mag);
% [tdsi_h, tdsi_p] = ttest2(wild_ave_DS_DSI, dscam_ave_DS_DSI);


% fprintf('t-test p value for DS cells over all cells minus duplicates only vs WT: %d \n', t_p);
% fprintf('ks-test p value for DS cells over all cells minus duplicates only vs WT: %d \n', ks_p);
% fprintf('t-test p value for Magnitude of DS cells vs WT: %d \n', tmag_p);


%% This is for DS_analysis_trial_1

% for i = 1 : length(textdata) %textdata should be DSnotes.xls imported
%
%
%     [data, classes, ind, ind4, date, datafile] =...
%         DS_analysis_trial1(textdata{i,1},textdata{i,2}, textdata{i,6}); % date, datafile,
%
%
%     if strncmp(textdata{i,4}, 'W', 1), % wildtype runs
%
%         wildindex = [wildindex length(ind)]; % number of wt neurons in each run
%         wildindexDS = [wildindexDS length(ind4)]; % number of wt DS neurons
%
%     elseif strncmp(textdata{i,4}, 'D', 1), % dscam runs
%
%         dscamindex = [dscamindex length(ind)]; % number of dscam neurons in each run
%         dscamindexDS = [dscamindexDS length(ind4)]; % number of dscam DS neurons
%
%     end
%
% end
%
% % histogram of 20 bins with fraction
%
% diffcells = mean(wildindex) - mean(dscamindex); % difference # of cells in average
% diffcellsDS = mean(wildindex4) - mean(dscamindex4); % diff in # DS cells
% percentwild = wildindex4/wildindex)*100;
% percentdscam = (dscamindex4/dscamindex)*100;
%
%
% fig5 = figure(5); clf;
% subplot(2,1,1);
% %hist(mag, (0.025:0.025:1));
% hist(percentwild, 'b', percentdscam, 'r');
% xlabel('# of Wildtype Neurons After Duplication Removal')
% ylabel('# of # of Neurons')
%
% subplot(2,1,2);
% hist(wildindexDS);
% xlabel('# of Wildtype DS Neurons After Duplication Removal')
% ylabel('# of # of DS Neurons')
% hold on;

%
% % Historgrams of WT neurons
% fig5 = figure(5); clf;
% subplot(2,1,1);
% %hist(mag, (0.025:0.025:1));
% hist(wildindex);
% xlabel('# of Wildtype Neurons After Duplication Removal')
% ylabel('# of # of Neurons')
%
% subplot(2,1,2);
% hist(wildindexDS);
% xlabel('# of Wildtype DS Neurons After Duplication Removal')
% ylabel('# of # of DS Neurons')
% hold on;
%
% % Histograms of DSCAM neurons
% fig6 = figure(6); clf;
% subplot(2,1,1);
% hist(dscamindex);
% xlabel('# of DSCAM Neurons After Duplication Removal')
% ylabel('# of # of Neurons')
%
% subplot(2,1,2);
% hist(dscamindexDS);
% xlabel('# of DSCAM DS Neurons After Duplication Removal')
% ylabel('# of # of DS Neurons')
% hold on;
%
% % Bargraphs of WT neurons
% fig7 = figure(7); clf;
% subplot(2,1,1);
% bar(wildindex);
% xlabel('WT Retinas'), ylabel('# of Neurons')
%
% subplot(2,1,2);
% bar(wildindexDS);
% xlabel('WT Retinas'), ylabel('# of DS Neurons')
% hold on;
%
% % Bargraphs of DSCAM neurons
% fig8 = figure(8); clf;
% subplot(2,1,1);
% %hist(mag, (0.025:0.025:1));
% bar(dscamindex);
% xlabel('DSCAM Retinas'), ylabel('# of Neurons')
%
% subplot(2,1,2);
% bar(dscamindexDS);
% xlabel('DSCAM Retinas'), ylabel('# of DS Neurons')
% hold on;
%
% fig9 = figure(9); clf;
% % subplot(2,1,1);
% %hist(mag, (0.025:0.025:1));
% hist(percentdscam, (1:2:20));
% hold on
% h = findobj(fig9,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w')
% hist(percentwild, (1:2:20));
%  xlim([-0 20]);
% % hist(percentwild, 'b', percentdscam, 'r');
% xlabel('Percent DS', 'FontSize', 15)
% ylabel('Retinas', 'FontSize', 15)