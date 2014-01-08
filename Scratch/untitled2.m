%% i-j indexing
counter = 1;
for i = [1,2,3,4]
    for j = [1,2,3,4]
        matrix(i,j) = counter;
         counter = counter+1;

    end
end




%% Matlab example of mesh plot
figure
[X,Y] = meshgrid(-8:.5:8);
R = sqrt(X.^2 +Y.^2) + eps;
Z = sin(R)./R;

mesh(X,Y,real(Z));

surf(spatial, temporal, spike_rate_all_orientation{q});


%%
[neuronIDs, classes, spiketimes, EIx, EIy] = import_neuron_info('/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
    '/Users/erinzampaglione/Documents/processed_data/', '2011-02-23-0', 'data000');

%% FFT of turning curve

true_spike_rate = [0.6600; 1.9200; 3.0800; 2.5600; 0.3200; 0.6000; 0.5000; 0.9600; 0.5200; 2.5354; 2.0200; 2.7400; 1.0400; 0.5200; 0.6200; 0.4600];
% true_spike_rate = [0.5200; 2.5354; 2.0200; 2.7400; 1.0400; 0.5200; 0.6200; 0.4600; 0.6600; 1.9200; 3.0800; 2.5600; 0.3200; 0.6000; 0.5000; 0.9600];

    true_spike_rate = [7.0600
    8.6400
    8.6200
    8.1200
    4.0400
    1.2800
    1.0600
    1.0000
    0.7600
    0.7326
    0.4600
    0.4200
    0.5000
    0.5800
    1.3000
    3.8000]; % DS cell
true_spike_rate = [1.64000000000000;1.46000000000000;2.86000000000000;2.46000000000000;2.78000000000000;2.02000000000000;2.30000000000000;1.78000000000000;...
    1.68000000000000;1.78233437018518;2.74000000000000;3.12000000000000;4.40000000000000;2.78000000000000;2.48000000000000;1.78000000000000;];


Y=fft(true_spike_rate);  %% all these index_of_mins were thebesttimulus
n=length(Y);
% if n>1;
%     Y(1)=[];
% end
power = abs(Y(1:floor(n/2))).^2
% % % % power = abs(Y(1:ceil(n/2)+1)).^2
nyquist = 1/2;
% % freq = (1:n/2)/(n/2)*nyquist*(length(segments)/10);




figure;
% plot(freq, power)
plot(power)
keyboard
%%
% % % % create random dataset
% % % x = randn(200,1);
% % % y = randn(200,1);
% % % z = randn(200,1);
% % % g = [1*ones(50,1); 2*ones(50,1); 3*ones(50,1); 4*ones(50,1); ];
% % % 
% % % % call GSCATTER and capture output argument (handles to lines)
% % % h = gscatter(x, y, g);
% % % 
% % % % for each unique group in 'g', set the ZData property appropriately
% % % gu = unique(g);
% % % for k = 1:numel(gu)
% % % set(h(k), 'ZData', z( g == gu(k) ));
% % % end
% % % view(3)
% % % 
% % % keyboard
%% celltypes
% % ONSBT = intersect(intersect(ON,small), intersect(brisk, transient));
% % ONSBS = intersect(intersect(ON,small), intersect(brisk, sustained));
% % ONSST = intersect(intersect(ON,small), intersect(sluggish, transient));
% % ONSSS = intersect(intersect(ON,small), intersect(sluggish, sustained));
% % 
% % ONMBT = intersect(intersect(ON,medium), intersect(brisk, transient));
% % ONMBS = intersect(intersect(ON,medium), intersect(brisk, sustained));
% % ONMST = intersect(intersect(ON,medium), intersect(sluggish, transient));
% % ONMSS = intersect(intersect(ON,medium), intersect(sluggish, sustained));
% % 
% % ONLBT = intersect(intersect(ON,large), intersect(brisk, transient));
% % ONLBS = intersect(intersect(ON,large), intersect(brisk, sustained));
% % ONLST = intersect(intersect(ON,large), intersect(sluggish, transient));
% % ONLSS = intersect(intersect(ON,large), intersect(sluggish, sustained));
% % 
% % OFFSBT = intersect(intersect(OFF,small), intersect(brisk, transient));
% % OFFSBS = intersect(intersect(OFF,small), intersect(brisk, sustained));
% % OFFSST = intersect(intersect(OFF,small), intersect(sluggish, transient));
% % OFFSSS = intersect(intersect(OFF,small), intersect(sluggish, sustained));
% % 
% % OFFMBT = intersect(intersect(OFF,medium), intersect(brisk, transient));
% % OFFMBS = intersect(intersect(OFF,medium), intersect(brisk, sustained));
% % OFFMST = intersect(intersect(OFF,medium), intersect(sluggish, transient));
% % OFFMSS = intersect(intersect(OFF,medium), intersect(sluggish, sustained));
% % 
% % OFFLBT = intersect(intersect(OFF,large), intersect(brisk, transient));
% % OFFLBS = intersect(intersect(OFF,large), intersect(brisk, sustained));
% % OFFLST = intersect(intersect(OFF,large), intersect(sluggish, transient));
% % OFFLSS = intersect(intersect(OFF,large), intersect(sluggish, sustained));

%% Tictoc test
% % % % tic
% % % % test1 = sum((1:100).*(1:100))
% % % % toc
% % % % test2=0;
% % % % tic
% % % % for j = 1:100
% % % %  test2 = test2 + j*j;
% % % % end
% % % % test2
% % % % toc

% for i = 1 : length(textdata)
%     
%     file = strcat('/Users/erinzampaglione/Documents/Second term 1styear/', strcat('DSCells', date), '/', strcat(datafile, '.mat'));
% load(file)
% 
% 
% end
%% 

%  percentwild = (wildindexDS./wildindex)*100;
% percentdscam = (dscamindexDS./dscamindex)*100;
% 
% figure % bar graph number of DS neurons
% ave_num_DS = [mean(wildindexDS) mean(dscamindexDS)];
% error = [std(wildindexDS)/sqrt(length(wildindexDS)) std(dscamindexDS)/sqrt(length(dscamindexDS))];
% bar(ave_num_DS);
% hold on
% errorbar(ave_num_DS, error, '.k')
% x = {'Wildtype' 'DSCAM'};
% ylabel('Number of DS cells', 'FontSize', 15);
% set(gca, 'XTickLabel', x, 'FontSize', 15);
% 
% figure % bar graph percent  DS neurons
% ave_per_DS = [mean(percentwild) mean(percentdscam)];
% error_per = [std(percentwild)/sqrt(length(percentwild)) std(percentdscam)/sqrt(length(percentdscam))];
% bar(ave_per_DS);
% hold on
% errorbar(ave_per_DS, error_per, '.k')
% x = {'Wildtype' 'DSCAM'};
% ylabel('Percent DS cells', 'FontSize', 15);
% set(gca, 'XTickLabel', x, 'FontSize', 15);
% 
% figure % bar graph percent DS neurons over Healthy cells
% ave_per_DS_Health = [mean(percentwildHealth) mean(percentdscamHealth)];
% error_per_Health = [std(percentwildHealth)/sqrt(length(percentwildHealth))...
%     std(percentdscamHealth)/sqrt(length(percentdscamHealth))];
% bar(ave_per_DS_Health);
% hold on
% errorbar(ave_per_DS_Health, error_per_Health, '.k')
% x = {'Wildtype' 'DSCAM'};
% ylabel('Percent DS cells of Healthy Cells', 'FontSize', 15);
% set(gca, 'XTickLabel', x, 'FontSize', 15);
% 
% 
% 
% figure % historgram number of DS neurons
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
% 
% figure % historgram percentage of DS neurons
% hist(percentdscam, (1:2:15));
% % hist(percentdscam, (0.5:1:14.5));
% hold on
% h = findobj(gcf,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w', 'facealpha', 0.75)
% hist(percentwild, (1:2:15));
% % hist(percentwild, (0.5:1:14.5));
% i = findobj(gcf, 'Type', 'patch');
% set(i, 'facealpha', 0.75)
%  xlim([-0 20]);
% xlabel('Percent DS', 'FontSize', 15)
% ylabel('Retinas', 'FontSize', 15)
% 
% 
% figure % historgram percentage of DS neurons for Healthy Cells
% hist(percentdscamHealth, (1:2:15));
% % hist(percentdscam, (0.5:1:14.5));
% hold on
% h = findobj(gcf,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w', 'facealpha', 0.75)
% hist(percentwildHealth, (1:2:15));
% % hist(percentwild, (0.5:1:14.5));
% i = findobj(gcf, 'Type', 'patch');
% set(i, 'facealpha', 0.75)
%  xlim([-0 20]);
% xlabel('Percent DS for Healthy Cells', 'FontSize', 15)
% ylabel('Retinas', 'FontSize', 15)
% 
% % [h, p] = ttest2(percentwild, percentdscam)
% [h, p] = kstest2(percentwild, percentdscam)
% figure
% plot(wildindex, wildindexDS, 'ko');
% hold on
% plot(dscamindex, dscamindexDS, 'rx');
% xlabel('Total neurons'); ylabel('Total DS neurons')




%%
% troubleindex = find(DSI(:) >0.5 & AllCells(:,2) <0.4);

 %%
% figure
% x_fake=[0 1 0 -1];
% y_fake=[1 0 -1 0];
% 
% h_fake=compass(x_fake,y_fake);
% hold on;
% % h=compass(x,y);
% 
% 
% for i = 1 : length(AllCells)
% %     if AllCells(i,2) >0.4 %normalized magnitude vector
%          if DSI(i) > 0.5 && ratio(i) >1 % Direction selectivity index
%               ch = compass((AllCells(i,2)*cos(AllCells(i,1))), (AllCells(i,2)*sin(AllCells(i,1))), '-k');
% %               xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
% %         set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
% % set(ch,'linewidth',4);
%         
% %          polar2(AllCells(i,1), AllCells(i,2),[0 1], '-k'); % this is in the Matlab folder - wtf polar
%         hold on
%     end
% end
% title('Direction and magnitude of DS vectors when DSI >0.5 and f2/f1 >1')
%  set(h_fake,'Visible','off')


 %% DSI vs Magnitude

% rgreat1ind = find(ratio(:) >1);
% figure
% plot(DSI(rgreat1ind), AllCells(rgreat1ind,2), 'ko')
% xlabel('DSI'); ylabel('Magnitude of DS vector')
% xlim([-0.2 1.2])

% 
 %% does it add to one?
% % check = 0;
% % for i  = 1:16;
% %     check = check +norm_true_spike_rate(i);
% %     
% % end
% % check
 %% Scatterplots of data from f2/f1 classification
% % close all
% 
% % mag4ind = find(data(:,4) >=0.4); % making an index based on a minimum magnitude
% 
% 
% % let's make a new datalist for the magnitude cutoff
% 
% dataMag = zeros(length(idList), 7); %[ID xOffDS yOffDS magnitude angle ratio IDcheck]
% for i = 1 : length(idList)
%     for j = 1 : length(data)
%         if data(j,1) == idList(i)
%             
%             dataMag(i,1) = data(j,1);
%             dataMag(i,2) = data(j,2);
%             dataMag(i,3) = data(j,3);
%             dataMag(i,4) = data(j,4);
%             dataMag(i,5) = data(j,5);
%             dataMag(i,6) = ratio(i);
%             dataMag(i,7) = idList(i);
%         end
%     end
% end
% 
% % indices for mag4 separated ratios
% 
% rless1ind = find(dataMag(:,6) <=1);
% r1to5ind = find(dataMag(:,6) >1 & dataMag(:,6) <5);
% rgreat5ind = find(dataMag(:,6) >=5);
% 
% rgreat1ind = find(dataMag(:,6) >1);
% 
% % for ratios of greater than 1
% fig4 = figure(4); clf;
% plot(dataMag(rgreat1ind,2), dataMag(rgreat1ind,3), 'ko');
% hold on;
% 
% xlim([-1 1]); ylim([-1,1]);
% axis square; set(gca,'XTick',-1:0.2:1); xlabel('X DS Vector Component', 'FontSize', 15);
% ylabel('Y DS Vector Component', 'FontSize', 15);
% title(strcat('For Ratio (f2/f1) >1, Datarun-', date), 'FontSize', 15, 'FontWeight', 'bold');
% 
% 
% % for ratios of less than 1
%     fig1 = figure(1); clf;
%     plot(dataMag(rless1ind,2), dataMag(rless1ind,3), 'ko');
%     hold on;
% 
%     xlim([-1 1]); ylim([-1,1]);
%     axis square; set(gca,'XTick',-1:0.2:1); xlabel('X DS Vector Component', 'FontSize', 15);
%     ylabel('Y DS Vector Component', 'FontSize', 15);
%     title(strcat('For Ratio (f2/f1) <=1, Datarun-', date), 'FontSize', 15, 'FontWeight', 'bold');
% 
%     % for ratios of between 1 and 5
%     fig2 = figure(2); clf;
%     plot(dataMag(r1to5ind,2), dataMag(r1to5ind,3), 'ko');
%     hold on;
% 
%     xlim([-1 1]); ylim([-1,1]);
%     axis square; set(gca,'XTick',-1:0.2:1); xlabel('X DS Vector Component', 'FontSize', 15);
%     ylabel('Y DS Vector Component', 'FontSize', 15);
%     title(strcat('For Ratio (f2/f1) >1 and <5, Datarun-', date), 'FontSize', 15, 'FontWeight', 'bold');
% 
%     % for ratios greater than 5
%     fig3 = figure(3); clf;
%     plot(dataMag(rgreat5ind,2), dataMag(rgreat5ind,3), 'ko');
%     hold on;
% 
%     xlim([-1 1]); ylim([-1,1]);
%     axis square; set(gca,'XTick',-1:0.2:1); xlabel('X DS Vector Component', 'FontSize', 15);
%     ylabel('Y DS Vector Component', 'FontSize', 15);
%     title(strcat('For Ratio (f2/f1) >=5, Datarun-', date), 'FontSize', 15, 'FontWeight', 'bold');
% 



%% Histogram of f2/f1 ratios for vision-selected DS cells
%     close all
%     figure
% 
%     hist(ratio, (0.2:0.4:48.8));
%     hold on
%      xlabel('Ratio of Freq(2) to Freq(1)'); ylabel('Neurons');
%      title(strcat('Datarun-', date), 'FontSize', 15, 'FontWeight', 'bold');
%     h = findobj(gcf,'Type','patch');
%     set(h,'FaceColor','b','EdgeColor','w')


%% Differences between TTL pulses
% % % % % % % TTLdiff = [];
% % % % % % % 
% % % % % % % for i = 2:800
% % % % % % % TTLdiff(i) = TTL(i) - TTL(i-1);
% % % % % % % 
% % % % % % % end
% % % % % % % 
% % % % % % % clf
% % % % % % % 
% % % % % % % for i = 2:700
% % % % % % %     asd(i-1) = TTL(i) - TTL(i-1);
% % % % % % % end
% % % % % % % asd = asd / 20;
% % % % % % % %plot(asd(1:20), 'ko')
% % % % % % % 
% % % % % % % nTTL = 2;
% % % % % % % for i = 1:2000000
% % % % % % % x(i) = i;
% % % % % % % if (i == TTL(nTTL))
% % % % % % %     y(i)=1;
% % % % % % %     nTTL = nTTL + 1;
% % % % % % % else y(i) = 0;
% % % % % % % end
% % % % % % % end
% % % % % % % plot (y, 'ko');

% % 10 ttl pulses for each display of gratings

% 
% plot(TTL(50:100), 'ko')
% 
% for i = 1:22000
% if i == TTL(i)
%     continue
% end
%   plot (TTL(i),1, 'ko');
%   hold on
% 
% end



%% Histograms of WT vs DSCAM DS vectors (from DS_runthru.m)
% diffcells = mean(wildindex) - mean(dscamindex); % difference # of cells in average
% diffcellsDS = mean(wildindexDS) - mean(dscamindexDS); % diff in # DS cells
% percentwild = (wildindexDS./wildindex)*100;
% percentdscam = (dscamindexDS./dscamindex)*100;
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
