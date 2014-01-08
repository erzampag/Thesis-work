
% DSparam_runthru_OLDCODE
%%

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

%% I want to make a 3D histogram with the count of how many neurons were DS in each SP/TEMP!
% figure
% h = bar3(hist_counter) %, 'FaceAlpha', 0.65);
% colorbar
% 
% for k = 1:length(h)
%     zdata = get(h(k), 'ZData');
%     set(h(k), 'CData', zdata, 'FaceColor', 'interp')
% end
% set(gcf,'renderer', 'opengl')
% 
% 
% % Actually, instead let's put each neuron in it's favorite SP/TEMP ONCE!!!
% hist_all = zeros(length(idList),2);
% 
% for i  = 1 : length(idList)
%     [spat_all, tem_all] = find(mean_rates{i} == max(max(mean_rates{i})));
%     if length(tem_all) > 1
%         fprintf('neuron with more than one max mean rate\n')
%         
%         continue % ignore all neurons that have more than one max spike rate
%     end
% 
%     hist_all(i,1) = spat_all; hist_all(i,2)=  tem_all;
% %     hist_all
% end
% 
% figure
% hist3(hist_all, [5,5])
% figure
% hist(hist_all(:,1))
% 
% 
% for i  = good_neurons_index
%     [spat_all, tem_all] = find(mean_rates{i} == max(max(mean_rates{i})));
%     if length(tem_all) > 1
%         fprintf('neuron with more than one max mean rate\n')
%         
%         continue % ignore all neurons that have more than one max spike rate
%     end
% 
%     hist_all(i,1) = spat_all; hist_all(i,2)=  tem_all;
% %     hist_all
% end

%%
% attempt at surface plot

% % % % %     figure
% % % % %     %     subplot(1,3,1)
% % % % %     [X,Y] = meshgrid(log(spatial), log(temporal));
% % % % %     surf(Y, X, max_rates{index});
% % % % %     colormap(gray); colorbar
% % % % %     ylabel('Spatial periods'); xlabel('Temporal Periods'); zlabel('Spike Rate');

%%
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

