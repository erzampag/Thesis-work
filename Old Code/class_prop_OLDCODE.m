%% 3D scatter plots for classes, this was converted into a function

% % % % % % figure(123);
% % % % % % figure(456);
% % % % % % for j = 1:length(unique_classes)
% % % % % %     
% % % % % %     index = find(strcmp(classes_only_cells(:), unique_classes(j)));
% % % % % %     color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
% % % % % %     
% % % % % %     if strcmp(unique_classes{j}(2), 'N')
% % % % % %         figure(123)
% % % % % %         scatter3(All_norm_DOT(index), All_norm_RL(index), All_norm_RF(index),150, color(color_index,:), 'filled', 'o');
% % % % % %         text(All_norm_DOT(index), All_norm_RL(index), All_norm_RF(index), date_only_cells(index));
% % % % % % 
% % % % % %         hold on
% % % % % %     else
% % % % % %         figure(456)
% % % % % %         scatter3(All_norm_DOT(index), All_norm_RL(index), All_norm_RF(index),150, color(color_index,:), 'filled', 'o');
% % % % % %         text(All_norm_DOT(index), All_norm_RL(index), All_norm_RF(index), date_only_cells(index));
% % % % % % 
% % % % % %         hold on
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % OFF_index = [];
% % % % % % for i = 1 : length(unique_classes)
% % % % % %     OFF_index = [OFF_index strcmp(unique_classes{i}(2), 'F')];
% % % % % % end
% % % % % % 
% % % % % % figure(123)
% % % % % % % legend(unique_classes(10:17));
% % % % % % % legend(unique_classes(8:12));
% % % % % % 
% % % % % % legend(unique_classes(logical(~OFF_index)));
% % % % % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % % % % title('Plot for ON Normalized Average Properties','FontSize', 25);
% % % % % % 
% % % % % % figure(456)
% % % % % % % legend(unique_classes(1:9));
% % % % % % % legend(unique_classes(1:7));
% % % % % % legend(unique_classes(logical(OFF_index)));
% % % % % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % % % % title('Plot for OFF Normalized Average Properties','FontSize', 25);
% % % % % % 
% % % % % % keyboard


%% Plot TC and ACF (and end of ACF) on top of each other

% figure
% for i = 1 : length(unique_classes)
%     
%     index = find(strcmp(classes_only_cells(:), unique_classes(i)));
%     color_index = find(strcmp(average_properties(:,1), unique_classes(i)));
%     
%     p = plot(All_TC(index,:)');
%     set(p, 'Color', color(color_index,:), 'LineWidth', 2);
%     
%     hold on
%     
% end
% ylabel('STA (normalized?)', 'FontSize', 20); xlabel('Time to Spike (a.u.)', 'FontSize', 20);
% % % % legend(unique_classes);
% title('All Time Courses','FontSize', 25);


% figure 
% for i = 1 : length(unique_classes)
%     
%     index = find(strcmp(classes_only_cells(:), unique_classes(i)));
%     color_index = find(strcmp(average_properties(:,1), unique_classes(i)));
%     
%     p = plot(All_ACF(index,:)');
%     set(p, 'Color', color(color_index,:), 'LineWidth', 2);
%     hold on
%     
%     
%     
% %     p = plot(All_ACF(index,101:200)');
% %     set(p, 'Color', color(color_index,:), 'LineWidth', 2);
% %     hold on
%     
% % % %     q = plot(log(All_ACF(index,101:200))');
% % % %     set(q, 'Color', color(color_index,:), 'LineWidth', 2);
% % % %     hold on
%     
% % % %     q = plot(log(All_ACF(index,:))');
% % % %     set(q, 'Color', color(color_index,:), 'LineWidth', 2);
% % % %     hold on
% 
% end
% ylabel('# of pairs (normalized?)', 'FontSize', 20); xlabel('Time Differences (a.u.)', 'FontSize', 20);
% % % % legend(unique_classes);
% title('All Autocorrelation Functions','FontSize', 25);
% 
% title('All Autocorrelation Functions (log)','FontSize', 25);

%% Old colors

% % % color = rand(24,3);
% color = [0.645551874972524,0.0938200267748656,0.0195776235533187;0.479463224948888,0.525404403859336,0.330857880214071;0.639316961040108,0.530344218392863,...
%     0.424309496833137;0.544716110526763,0.861139811393332,0.270270423432065;0.647311480293128,0.484853333552102,0.197053798095456;0.543885933999639,...
%     0.393456361215266,0.821721184961310;0.721046620579811,0.671431139674026,0.429921409383266;0.522495305777102,0.741257943454207,0.887770954256354;...
%     0.993704624120852,0.520052467390387,0.391182995461163;0.218676632399634,0.347712671277525,0.769114387388296;0.105798273250228,0.149997253831683,...
%     0.396791517013617;0.109697464523194,0.586092067231462,0.808514095887345;0.0635913709751057,0.262145317727807,0.755077099007084;0.404579995857626,...
%     0.0444540922782385,0.377395544835103;0.448372912066495,0.754933267231179,0.216018915961394;0.365816176838171,0.242785357820962,0.790407217966913;...
%     0.763504640848813,0.442402313001943,0.949303911849797;0.627896379614169,0.687796085120107,0.327565434075205;0.771980385554245,0.359228210401861,...
%     0.671264370451740;0.932853570278820,0.736340074301202,0.438644982586956;0.972740854003014,0.394707475278763,0.833500595588975;0.192028349427775,...
%     0.683415866967978,0.768854252429615;0.138874202829155,0.704047430334266,0.167253545494722;0.696266337082995,0.442305413383371,0.861980478702072;];


%% 3D scatter plots, old versions

% % % scatter3(All_DOT, All_RL, All_RF, 150, 'filled', 'o');
% % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % % % % % legend(classes_only_cells);

% % % figure
% % % for j = 1:length(unique_classes)
% % %     
% % % index = find(strcmp(classes_only_cells(:), unique_classes(j)));
% % % color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
% % % 
% % %    scatter3(All_DOT(index), All_RL(index), All_RF(index),150, color(color_index,:), 'filled', 'o');
% % %     hold on
% % % end
% % % 
% % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % legend(unique_classes);
% % % title('Plot for Average Properties','FontSize', 25);
% % % 
% % % 
% % % 
% % % 
% % % figure
% % % for j = 1:length(unique_classes)
% % %     
% % % index = find(strcmp(classes_only_cells(:), unique_classes(j)));
% % % color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
% % % 
% % %    scatter3(All_DOT(index), All_RL_norm(index), All_norm_RF(index),150, color(color_index,:), 'filled', 'o');
% % %     hold on
% % % end
% % % 
% % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % legend(unique_classes);
% % % title('Plot for Average Properties','FontSize', 25);

% % % % % % % % % % % % 

% % % % % % % % % % % % % figure
% % % % % % % % % % % % % for j = 1:length(unique_classes)
% % % % % % % % % % % % %     
% % % % % % % % % % % % % index = find(strcmp(classes_only_cells(:), unique_classes(j)));
% % % % % % % % % % % % % color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
% % % % % % % % % % % % % 
% % % % % % % % % % % % %    scatter3(All_norm_DOT(index), All_norm_RL(index), All_norm_RF(index),150, color(color_index,:), 'filled', 'o');
% % % % % % % % % % % % %       text(All_norm_DOT(index), All_norm_RL(index), All_norm_RF(index), date_only_cells(index));
% % % % % % % % % % % % % 
% % % % % % % % % % % % %     hold on
% % % % % % % % % % % % % end
% % % % % % % % % % % % % 
% % % % % % % % % % % % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % % % % % % % % % % % legend(unique_classes);
% % % % % % % % % % % % % title('Plot for Normalized Average Properties','FontSize', 25);



% % % figure
% % % for j = 1:length(unique_classes)
% % %     
% % % index = find(strcmp(classes_only_cells(:), unique_classes(j)));
% % % color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
% % % 
% % %    scatter3(All_RF(index), All_DOT(index), All_ACF_SS(index),150, color(color_index,:), 'filled', 'o');
% % %    text(All_RF(index), All_DOT(index), All_ACF_SS(index), date_only_cells(index))
% % %     hold on
% % % end
% % % 
% % % xlabel('Mean RF', 'FontSize', 20); ylabel('Mean DOT x Polarity', 'FontSize', 20); zlabel('Mean ACF steady state', 'FontSize', 20)
% % % legend(unique_classes);
% % % title('Plot for Average Properties','FontSize', 25);
% % % 
% % % figure
% % % for j = 1:length(unique_classes)
% % %     
% % % index = find(strcmp(classes_only_cells(:), unique_classes(j)));
% % % color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
% % % 
% % %    scatter3(All_norm_RF(index), All_norm_DOT(index), All_norm_ACF_SS(index),150, color(color_index,:), 'filled', 'o');
% % %    text(All_norm_RF(index), All_norm_DOT(index), All_norm_ACF_SS(index), date_only_cells(index));
% % %     hold on
% % % end
% % % 
% % % xlabel('Mean RF', 'FontSize', 20); ylabel('Mean DOT x Polarity', 'FontSize', 20); zlabel('Mean ACF steady state', 'FontSize', 20)
% % % legend(unique_classes);
% % % title('Plot for Normalized Average Properties','FontSize', 25);


%% 3D scatter plots, more old versions

% % % scatter3(All_DOT, All_RL, All_RF, 150, 'filled', 'o');
% % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % % % % % legend(classes_only_cells);

% % % % % % % % % % % % figure
% % % % % % % % % % % % for j = 1:length(unique_classes)
% % % % % % % % % % % %     
% % % % % % % % % % % % index = find(strcmp(classes_only_cells(:), unique_classes(j)));
% % % % % % % % % % % % color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
% % % % % % % % % % % % 
% % % % % % % % % % % %    scatter3(All_DOT(index), All_RL(index), All_RF(index),150, color(color_index,:), 'filled', 'o');
% % % % % % % % % % % %     hold on
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % % % % % % % % % % legend(unique_classes);
% % % % % % % % % % % % title('Plot for Average Properties','FontSize', 25);

% % % % % % % % % % % % 

% % % % % % % % % % % % % figure
% % % % % % % % % % % % % for j = 1:length(unique_classes)
% % % % % % % % % % % % %     
% % % % % % % % % % % % % index = find(strcmp(classes_only_cells(:), unique_classes(j)));
% % % % % % % % % % % % % color_index = find(strcmp(average_properties(:,1), unique_classes(j)));
% % % % % % % % % % % % % 
% % % % % % % % % % % % %    scatter3(All_norm_DOT(index), All_norm_RL(index), All_norm_RF(index),150, color(color_index,:), 'filled', 'o');
% % % % % % % % % % % % %       text(All_norm_DOT(index), All_norm_RL(index), All_norm_RF(index), date_only_cells(index));
% % % % % % % % % % % % % 
% % % % % % % % % % % % %     hold on
% % % % % % % % % % % % % end
% % % % % % % % % % % % % 
% % % % % % % % % % % % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % % % % % % % % % % % legend(unique_classes);
% % % % % % % % % % % % % title('Plot for Normalized Average Properties','FontSize', 25);

%% PCA plots only
% % % % % % % 
% % % % % % %       [ACF_coeff,ACF_score] = princomp(All_ACF);
% % % % % % %       [TC_coeff,TC_score] = princomp(All_TC);
% % % % % % %       
% % % % % % %       ON_index = [];
% % % % % % %       OFF_index = [];
% % % % % % %       
% % % % % % %       for i = 1: length(classes_only)
% % % % % % %           if strcmp(classes_only(i,2), 'N')
% % % % % % %               ON_index = [ON_index i];
% % % % % % %           else
% % % % % % %               OFF_index = [OFF_index i];
% % % % % % %           end
% % % % % % %       end
% % % % % % %       
% % % % % % %       
% % % % % % %       [ACF_ON_coeff,ACF_ON_score] = princomp(All_ACF(ON_index,:));
% % % % % % %       [TC_ON_coeff,TC_ON_score] = princomp(All_TC(ON_index,:));
% % % % % % %       
% % % % % % %       [ACF_OFF_coeff,ACF_OFF_score] = princomp(All_ACF(OFF_index,:));
% % % % % % %       [TC_OFF_coeff,TC_OFF_score] = princomp(All_TC(OFF_index,:));
% % % % % % %       
%       keyboard

% Vision Plots
% % % % % % % figure %ON and OFF
% % % % % % % subplot(2,3,1); scatter(TC_score(:,1),TC_score(:,2))
% % % % % % % % gscatter(TC_score(:,1),TC_score(:,2),classes_only)
% % % % % % % xlabel('TF1', 'FontSize', 20); ylabel('TF2', 'FontSize', 20)
% % % % % % % 
% % % % % % % subplot(2,3,2); scatter(All_norm_RF, TC_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF1', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,4); scatter(ACF_score(:,1),ACF_score(:,2));
% % % % % % % xlabel('Autocorrelation 1', 'FontSize', 20); ylabel('Autocorrelation 2', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,5); scatter(All_norm_RF,ACF_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,6); scatter(TC_score(:,1),ACF_score(:,1))
% % % % % % % xlabel('TF1', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
% % % % % % % 
% % % % % % % 
% % % % % % % %Vision Plots
% % % % % % % figure % ON classes
% % % % % % % subplot(2,3,1)
% % % % % % % % gscatter(TC_ON_score(:,1),TC_ON_score(:,2),classes_only_cells(ON_index))
% % % % % % % scatter(TC_ON_score(:,1),TC_ON_score(:,2))
% % % % % % % xlabel('TF1', 'FontSize', 20); ylabel('TF2', 'FontSize', 20)
% % % % % % % 
% % % % % % % subplot(2,3,2); scatter(All_norm_RF(ON_index), TC_ON_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF1', 'FontSize', 20);
% % % % % % % title('Vision PCA plots for ON classes')
% % % % % % % 
% % % % % % % subplot(2,3,4); scatter(ACF_ON_score(:,1),ACF_ON_score(:,2));
% % % % % % % xlabel('Autocorrelation 1', 'FontSize', 20); ylabel('Autocorrelation 2', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,5); scatter(All_norm_RF(ON_index),ACF_ON_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,6); scatter(TC_ON_score(:,1),ACF_ON_score(:,1))
% % % % % % % xlabel('TF1', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
% % % % % % % 
% % % % % % % % gname(full_description(ON_index,:))
% % % % % % % % gname(classes_only(ON_index,:))
% % % % % % % 
% % % % % % %     
% % % % % % % figure % OFF classes
% % % % % % % subplot(2,3,1)
% % % % % % % % gscatter(TC_ON_score(:,1),TC_ON_score(:,2),classes_only_cells(ON_index))
% % % % % % % scatter(TC_OFF_score(:,1),TC_OFF_score(:,2))
% % % % % % % xlabel('TF1', 'FontSize', 20); ylabel('TF2', 'FontSize', 20)
% % % % % % % 
% % % % % % % subplot(2,3,2); scatter(All_norm_RF(OFF_index), TC_OFF_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF1', 'FontSize', 20);
% % % % % % % title('Vision PCA plots for OFF classes')
% % % % % % % 
% % % % % % % subplot(2,3,4); scatter(ACF_OFF_score(:,1),ACF_OFF_score(:,2));
% % % % % % % xlabel('Autocorrelation 1', 'FontSize', 20); ylabel('Autocorrelation 2', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,5); scatter(All_norm_RF(OFF_index),ACF_OFF_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,6); scatter(TC_OFF_score(:,1),ACF_OFF_score(:,1))
% % % % % % % xlabel('TF1', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
% % % % % % % 
% % % % % % % % gname(full_description(OFF_index,:))
% % % % % % % % gname(classes_only(OFF_index,:))
% % % % % % % 
% % % % % % % 
% % % % % % % % Anishchenko Plots
% % % % % % % figure % ON classes
% % % % % % % subplot(2,3,1); scatter(All_norm_RF(ON_index), TC_ON_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF1', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,2); scatter(All_norm_RF(ON_index), TC_ON_score(:,2));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF2', 'FontSize', 20);
% % % % % % % title('Anishchenko PCA plots for ON classes')
% % % % % % % 
% % % % % % % subplot(2,3,3); scatter(All_norm_RF(ON_index), TC_ON_score(:,3));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF3', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,4); scatter(All_norm_RF(ON_index),ACF_ON_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,5); scatter(All_norm_RF(ON_index),ACF_ON_score(:,2));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 2', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,6); scatter(All_norm_RF(ON_index),ACF_ON_score(:,3));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 3', 'FontSize', 20);
% % % % % % % 
% % % % % % % % gname(full_description(ON_index,:))
% % % % % % % % gname(classes_only(ON_index,:))
% % % % % % % 
% % % % % % % 
% % % % % % % figure % OFF classes
% % % % % % % subplot(2,3,1); scatter(All_norm_RF(OFF_index), TC_OFF_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF1', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,2); scatter(All_norm_RF(OFF_index), TC_OFF_score(:,2));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF2', 'FontSize', 20);
% % % % % % % title('Anishchenko PCA plots for OFF classes')
% % % % % % % 
% % % % % % % subplot(2,3,3); scatter(All_norm_RF(OFF_index), TC_OFF_score(:,3));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('TF3', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,4); scatter(All_norm_RF(OFF_index), ACF_OFF_score(:,1));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 1', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,5); scatter(All_norm_RF(OFF_index),ACF_OFF_score(:,2));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 2', 'FontSize', 20);
% % % % % % % 
% % % % % % % subplot(2,3,6); scatter(All_norm_RF(OFF_index),ACF_OFF_score(:,3));
% % % % % % % xlabel('RF (normalized)', 'FontSize', 20); ylabel('Autocorrelation 3', 'FontSize', 20);
% % % % % % % 
% % % % % % % % gname(full_description(OFF_index,:))
% % % % % % % % gname(classes_only(OFF_index,:))
% % % % % % % 
% % % % % % % % keyboard

%% Attempt at 3D gscatter

% % % figure
% % % % call GSCATTER and capture output argument (handles to lines)
% % % h = gscatter(All_DOT, All_RF, classes_only,[],[], 20);
% % % 
% % % gname(full_description)
% % % 
% % % % for each unique group in 'g', set the ZData property appropriately
% % % gu = unique(classes_only_cells);
% % % gu=gu(end:-1:1);
% % % keyboard
% % % 
% % % for k = 1:numel(gu)
% % %     set(h(k), 'ZData', All_RL(find(strcmp(classes_only_cells(:), gu(k)))));
% % % end
% % % view(3)
% % % keyboard


%% START OF OLD CODE - Loop and Graph

% % % % % % h = zeros(24,1);
% % % % % % dummy_class = {};
% % % % % % for i = [2 4 5 10 23]  % WT and HET runs only WITH OFFLBT (for normalization)
% % % % % % 
% % % % % %     
% % % % % %     if i == 16 || i== 17 || i == 18 || i == 19
% % % % % %         average_properties = properties_by_class(textdata{i,1}, 'data001');
% % % % % %     else
% % % % % %         [celltypes, average_properties, norm_ave_properties] = properties_by_class(textdata{i,1}, 'data000');
% % % % % %        
% % % % % %     end
% % % % % %     
% % % % % %     % %     hold on
% % % % % %     %% Put them Up on the graph (average properties)
% % % % % %     for j = 1:length(average_properties)
% % % % % %         
% % % % % %         % % %         if strcmp(average_properties{j,1}(1:5),'ONSBT') ||  strcmp(average_properties{j,1}(1:6),'OFFSBT')
% % % % % %         if strcmp(average_properties{j,1}(1:2),'ON')
% % % % % %             
% % % % % %             scatter3(average_properties{j,2}, average_properties{j,3}, average_properties{j,4}, 150, color(j,:), 'filled', 'o');
% % % % % %             hold on
% % % % % %             
% % % % % %         else
% % % % % %             scatter3(average_properties{j,2}, average_properties{j,3}, average_properties{j,4}, 150, color(j,:), 'filled', 'd');
% % % % % %             hold on
% % % % % %             
% % % % % %         end
% % % % % %         
% % % % % %         if ~isempty(average_properties{j,2})
% % % % % %             dummy_class{end+1} = average_properties{j,1};
% % % % % %         end
% % % % % %     end
% % % % % % title('Plot for Average Properties','FontSize', 25);

    %% Graph (normalized average properties)
% % % % % %     for j = 1:length(norm_ave_properties)
% % % % % %         
% % % % % %         % % %         if strcmp(average_properties{j,1}(1:5),'ONSBT') ||  strcmp(average_properties{j,1}(1:6),'OFFSBT')
% % % % % %         if strcmp(norm_ave_properties{j,1}(1:2),'ON')
% % % % % %             
% % % % % %             scatter3(norm_ave_properties{j,2}, norm_ave_properties{j,3}, norm_ave_properties{j,4}, 150, color(j,:), 'filled', 'o');
% % % % % %             hold on
% % % % % %             
% % % % % %         else
% % % % % %             scatter3(norm_ave_properties{j,2}, norm_ave_properties{j,3}, norm_ave_properties{j,4}, 150, color(j,:), 'filled', 'd');
% % % % % %             hold on
% % % % % %             
% % % % % %         end
% % % % % %         
% % % % % %         if ~isempty(norm_ave_properties{j,2})
% % % % % %             dummy_class{end+1} = norm_ave_properties{j,1};
% % % % % %         end
% % % % % %     end
% % % % % %     
% % % % % %     title('Plot for Normalized Average Properties', 'FontSize', 25);


    %% Diagnostics
% % %     i
% % %     norm_ave_properties
% % %     
% % %     
    %     hold on
% % % % %     
% % % % % end
% % % % % % % % keyboard
% % % % % 
% % % % % 
% % % % % % % % %         legend([celltypes])
% % % % % 
% % % % % % % % % [h, o] = legend(celltypes);
% % % % % 
% % % % % % % % legend([celltypes]);
% % % % % 
% % % % % xlabel('Mean DOT x Polarity', 'FontSize', 20); ylabel('Mean RL', 'FontSize', 20); zlabel('Mean RF', 'FontSize', 20)
% % % % % 
% % % % % legend(dummy_class);
% % % % % 
% % % % % keyboard

%% attempts to change legend
% % % % % for k = 1: length(o)
% % % % %     obj = o(k);
% % % % %     message = sprintf('Click to turn off object handle %f', obj);
% % % % %     uiwait(msgbox(message));
% % % % %     set(obj,'visible', 'off')
% % % % %
% % % % % end

% % % % %
% % % % %         legend_markers = findobj(get(h, 'Children'));
% % % % % for k = 1: length(legend_markers)
% % % % %     obj = legend_markers(k);
% % % % %         message = sprintf('Click to turn off object handle %f', obj);
% % % % %     uiwait(msgbox(message));
% % % % %     set(obj,'visible', 'off')
% % % % % end

% % legend_markers = o;

% % %         for i = 26:25+24
% % %         set(legend_markers(i), 'FaceColor', color(i-24,:));
% % % % % %                 set(o(i), 'Color', color(i-24,:));
% % % i
% % %
% % %         end

% % % %         for i = 1:24
% % % %         set(legend_markers((i-1)*3+2), 'FaceColor', color(end+1-i,:));
% % % % % % %                 set(o(i), 'Color', color(i-24,:));
% % % %
% % % %
% % % %         end