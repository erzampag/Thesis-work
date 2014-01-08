% This script takes a list of experimental runs (organized like
% DSnotes.xls), and runs properties_by_class.m, making several arrays for
% all classes, all classes plus which retina they were found in, and all
% DOT, RL, RF, TC, ACF, and ACF_SS.  It then graphs them in a variety of
% manners.
% 
% ENZ, Fall 2012


% import textdata

color = zeros(27,3);
n = 0;
for i = 0:0.5:1
    for j = 0:0.5:1
        for k = 0:0.5:1
            n = n+1;
        color(n,:) = [i, j, k];  
        end
    end
end

color2 = zeros(27,3);
n = 0;
for i = 0:0.5:1
    for j = 0:0.5:1
        for k = 0:0.5:1
            n = n+1;
        color2(n,:) = [k, i, j];  
        end
    end
end

color3 = ones(27,3);
color3([5 8 9 15 16 17 21 24],:) = [0 0 1; 0 1 0; 1 0 0; 0 0 0; 0 1 0; 0 0 1; 1 0 0; 1 0 1];

color = zeros(27,3);
n = 0;
for i = [0,0.5,1]
    for j = [0,0.5,0.75]
        for k = [0,0.5,1]
            n = n+1;
        color(n,:) = [i, j, k];  
        end
    end
end


figure
hold on
for i = 1: length(color)
    scatter(i,i,100, color(i,:), 'filled', 'o')
end

%% accumulate data
full_description = [];
classes_only = [];

All_DOT_OLD = [];
All_RL_OLD = [];
% All_DOT = [];
% All_RL = [];
All_RF = [];
All_ACF = [];
All_TC = [];
All_ACF_SS = [];
All_RF_med = [];

% All_norm_DOT = [];
% All_norm_RL = [];
% All_norm_RF = [];
All_norm_ACF = [];
% All_norm_TC = [];
All_norm_ACF_SS = [];

All_FR = [];
All_FR_bin = [];
All_NLI = [];


% % for i = [1 2 4 5 6 9 10 15 16 17 18 19 20 21 22 23]  % WT and HET runs only
% % for i = [2 4 5 9 10 22 23]  % WT and HET runs only
% %  for i = [2 4 5 7 8 9 10 11 12 13 14 22 23]  % WT, HET, DSCAM-/-
% for i = [2 4 5 10 23 24 30] % Loop over Retinas
% % for i = 31

% for i = [1 2 4 5 10 15 17 19 21 23] % Loop over Retinas (Useable data: has OFFLBT, not nec. CRGs)    
% for i = [1 2 4 5 9 10 15 17 19 21 23 34] % Loop over Retinas (Useable data: has OFFLBT, not nec. CRGs)

% for i = [4 5 9 10 15 21 23 34]; % retinas WITH CRGs and FFFs (9 has no CRG))
for i = [4 5 9 10 15 34]; % retinas IN POSTER WITH CRGs and FFFs (9 has no CRG))

    %  i
    if i == 16 || i== 17 || i == 18 || i == 19
        [data, indices, celltypes, average_properties, norm_ave_properties] = properties_by_class(textdata{i,1}, 'data001');
    else
        [data, indices, celltypes, average_properties, norm_ave_properties] = properties_by_class(textdata{i,1}, 'data000');
    end
    
    for j = 1:length(norm_ave_properties) % Loop over Celltypes
        if ~isempty(norm_ave_properties{j,2}) % if this class exists

            if isempty(full_description) % for first iteration
                full_description = char(strcat(textdata{i,1}, '-', average_properties{j,1}));
                classes_only = char(average_properties{j,1});
                date_only = char(textdata{i,1});
            else
                full_description = char(full_description, strcat(textdata{i,1}, '-', average_properties{j,1}));
                classes_only = char(classes_only, average_properties{j,1});
                date_only = char(date_only, textdata{i,1});
            end
            
            All_DOT_OLD = [All_DOT_OLD; average_properties{j,2}]; % these are from average of the DOTs
            All_RL_OLD = [All_RL_OLD; average_properties{j,3}]; % these are from average of the RLs
            All_RF = [All_RF; average_properties{j,4}];
            All_ACF = [All_ACF; average_properties{j,5}];
            All_TC = [All_TC; average_properties{j,6}];
            
            All_ACF_SS = [All_ACF_SS; average_properties{j,7}];
            
            All_RF_med = [All_RF_med; average_properties{j,11}];
            
            
          if ~isempty(average_properties{j,8})
              All_FR = [All_FR; average_properties{j,8}];
              All_FR_bin = [All_FR_bin; average_properties{j,9}];
%               All_NLI = [All_NLI; average_properties{j,10}];
          end

            
%             All_norm_DOT = [All_norm_DOT; norm_ave_properties{j,2}]; % these come from properties_by_class.m
%             All_norm_RL = [All_norm_RL; norm_ave_properties{j,3}];
%             All_norm_RF = [All_norm_RF; norm_ave_properties{j,4}];
%             All_norm_ACF = [All_norm_ACF; norm_ave_properties{j,5}];
%             All_norm_TC = [All_norm_TC; norm_ave_properties{j,6}];
%             
%             All_norm_ACF_SS = [All_norm_ACF_SS; norm_ave_properties{j,7}];
        end
    end
end
% keyboard
% convert charater arrays to cell arrays
full_description_cells = cellstr(full_description);
classes_only_cells = cellstr(classes_only);
date_only_cells = cellstr(date_only);
unique_classes = unique(classes_only_cells);
unique_dates = unique(date_only_cells);

% fitted y values for TC
All_TC_fit = [];
All_RL = zeros(length(date_only_cells), 1);
All_DOT = zeros(length(date_only_cells), 1);
All_TTP = zeros(length(date_only_cells), 1);
% [y_new(1,:), All_RL(1), All_DOT(1)] = fit_time_course((1:25)', All_TC(3,:)');

for i = 1 : length(date_only_cells)
    [All_TC_fit(i,:), All_RL(i), All_DOT(i), All_TTP(i)] = fit_time_course((1:25)', All_TC(i,:)');
% i
end
% keyboard

% New method of normalization: normalize to one retina (OFFLBT 2011-02-17)
master_retina = find(strcmp(full_description_cells(:), '2011-02-17-0-OFFLBT'));

x_axis = linspace(-400,0,length(All_TC_fit));
x_axis = x_axis([ones(length(classes_only_cells),1)],:); % repeat that vector for all time courses (y values)

x_axis_raw = linspace(-400,0,length(All_TC(1,:)));
x_axis_raw = x_axis_raw([ones(length(classes_only_cells),1)],:); % repeat that vector for all time courses (y values)

All_norm_RF = zeros(length(date_only_cells), 1);
All_norm_RL = zeros(length(date_only_cells), 1);
All_norm_DOT = zeros(length(date_only_cells), 1);
x_axis_norm = zeros(length(date_only_cells), length(All_TC_fit));

x_axis_norm_raw = zeros(length(date_only_cells), length(All_TC(1,:)));

for i = 1 : length(unique_dates)
    for j = 1 : length(date_only_cells)
        if strcmp(date_only_cells{j}, unique_dates{i})
            that_retina = find(strcmp(full_description_cells(:), strcat(unique_dates{i}, '-OFFLBT')));
            
            % % %             All_norm_RF(j) = All_RF(j).*(All_RF(master_retina) / All_RF(that_retina));
            % % %
            % % %             All_norm_RL(j) = All_RL(j).*(All_RL(master_retina) / All_RL(that_retina)); % these are now from zero of the average
            % % %             All_norm_DOT(j) = All_DOT(j).*(All_DOT(master_retina) / All_DOT(that_retina)); % these are now from DOT of the average
            % % %
            % % %             x_axis_norm(j,:) = x_axis(j,:).*(All_RL(master_retina) / All_RL(that_retina));
            % % %             x_axis_norm_raw(j,:) = x_axis_raw(j,:).*(All_RL(master_retina) / All_RL(that_retina));
            
%             All_norm_RF_true(j) = 
%             All_norm_RL_true(j) = 
%             All_norm_DOT_true(j) = 
%             
            All_norm_RF(j) = All_RF(j).*(All_RF(master_retina) / All_RF(that_retina));
            
            All_norm_RL(j) = All_RL(j).*(All_RL(master_retina) / All_RL(that_retina)); % these are now from zero of the average
            All_norm_DOT(j) = All_DOT(j).*(All_DOT(master_retina) / All_DOT(that_retina)); % these are now from DOT of the average
            
% % % % % %             % Scaling to RL
% % % % % %             x_axis_norm(j,:) = x_axis(j,:).*(All_RL(master_retina) / All_RL(that_retina));
% % % % % %             x_axis_norm_raw(j,:) = x_axis_raw(j,:).*(All_RL(master_retina) / All_RL(that_retina));
            
            % Scaling to TTP
            x_axis_norm(j,:) = x_axis(j,:).*(All_TTP(master_retina) / All_TTP(that_retina));
            x_axis_norm_raw(j,:) = x_axis_raw(j,:).*(All_TTP(master_retina) / All_TTP(that_retina));
            
        end
    end
end
%  keyboard

% Re-fit to the scaled time course
All_TC_fit_norm = [];
All_RL_norm = zeros(length(date_only_cells), 1);
All_DOT_norm = zeros(length(date_only_cells), 1);
All_TTP_norm = zeros(length(date_only_cells), 1);
%     [All_TC_fit_norm(i,:), All_RL_norm(i), All_DOT_norm(i)] = fit_time_course((x_axis_norm(21,:)'+400*25/24)*24/400, All_TC_fit(21,:)');

for i = 1 : length(date_only_cells)
%     [All_TC_fit_norm(i,:), All_RL_norm(i), All_DOT_norm(i)] = fit_time_course((x_axis_norm(i,:)'+400*25/24)*24/400, All_TC_fit(i,:)');
    
        [All_TC_fit_norm(i,:), All_RL_norm(i), All_DOT_norm(i), All_TTP_norm(i)] = fit_time_course((x_axis_norm_raw(i,:)'+400*25/24)*24/400, All_TC(i,:)');

% i
end

keyboard

% Use scaled TC to make Y-scaling as well!!

All_TC_xynorm = TC_xy_scaling(All_TC_fit, x_axis, All_RL, date_only_cells, unique_classes, ...
    unique_dates, classes_only_cells, full_description_cells, color);

% % % % % Re-fit to scaled raw data
% % % % for i = 1 : length(date_only_cells)
% % % %     fit_time_course((x_axis_norm_raw(i,:)'+400*25/24)*24/400, All_TC(i,:)');
% % % % % i
% % % % end

% new interpolated y-values for TC

All_TC_spline_norm = zeros(length(date_only_cells), 241);

% spline_x_axis = linspace(-330,0,25); % chosen from min of x_axis_norm
spline_x_axis = linspace(-330,0,241); % chosen from min of x_axis_norm
spline_x_axis = spline_x_axis([ones(length(classes_only_cells),1)],:); % repeat that vector for all time courses (y values)

for i = 1 : length(date_only_cells);
    All_TC_spline_norm(i,:) = interp1(x_axis_norm_raw(i,:), All_TC(i,:), spline_x_axis(i,:), 'spline');
end


All_norm_DOT_true= All_norm_DOT./All_DOT(master_retina);
All_RL_norm_true = All_RL_norm./-All_RL(master_retina); % from fit
All_norm_RF_true = All_norm_RF./All_RF(master_retina);

keyboard


plot_timecourses_all_retinas(x_axis_raw, All_TC, unique_classes, unique_dates, classes_only_cells, date_only_cells, color)

%% variances

for i  = 1: length(unique_classes)
    index = find(strcmp(classes_only_cells(:), unique_classes(i)));
    
    % variance
%     old_RL_var(i) = var(All_RL(index));
%     old_RF_var(i) = var(All_RF(index));
%     old_DOT_var(i) = var(All_DOT(index));
%     
%     norm_RL_var(i) = var(All_RL_norm(index));
%     norm_RF_var(i) = var(All_norm_RF(index));

% % % % %     % mean NORMALIZED RF, RL, and DOT
% % % % %     mean_RF_all(i) = mean(All_norm_RF(index));
% % % % %     mean_RL_all(i) = mean(All_RL_norm(index));
% % % % %     mean_DOT_all(i) = mean(All_DOT_norm(index));
% % % % % 
% % % % %     std_RF_all(i) = std(All_norm_RF(index));
% % % % %     std_RL_all(i) = std(All_RL_norm(index));
% % % % %     std_DOT_all(i) = std(All_DOT_norm(index));
% % % % %     
    % NOT the normalized values
    mean_RF_all(i) = mean(All_RF(index));
    mean_RL_all(i) = mean(All_RL(index));
    mean_DOT_all(i) = mean(All_DOT(index));

    std_RF_all(i) = std(All_RF(index));
    std_RL_all(i) = std(All_RL(index));
    std_DOT_all(i) = std(All_DOT(index));
    
    sem_RF_all(i) = std(All_RF(index))/sqrt(length(All_RF(index)));
    sem_RL_all(i) = std(All_RL(index))/sqrt(length(All_RL(index)));
    sem_DOT_all(i) = std(All_DOT(index))/sqrt(length(All_DOT(index)));
   

    % conformity ratio
    old_RL_cr(i) = mean(All_RL(index))/std(All_RL(index));
    old_RF_cr(i) = mean(All_RF(index))/std(All_RF(index));
    
    norm_RL_cr(i) = mean(All_RL_norm(index))/std(All_RL_norm(index));
    norm_RF_cr(i) = mean(All_norm_RF(index))/std(All_norm_RF(index));
    
    % log variance
    old_RL_var(i) = var(log(All_RL(index))); % from the fit
    old_RF_var(i) = var(log(All_RF(index)));
    old_DOT_var(i) = var(log(All_DOT(index)));
    
    norm_RL_var(i) = var(log(All_RL_norm(index))); % from the fit
%     norm_RL_var(i) = var(log(All_norm_RL(index))); % from vision
    norm_RF_var(i) = var(log(All_norm_RF(index)));
   
    % dynamic range
    old_RL_dr(i) = max(All_RL(index))/min(All_RL(index));
    old_RF_dr(i) = max(All_RF(index))/min(All_RF(index));
    
    norm_RL_dr(i) = max(All_RL_norm(index))/min(All_RL_norm(index));
    norm_RF_dr(i) = max(All_norm_RF(index))/min(All_norm_RF(index));

end

figure
scatter(old_RL_var, norm_RL_var)
hold on
plot([0,.02],[0,.02], '-')
xlabel('Variance of the log RL before norm'); ylabel('Variance of the log RL after norm');

figure
scatter(old_RF_var, norm_RF_var)
hold on
plot([0,.035],[0,.035], '-')

xlabel('Variance of the log RF before norm'); ylabel('Variance of the log RF after norm');

figure
scatter(abs(old_RL_cr(2:end)), abs(norm_RL_cr(2:end)))
hold on
% plot([-20,-6],[-20,-6], '-')
plot([6,20],[6,20], '-')

xlabel('Conformity ratio for RL before norm'); ylabel('Conformity ratio for RL after norm');

figure
scatter(old_RF_cr(2:end), norm_RF_cr(2:end))
hold on
plot([5,35],[5,35], '-')
xlabel('Conformity ratio for RF before norm'); ylabel('Conformity ratio for RF after norm');






% ttests
index_OFFLBT = find(strcmp(classes_only_cells(:), 'OFFLBT'));
index_OFFMBT = find(strcmp(classes_only_cells(:), 'OFFMBT'));
index_OFFSSS = find(strcmp(classes_only_cells(:), 'OFFSSS'));
index_ONLBT = find(strcmp(classes_only_cells(:), 'ONLBT'));
index_ONMBT = find(strcmp(classes_only_cells(:), 'ONMBT'));

% ONLBT vs ONMBT: RL insig -> sig, RF gets worse
[h p] = ttest2(All_RL(index_ONLBT), All_RL(index_ONMBT))
[h p] = ttest2(All_RL_norm(index_ONLBT), All_RL_norm(index_ONMBT))
[h p] = ttest2(All_RL_norm_true(index_ONLBT), All_RL_norm_true(index_ONMBT))

[h p] = ttest2(All_RF(index_ONLBT), All_RF(index_ONMBT))
[h p] = ttest2(All_norm_RF(index_ONLBT), All_norm_RF(index_ONMBT))

[h p] = ttest2(All_TTP(index_ONLBT), All_TTP(index_ONMBT))
[h p] = ttest2(All_TTP_norm(index_ONLBT), All_TTP_norm(index_ONMBT))

figure
% ave_RL = [mean(All_RL(index_ONLBT)),mean(All_RL(index_ONMBT)); mean(All_RL_norm_true(index_ONLBT)), mean(All_RL_norm_true(index_ONMBT))];
% ave_RL = [mean(All_RL(index_ONLBT)),mean(All_RL(index_ONMBT)); mean(All_RL_norm(index_ONLBT)), mean(All_RL_norm(index_ONMBT))];
% bar(-ave_RL, 'grouped')
% hold on
% std_RL = [std(All_RL(index_ONLBT)),std(All_RL(index_ONMBT)); std(All_RL_norm(index_ONLBT)), std(All_RL_norm(index_ONMBT))];
% errorbar(-ave_RL, std_RL, '.k')



ave_RL = [mean(All_RL(index_ONLBT)) mean(All_RL(index_ONMBT))];
std_RL = [std(All_RL(index_ONLBT)) std(All_RL(index_ONMBT))];
error_RL = std_RL./[sqrt(length(All_RL(index_ONLBT))) sqrt(length(All_RL(index_ONMBT)))];
figure
bar(-ave_RL)
hold on; box off
% errorbar(-ave_RL, std_RL, '.k')
errorbar(-ave_RL, error_RL, '.k')
x = {'ONLBT' 'ONMBT'};
set(gca, 'XTickLabel', x, 'FontSize', 20);


% ave_RL_norm = [mean(All_RL_norm_true(index_ONLBT)), mean(All_RL_norm_true(index_ONMBT))];
% std_RL_norm = [std(All_RL_norm_true(index_ONLBT)), std(All_RL_norm_true(index_ONMBT))];
ave_RL_norm = [mean(All_RL_norm(index_ONLBT)), mean(All_RL_norm(index_ONMBT))];
std_RL_norm = [std(All_RL_norm(index_ONLBT)), std(All_RL_norm(index_ONMBT))];
error_RL_norm = std_RL_norm./[sqrt(length(All_RL_norm(index_ONLBT))) sqrt(length(All_RL_norm(index_ONMBT)))];

figure
bar(-ave_RL_norm)
hold on; box off
% errorbar(-ave_RL_norm, std_RL_norm, '.k')
errorbar(-ave_RL_norm, error_RL_norm, '.k')
x = {'ONLBT' 'ONMBT'};
set(gca, 'XTickLabel', x, 'FontSize', 20);


% OFFMBT vs OFFSSS RL p-value goes down, RF gets worse
[h p] = ttest2(All_RL(index_OFFMBT), All_RL(index_OFFSSS))
[h p] = ttest2(All_RL_norm(index_OFFMBT), All_RL_norm(index_OFFSSS))

[h p] = ttest2(All_RF(index_OFFMBT), All_RF(index_OFFSSS))
[h p] = ttest2(All_norm_RF(index_OFFMBT), All_norm_RF(index_OFFSSS))

[h p] = ttest2(All_TTP(index_OFFMBT), All_TTP(index_OFFSSS))
[h p] = ttest2(All_TTP_norm(index_OFFMBT), All_TTP_norm(index_OFFSSS))

% OFFLBT vs OFFMBT, RL pvalue goes down, RF goes from insig to sig
[h p] = ttest2(All_RL(index_OFFMBT), All_RL(index_OFFLBT))
[h p] = ttest2(All_RL_norm(index_OFFMBT), All_RL_norm(index_OFFLBT))
[h p] = ttest2(All_RL_norm_true(index_OFFMBT), All_RL_norm_true(index_OFFLBT))


[h p] = ttest2(All_RF(index_OFFMBT), All_RF(index_OFFLBT))
[h p] = ttest2(All_norm_RF(index_OFFMBT), All_norm_RF(index_OFFLBT))

[h p] = ttest2(All_TTP(index_OFFMBT), All_TTP(index_OFFLBT))
[h p] = ttest2(All_TTP_norm(index_OFFMBT), All_TTP_norm(index_OFFLBT))


ave_RF = [mean(All_RF(index_OFFLBT)) mean(All_RF(index_OFFMBT))];
std_RF = [std(All_RF(index_OFFLBT)) std(All_RF(index_OFFMBT))];
error_RF = std_RF./[sqrt(length(All_RF(index_OFFLBT))) sqrt(length(All_RF(index_OFFMBT)))];
figure
bar(ave_RF)
hold on; box off
% errorbar(-ave_RL, std_RL, '.k')
errorbar(ave_RF, error_RF, '.k')
x = {'OFFLBT' 'OFFMBT'};
set(gca, 'XTickLabel', x, 'FontSize', 20);


% ave_RL_norm = [mean(All_RL_norm_true(index_ONLBT)), mean(All_RL_norm_true(index_ONMBT))];
% std_RL_norm = [std(All_RL_norm_true(index_ONLBT)), std(All_RL_norm_true(index_ONMBT))];
ave_RF_norm = [mean(All_norm_RF(index_OFFLBT)), mean(All_norm_RF(index_OFFMBT))];
std_RF_norm = [std(All_norm_RF(index_OFFLBT)), std(All_norm_RF(index_OFFMBT))];
error_RF_norm = std_RF_norm./[sqrt(length(All_norm_RF(index_ONLBT))) sqrt(length(All_norm_RF(index_ONMBT)))];

figure
bar(ave_RF_norm)
hold on; box off
% errorbar(-ave_RL_norm, std_RL_norm, '.k')
errorbar(ave_RF_norm, error_RF_norm, '.k')
x = {'OFFLBT' 'OFFMBT'};
set(gca, 'XTickLabel', x, 'FontSize', 20);

 
%% correlations

OFF_LBT_index = [];
OFF_MBT_index = [];
OFF_SSS_index = [];
ON_LBT_index = [];

for i = 1 : length(unique_dates) 
    for j = 1 : length(classes_only_cells)
        if strcmp(date_only_cells{j}, unique_dates{i})
            if strcmp(classes_only_cells{j}, 'OFFLBT')
                OFF_LBT_index(i) = j;
            end
            if strcmp(classes_only_cells{j}, 'OFFMBT')
                OFF_MBT_index(i) = j;
            end
            if strcmp(classes_only_cells{j}, 'OFFSSS')
                OFF_SSS_index(i) = j;
            end
            if strcmp(classes_only_cells{j}, 'ONLBT')
                ON_LBT_index(i) = j;
            end
        end
    end
end

removal_index = [];
for i = 1 : length(OFF_LBT_index)
        if OFF_LBT_index(i) == 0 || OFF_MBT_index(i) == 0
    %     if OFF_LBT_index(i) == 0 || OFF_SSS_index(i) == 0
%     if OFF_LBT_index(i) == 0 || ON_LBT_index(i) == 0
        
        removal_index = [removal_index i];
    end
end

OFF_LBT_index(removal_index) = []; OFF_MBT_index(removal_index) = []; OFF_SSS_index(removal_index) = []; ON_LBT_index(removal_index) = [];

% figure
% scatter(All_RL(OFF_LBT_index),All_RL(ON_LBT_index));
% text(All_RL(OFF_LBT_index),All_RL(ON_LBT_index), date_only_cells(OFF_LBT_index));
% title('Response Latencies','FontSize', 25); xlabel('OFFLBT','FontSize', 20); ylabel('ONLBT','FontSize', 20);




figure
scatter(All_RL(OFF_LBT_index),All_RL(OFF_MBT_index));
text(All_RL(OFF_LBT_index),All_RL(OFF_MBT_index), date_only_cells(OFF_LBT_index));
title('Response Latencies','FontSize', 25); xlabel('OFFLBT','FontSize', 20); ylabel('OFFMBT','FontSize', 20);
hold on
plot([-140, -90], [-140,-90], '-')
% 
figure
scatter(All_RF(OFF_LBT_index),All_RF(OFF_MBT_index));
text(All_RF(OFF_LBT_index),All_RF(OFF_MBT_index), date_only_cells(OFF_LBT_index));
title('Receptive Fields','FontSize', 25); xlabel('OFFLBT','FontSize', 20); ylabel('OFFMBT','FontSize', 20);
hold on
plot([120, 190], [120, 190], '-')

keyboard

% figure
% scatter(All_RL(OFF_LBT_index),All_RL(OFF_SSS_index));
% title('Response Latencies','FontSize', 25); xlabel('OFFLBT','FontSize', 20); ylabel('OFFSSS','FontSize', 20);
% 
% figure
% scatter(All_RF(OFF_LBT_index),All_RF(OFF_SSS_index));
% title('Receptive Fields','FontSize', 25); xlabel('OFFLBT','FontSize', 20); ylabel('OFFSSS','FontSize', 20);

%% Time Courses

figure
for i = 1 : length(unique_classes)
    
    index = find(strcmp(classes_only_cells(:), unique_classes(i)));
    subplot(4,4,i); % suplot(5,4,i) if there's more than 16 cell types
    hold on
    
    for j = 1 : length(index)
        
        color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
         
%         plot(x_axis_norm(index(j),:), All_TC(index(j),:)', 'b');
%         plot(x_axis(index(j),:), All_TC(index(j),:)'+0.4, 'r');
        
        plot(spline_x_axis(index(j),:), All_TC_spline_norm(index(j),:)' -0.5, 'k') % timecourses from normalized spline
        
% %         plot(linspace(-400,0,length(All_TC_fit_norm(index(j),:))), All_TC_fit_norm(index(j),:)'+0.5, 'k'); % time course from normalized fit
        
                plot(linspace(-400,0,length(All_TC_xynorm(index(j),:))), All_TC_xynorm(index(j),:)'+0.5, 'k'); % time course from normalized fit

        
        
        plot(linspace(-400,0,length(All_TC_fit(index(j),:))), All_TC_fit(index(j),:)', 'b') % fitted time courses (translated up)
        
        plot(x_axis_raw, All_TC(index(j),:)', 'r.'); % raw data (translated up)   

    end
%     ylim([-0.7 0.7]);
    ylim([-1.1 1.1]);
    plot((-400:-1),zeros(1,400),':');
    plot((-400:-1),zeros(1,400)+0.5,':');
    plot((-400:-1),zeros(1,400)-0.5,':');

    title(unique_classes(i));
    
end

figure
for i = 1 : length(unique_classes)
    index = find(strcmp(classes_only_cells(:), unique_classes(i)));
    
%     subplot(4,4,i); % suplot(5,4,i) if there's more than 16 cell types
    subplot(2,3,i); % suplot(5,4,i) if there's more than 16 cell types
    hold on
    for j = 1 : length(index)
        color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
        
%         plot(x_axis_raw, All_TC(index(j),:)'+ (((color_index/2)-1)/5), 'Color', color(color_index,:), 'LineWidth', 1.25); %color by date

% %         plot(linspace(-400,0,length(All_TC_fit_norm(index(j),:))), All_TC_fit_norm(index(j),:)'+ (((color_index/2)-1)/5),...
% %             'Color', color(color_index,:), 'LineWidth', 1.25); % SCALED color by date FIT
        
                plot(linspace(-400,0,length(All_TC_xynorm(index(j),:))), All_TC_xynorm(index(j),:)'+ (((color_index/2)-1)/5),...
            'Color', color(color_index,:), 'LineWidth', 1.25); % SCALED X AND Y color by date FIT
        
        
%         plot(linspace(-330,0,length(All_TC_spline_norm(index(j),:))), All_TC_spline_norm(index(j),:)'+ (((color_index/2)-1)/5),...
%             'Color', color(color_index,:), 'LineWidth', 1.25); % SCALED color by date SPLINE

        plot((-400:-1),zeros(1,400)+(((color_index/2)-1)/5),':');
        
%         text(-400,(((color_index/2)-1)/5), unique_dates{color_index/2});All_DOT        
    end
    title(unique_classes(i), 'FontSize', 20);
    ylim([-0.7 (length(unique_dates)-1)/5+0.7]);

end
 %% ACFs

figure
for i = 1 : length(unique_classes)
    
    index = find(strcmp(classes_only_cells(:), unique_classes(i)));
    %subplot(5,4,i);
%     subplot(4,4,i);
    subplot(2,3,i);
    hold on
    
% %     text(50,length(unique_dates)-1/10, 'test')
%     RF_string = strcat('Ave RF =', num2str(mean_RF_all(i)), ' +/- ', num2str(std_RF_all(i))); 
%     RL_string = strcat('Ave RL =', num2str(mean_RL_all(i)), ' +/- ', num2str(std_RL_all(i))); 
%     DOT_string = strcat('Ave DOT =', num2str(mean_DOT_all(i)), ' +/- ', num2str(std_DOT_all(i))); 

    RF_string = ['Ave RF = ', num2str(mean_RF_all(i)), ' +/- ', num2str(sem_RF_all(i))];
    RL_string = ['Ave RL = ', num2str(mean_RL_all(i)), ' +/- ', num2str(sem_RL_all(i))]; 
    DOT_string = ['Ave DOT = ', num2str(mean_DOT_all(i)), ' +/- ', num2str(sem_DOT_all(i))]; 
    
    
    text(40,0.9, RF_string, 'FontSize', 15)
    text(40,0.8, RL_string, 'FontSize', 15)
    text(40,0.7, DOT_string, 'FontSize', 15)
 
    
    
    for j = 1 : length(index)
        
        color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
        
        %plot(All_ACF(index(j),:)', 'b');
        plot(linspace(0,100, length(All_ACF(index(j),:))), All_ACF(index(j),:)' + (((color_index/2)-1)/10),...
            'Color', color(color_index,:), 'LineWidth', 1.25);
        plot((0:99),zeros(1,100)+(((color_index/2)-1)/10),':');
        
%         text(150,(((color_index/2)-1)/10+.01), unique_dates{color_index/2});
    end
    ylim([0 (length(unique_dates)-1)/10+0.6]);
    title(unique_classes(i), 'FontSize', 20);
end


% to make legend
figure
index = find(strcmp(classes_only_cells(:), unique_classes(1)));
legend_index = [];
hold on
for j = 1 : length(index)
    color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
    legend_index = [legend_index color_index/2];
    plot(All_ACF(index(j),:)' + (((color_index/2)-1)/10), 'Color', color(color_index,:), 'LineWidth', 2);
    
    title(unique_classes(1));
end
legend(unique_dates(legend_index))


%% RFs

figure
for i = 1 : length(unique_classes)
    
    index = find(strcmp(classes_only_cells(:), unique_classes(i)));
    %subplot(5,4,i);
    subplot(4,4,i);
    hold on
    for j = 1 : length(index)
        
        color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
        
        %plot(All_ACF(index(j),:)', 'b');
        rectangle('Position', [0.8, (((color_index/2)-1)/5), All_RF(index(j))/1000, All_RF(index(j))/1000],...
            'Curvature', [1, 1], 'FaceColor', color(color_index,:))
                
        rectangle('Position', [1, (((color_index/2)-1)/5), All_norm_RF(index(j))/1000, All_norm_RF(index(j))/1000],...
            'Curvature', [1, 1], 'FaceColor', color(color_index,:))

       plot((0.6:0.1:1.4),zeros(1,length(0.6:0.1:1.4))+(((color_index/2)-1)/5),':');
       
%        text(1,(((color_index/2)-1)/5), unique_dates{color_index/2});
    end
    ylim([0 (length(unique_dates)-1)/5+0.2]);
    title(unique_classes(i));
        axis equal
    %     daspect([1,1,1])
%     axis auto
end

keyboard
%% Flash Response etc
% % % for unnormalized flash responses
% if ~isempty(All_FR)
%     figure
%     for i = 1 : length(unique_classes)
%         
%         index = find(strcmp(classes_only_cells(:), unique_classes(i)));
%         %subplot(5,4,i);
%         subplot(4,4,i);
%         hold on
%         for j = 1 : length(index)
%             
%             color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
%             
%             %plot(All_ACF(index(j),:)', 'b');
%             plot(All_FR(index(j),:)' + (((color_index/2)-1)/.01), 'Color', color(color_index,:), 'LineWidth', 1.25);
%             plot((0:200),zeros(1,201)+(((color_index/2)-1)/.01),':');
%             
%             text(0,(((color_index/2)-1)/.01), unique_dates{color_index/2});
%         end
%         ylim([0 (length(unique_dates)-1)/.01+200]);
%         xlim([0 201])
%         title(unique_classes(i));
%     end    
% end

%%% normalized flash responses
if ~isempty(All_FR)
    figure
    for i = 1 : length(unique_classes)
        
        index = find(strcmp(classes_only_cells(:), unique_classes(i)));
        %subplot(5,4,i);
%         subplot(4,4,i);
        subplot(2,3,i);
        hold on
        for j = 1 : length(index)
            
            color_index = 2*find(strcmp(date_only_cells(index(j)), unique_dates));
            
            %plot(All_ACF(index(j),:)', 'b');
            plot(linspace(0,8, length(All_FR(index(j),:))), All_FR(index(j),:)' + (((color_index/2)-1)/5),...
                'Color', color(color_index,:), 'LineWidth', 1.25);
            plot(linspace(0,8,200),zeros(1,200)+(((color_index/2)-1)/5),':');
            
%             text(0,(((color_index/2)-1)/5), unique_dates{color_index/2});
            
        end
        ylim([0 (length(unique_dates)-1)/5+0.5]);
        xlim([0 8])
        title(unique_classes(i), 'FontSize', 20);
    end
    
end

keyboard
%% 3D Scatter plots

% Spike Triggered Average plot
scatter3_by_classes(101, 201, 'ON Normalized Average STA Properties', 'OFF Normalized Average STA Properties',...
    All_norm_DOT, 'Mean DOT x Polarity', All_norm_RL, 'Mean RL', All_norm_RF, 'Mean RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(101, 201, 'ON Normalized Average STA Properties', 'OFF Normalized Average STA Properties',...
    All_norm_DOT./All_DOT(master_retina), 'Mean DOT x Polarity', All_norm_RL./-All_RL(master_retina), 'Mean RL',...
    All_norm_RF./All_RF(master_retina), 'Mean RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);


% Spike Triggered Average plot NOT NORMALIZED
scatter3_by_classes(105, 205, 'ON Average STA Properties', 'OFF Average STA Properties',...
    All_DOT, 'Mean DOT x Polarity', All_RL, 'Mean RL', All_RF, 'Mean RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

% Anishchenko plot
scatter3_by_classes(102, 202, 'ON Normalized Anishchenko plots', 'OFF Normalized Anishchenko plots',...
    All_norm_RF, 'Mean RF', All_norm_DOT, 'Mean DOT x Polarity', All_ACF_SS, 'ACF steady rate',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color);

% using nonlinearity index?
% scatter3_by_classes(103, 203, 'ON Normalized Average Properties', 'OFF Normalized Average Properties',...
%     All_norm_RF, 'Mean RF', All_norm_DOT, 'Mean DOT x Polarity', All_ACF_SS, 'ACF steady rate',...
%     unique_classes, classes_only_cells, average_properties, date_only_cells, color);

% PCA?
supervector = cat(2, All_TC_spline_norm, All_ACF, All_FR, All_norm_RF);
% supervector = cat(2, All_TC_spline_norm, All_ACF, All_FR);
% supervector = cat(2, All_ACF, All_FR);
% supervector = All_TC_spline_norm;
[coeff,score] = princomp(zscore(supervector));
[coeff,score] = princomp(supervector);

figure
scatter3(score(:,1), score(:,2), score(:,3))
xlabel('PC1', 'FontSize', 20); ylabel('PC2', 'FontSize', 20); zlabel('PC3', 'FontSize', 20);
title('PCA on Concatinated Vector', 'FontSize', 25);

scatter3_by_classes(104, 204, 'Plot for ON Normalized Average Properties', 'Plot for OFF Normalized Average Properties',...
    score(:,1), 'PC1', score(:,2), 'PC2', score(:,3), 'PC3',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);


% PCA ON vs OFF
ON_index = [];
OFF_index = [];

for i = 1: length(classes_only)
    if strcmp(classes_only(i,2), 'N')
        ON_index = [ON_index i];
    else
        OFF_index = [OFF_index i];
    end
end

[coeff_ON_TC,score_ON_TC] = princomp(All_TC_spline_norm(ON_index, :));
[coeff_OFF_TC,score_OFF_TC] = princomp(All_TC_spline_norm(OFF_index, :));

[coeff_ON_ACF,score_ON_ACF] = princomp(All_ACF(ON_index, :));
[coeff_OFF_ACF,score_OFF_ACF] = princomp(All_ACF(OFF_index, :));

[coeff_ON_FR,score_ON_FR] = princomp(All_FR(ON_index, :));
[coeff_OFF_FR,score_OFF_FR] = princomp(All_FR(OFF_index, :));

% re-assemble ONs and OFFs
score_ALL_TC(ON_index,:) = score_ON_TC;
score_ALL_TC(OFF_index,:) = score_OFF_TC;

score_ALL_ACF(ON_index,:) = score_ON_ACF;
score_ALL_ACF(OFF_index,:) = score_OFF_ACF;

score_ALL_FR(ON_index,:) = score_ON_FR;
score_ALL_FR(OFF_index,:) = score_OFF_FR;

scatter3_by_classes(106, 206, 'Plot for ON TC PCs', 'Plot for OFF TC PCs',...
    score_ALL_TC(:,1), 'PC1', score_ALL_TC(:,2), 'PC2', score_ALL_TC(:,3), 'PC3',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(107, 207, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,1), 'TC PC1', score_ALL_ACF(:,1), 'ACF PC1', score_ALL_FR(:,1), 'FR PC1',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(108, 208, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,2), 'TC PC2', score_ALL_ACF(:,1), 'ACF PC1', score_ALL_FR(:,2), 'FR PC2',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(109, 209, 'Plot for ON PCs and RF', 'Plot for OFF PCs and RF',...
    score_ALL_TC(:,2), 'TC PC2', score_ALL_ACF(:,1), 'ACF PC1', All_norm_RF, 'Mean RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(110, 210, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,2), 'TC PC2', score_ALL_FR(:,1), 'FR PC1', All_norm_RF, 'Mean RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);



scatter3_by_classes(111, 211, 'ON Average STA Properties', 'OFF Average STA Properties',...
    All_DOT_OLD, 'Mean OLD DOT x Polarity', All_RL, 'Mean RL', All_RF, 'Mean RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(112, 212, 'Plot for ON PCs', 'Plot for OFF PCs',...
    All_norm_RL, 'Mean RL', score_ALL_ACF(:,1), 'ACF PC1', score_ALL_FR(:,1), 'FR PC1',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(113, 213, 'Plot for ON PCs', 'Plot for OFF PCs',...
    All_norm_RL, 'Mean RL', All_norm_RF, 'Mean NORM RF', score_ALL_FR(:,1), 'FR PC1',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(114, 214, 'Plot for ON PCs', 'Plot for OFF PCs',...
    All_norm_RL, 'Mean RL', All_norm_RF, 'Mean RF', score_ALL_FR(:,2), 'FR PC2',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(115, 215, 'Plot for ON PCs', 'Plot for OFF PCs',...
    All_n, 'Mean RL', All_norm_RF, 'Mean RF', score_ALL_FR(:,2), 'FR PC2',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);




scatter3_by_classes(116, 216, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,1), 'TC PC1', score_ALL_ACF(:,1), 'ACF PC1', All_norm_RF, 'RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(117, 217, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,1), 'TC PC1', score_ALL_TC(:,2), 'TC PC2', All_norm_RF, 'RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(118, 218, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,1), 'TC PC1', score_ALL_FR(:,1), 'FR PC1', All_norm_RF, 'RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(119, 219, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,1), 'TC PC1', score_ALL_TC(:,2), 'TC PC2', score_ALL_FR(:,1), 'FR PC1',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(120, 220, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,1), 'TC PC1', score_ALL_TC(:,2), 'TC PC2', score_ALL_FR(:,2), 'FR PC2',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(121, 221, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,1), 'TC PC1', score_ALL_FR(:,1), 'FR PC1', score_ALL_FR(:,2), 'FR PC2',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(122, 222, 'Plot for ON PCs', 'Plot for OFF PCs',...
    score_ALL_TC(:,1), 'TC PC1', score_ALL_FR(:,1), 'FR PC1', All_norm_RF, 'norm RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

% STA with TTP instead of RL
scatter3_by_classes(123, 223, 'ON Normalized Average STA Properties', 'OFF Normalized Average STA Properties',...
    All_DOT_norm, 'Mean DOT x Polarity', All_TTP_norm, 'Mean TTP', All_norm_RF, 'Mean RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);

scatter3_by_classes(124, 224, 'ON Average STA Properties', 'OFF Average STA Properties',...
    All_DOT, 'Mean DOT x Polarity', All_TTP, 'Mean TTP', All_RF, 'Mean RF',...
    unique_classes, classes_only_cells, average_properties, date_only_cells, color3);



%% Clustering by PCA

PCA_select(classes_only, date_only_cells, All_TC_spline_norm, All_ACF, All_norm_RF,...
    All_FR, All_norm_DOT, All_norm_RL, color, unique_dates) % normalize TC with spline

PCA_select(classes_only, date_only_cells, All_TC_spline_norm, All_FR, All_norm_RF,...
    All_ACF, All_norm_DOT, All_norm_RL, color, unique_dates) % normalize TC with spline AND using FR in PCA

% PCA_select(classes_only, date_only_cells, All_TC_fit_norm, All_ACF, All_norm_RF,...
%     All_FR, All_norm_DOT, All_norm_RL, color, unique_dates) % normalize TC with fit


% PCA_select(classes_only, date_only_cells, All_TC_spline_norm, All_ACF, All_norm_RF,...
%     [], All_norm_DOT, All_norm_RL, color, unique_dates) % normalize TC with spline

% PCA_select(classes_only, date_only_cells, All_TC_fit_norm, All_ACF, All_norm_RF,...
%     All_norm_DOT, All_norm_RL, color, unique_dates) % normalized TC with fit

keyboard

