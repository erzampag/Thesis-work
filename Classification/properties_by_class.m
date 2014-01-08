function [data, indices, celltypes, average_properties, norm_ave_properties] = properties_by_class(date, datafile)
% This function takes a single retina white noise run that has been
% classified in Vision into ON/OFF, S/M/L, B/S, T/S classes.  It pulls in
% data from the parameters file and averages them by class.  It outputs a cell
% array with class name, DOT, RL, RF, TC, ACF, ACF_SS and an array with those
% properties normalized to OFFLBT.
%
% ENZ, Fall 2012

%%

other_stim = 1;

javaaddpath('/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar');

% Define Paths
neuronfile=strcat('/Users/erinzampaglione/Documents/processed_data/', date, '/', datafile, '/' , datafile, '.neurons');
neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(neuronfile);

param_path = strcat('/Users/erinzampaglione/Documents/processed_data/', date, '/', datafile, '/' , datafile,'.params');
pfile = edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);

sta_path = strcat('/Users/erinzampaglione/Documents/processed_data/', date, '/', datafile, '/' , datafile,'.sta');
stafile = edu.ucsc.neurobiology.vision.io.STAFile(sta_path);

% Get Neuron IDs
idListPreD =  import_neuronIDs(neuronFile);

% Get Classes
classes = import_classes(pfile, idListPreD);

%Get Parameters
data = {};
test = zeros(length(idListPreD), 2);


for i = 1: length(classes);
    
    data{i,1} = classes{i,1}; % neuron ID
    data{i,2} = classes{i,2}; % classes
    
    data{i,3} = pfile.getDoubleCell(idListPreD(i),'dot'); % degree of transience
    data{i,4} = pfile.getDoubleCell(idListPreD(i),'rl'); % response latency
    
    x = pfile.getDoubleCell(idListPreD(i),'SigmaX');
    y = pfile.getDoubleCell(idListPreD(i),'SigmaY');
    data{i,5} = sqrt(x*90*y*90)*2; % calculating RF diameter in microns (assuming these are always 10x10 pixels) 9 microns - 1 pixel
   
    acf(i,:) = pfile.getArrayCell(idListPreD(i),'Auto')';
    acf_norm(i,:) = acf(i,:)/sqrt(sum(acf(i,:).^2)); % normalized ACF (from MK scripts)

    data{i,6} =acf_norm(i,:);
    
 
    
    % Get timefilters.
    blue(:,i) =pfile.getArrayCell(idListPreD(i),'BlueTimeCourse');
    
    if other_stim == 1;
        % Responses to full field flashes

        data{i,8} = pfile.getArrayCell(idListPreD(i), 'flashResponse')';
         data{i,8} = data{i,8}/sqrt(sum(data{i,8}.^2)); % normalized flash response
        
        data{i,9} = pfile.getDoubleCell(idListPreD(i),'flashBinSize'); %ans
        
        % Reversing gratings?
%         data{i,10} = pfile.getDoubleCell(idListPreD(i),'nonLinIndex');
    end
    
end

% % keyboard

% Normalize TCs (from MK scripts) - normally does it for all 3 filters
nrm = sqrt(sum(blue.^2));
blue = bsxfun(@rdivide, blue, nrm)'; % This divides each column of the matrix by corresponding element of the nrm vector (which is a row vector).

for i = 1: length(classes)
    data{i,7} = blue(i,:);
end


% STA depth from MK sta.m
neuronID=data{1,1}(1); % get a token "neuronID"
staObj = stafile.getSTA(neuronID);
STADepth=staObj.getSTADepth;
STARefreshTime=staObj.getRefreshTime; % movie resolution
xlimit = STADepth * STARefreshTime;

% Get time dimension. from MK java_TC.m
[sz_y ~]=size(blue'); 
t=linspace(-xlimit,0,sz_y)'; % time vector


%% make graphs for poster
poster_PCA = 0;
if poster_PCA == 1
    keyboard
    index = [];
    for i = 1:length(idListPreD)
        
        if strcmp(data{i,2}(1:7), 'All/PCA');
            index = [index i];
        end
    end
    index = index';
    
    [coeff_TC,score_TC] = princomp(blue(index, :));
    
    figure
    hold on
    for i = 1: length(index)
        if strcmp(data{index(i),2}, 'All/PCA/OFF/SSS')
            scatter(score_TC(i,1), score_TC(i,2),'ok', 'MarkerFaceColor', 'r');
        else
            scatter(score_TC(i,1), score_TC(i,2),'ok', 'MarkerFaceColor', 'b');
        end
    end
    xlabel('Time Course PC1','FontSize', 20); ylabel('Time Course PC2','FontSize', 20);
    
    
    
    % Other OFF cells
    index_OFF = [];
    indexed_RFs = [];
    for i = 1:length(idListPreD)
        if length(char(data{i,2})) < 15
            continue
        elseif strcmp(data{i,2}(14:15), 'BT');
            index_OFF = [index_OFF i];
            indexed_RFs = [indexed_RFs data{i,5}];
        end
    end
    index_OFF = index_OFF';
    indexed_RFs = indexed_RFs';
    
    [coeff_TC_OFF,score_TC_OFF] = princomp(blue(index_OFF, :));
    
    [coeff_ACF_OFF,score_ACF_OFF] = princomp(acf_norm(index_OFF, :));
    
    
    figure
    hold on
    for i = 1: length(index_OFF)
        if strcmp(data{index_OFF(i),2}, 'All/PCA/OFF/LBT')
            scatter(score_ACF_OFF(i,1), score_TC_OFF(i,1),'ok', 'MarkerFaceColor', 'r');
            hold on
        else
            scatter(score_ACF_OFF(i,1), score_TC_OFF(i,1),'ok', 'MarkerFaceColor', 'b');
            hold on
        end
        %     keyboard
    end
    xlabel('AutoCorrelation PC1','FontSize', 20); ylabel('Time Course PC1','FontSize', 20);
    
    
    % figure
    % hold on
    % for i = 1: length(index_OFF)
    %     if strcmp(data{index_OFF(i),2}, 'All/PCA/OFF/LBT')
    %         scatter(score_ACF_OFF(i,1), score_TC_OFF(i,2),'ok', 'MarkerFaceColor', 'r');
    %         hold on
    %     else
    %         scatter(score_ACF_OFF(i,1), score_TC_OFF(i,2),'ok', 'MarkerFaceColor', 'b');
    %         hold on
    %     end
    % %     keyboard
    % end
    %
    % figure
    %  scatter3(score_TC_OFF(:,1), score_ACF_OFF(:,1), indexed_RFs(:), 150, 'o', 'filled', 'r');
    %
    % figure
    % hold on
    % for i = 1: length(index_OFF)
    %     if strcmp(data{index_OFF(i),2}, 'All/PCA/OFF/LBT')
    %         scatter3(score_TC_OFF(i,1), score_ACF_OFF(i,1), indexed_RFs(i), 150, 'o', 'filled', 'r');
    %         hold on
    %     else
    %         scatter3(score_TC_OFF(i,1), score_ACF_OFF(i,1), indexed_RFs(i), 150, 'o', 'filled', 'b');
    %         hold on
    %     end
    % %     keyboard
    % end
    % xlabel('Time Course PC1','FontSize', 20); ylabel('Time Course PC2','FontSize', 20); zlabel('Receptive Fields', 'FontSize', 20);
    
    
    
    % ON cells
    index_ON = [];
    for i = 1:length(idListPreD)
        if length(char(data{i,2})) < 14
            continue
        elseif strcmp(data{i,2}(13:14), 'BT');
            index_ON = [index_ON i];
        end
    end
    index_ON = index_ON';
    
    [coeff_TC_ON,score_TC_ON] = princomp(blue(index_ON, :));
    [coeff_ACF_ON,score_ACF_ON] = princomp(acf_norm(index_ON, :));
    
    figure
    hold on
    for i = 1: length(index_ON)
        if strcmp(data{index_ON(i),2}, 'All/PCA/ON/LBT')
            scatter(score_ACF_ON(i,1), score_TC_ON(i,1),'ok', 'MarkerFaceColor', 'r');
            hold on
        else
            scatter(score_ACF_ON(i,1), score_TC_ON(i,1),'ok', 'MarkerFaceColor', 'b');
            hold on
        end
        %     keyboard
    end
    xlabel('AutoCorrelation PC1','FontSize', 20); ylabel('Time Course PC1','FontSize', 20);
end
%% Make Indices
ON = [];
OFF = [];
small = [];
medium = [];
large = [];
brisk = [];
sluggish = [];
transient = [];
sustained = [];

% % % keyboard

for i = 1:length(idListPreD) % Creates indices of all assigned parameters (ON, OFF, S, M, L, B, S, T, S)
    if strcmp(data{i,2}(5:6), 'ON')
        ON = [ON i]; % making index
        
        if strcmp(data{i,2}(8), 'S')
            small = [small i];
        elseif strcmp(data{i,2}(8), 'M')
            medium = [medium i];
        elseif     strcmp(data{i,2}(8), 'L')
            large = [large i];
        end
        
        if strcmp(data{i,2}(9), 'B')
            brisk = [brisk, i];
        elseif strcmp(data{i,2}(9), 'S')
            sluggish = [sluggish i];
        end
        
        if strcmp(data{i,2}(10), 'T')
            transient = [transient i];
        elseif strcmp(data{i,2}(10), 'S')
            sustained = [sustained i];
        end
        
    elseif strcmp(data{i,2}(5:7), 'OFF')
        OFF = [OFF i];
        if strcmp(data{i,2}(9), 'S')
            small = [small, i];
        elseif strcmp(data{i,2}(9), 'M')
            medium = [medium i];
        elseif strcmp(data{i,2}(9), 'L')
            large = [large i];
        end
        
        if strcmp(data{i,2}(10), 'B')
            brisk = [brisk i];
        elseif strcmp(data{i,2}(10), 'S')
            sluggish = [sluggish i];
        end
        
        if strcmp(data{i,2}(11), 'T')
            transient = [transient i];
        elseif strcmp(data{i,2}(11), 'S')
            sustained = [sustained i];
        end
    end
end

% intersect(intersect(ON,medium), intersect(brisk, transient));

indices = {}; % Creates indices of all possible celltypes

indices{1} = intersect(intersect(ON,small), intersect(brisk, transient));
indices{2} = intersect(intersect(ON,small), intersect(brisk, sustained));
indices{3} = intersect(intersect(ON,small), intersect(sluggish, transient));
indices{4} = intersect(intersect(ON,small), intersect(sluggish, sustained));

indices{5} = intersect(intersect(ON,medium), intersect(brisk, transient));
indices{6} = intersect(intersect(ON,medium), intersect(brisk, sustained));
indices{7} = intersect(intersect(ON,medium), intersect(sluggish, transient));
indices{8} = intersect(intersect(ON,medium), intersect(sluggish, sustained));

indices{9} = intersect(intersect(ON,large), intersect(brisk, transient));
indices{10} = intersect(intersect(ON,large), intersect(brisk, sustained));
indices{11} = intersect(intersect(ON,large), intersect(sluggish, transient));
indices{12} = intersect(intersect(ON,large), intersect(sluggish, sustained));

indices{13} = intersect(intersect(OFF,small), intersect(brisk, transient));
indices{14} = intersect(intersect(OFF,small), intersect(brisk, sustained));
indices{15} = intersect(intersect(OFF,small), intersect(sluggish, transient));
indices{16} = intersect(intersect(OFF,small), intersect(sluggish, sustained));

indices{17} = intersect(intersect(OFF,medium), intersect(brisk, transient));
indices{18} = intersect(intersect(OFF,medium), intersect(brisk, sustained));
indices{19} = intersect(intersect(OFF,medium), intersect(sluggish, transient));
indices{20} = intersect(intersect(OFF,medium), intersect(sluggish, sustained));

indices{21} = intersect(intersect(OFF,large), intersect(brisk, transient));
indices{22} = intersect(intersect(OFF,large), intersect(brisk, sustained));
indices{23} = intersect(intersect(OFF,large), intersect(sluggish, transient));
indices{24} = intersect(intersect(OFF,large), intersect(sluggish, sustained));

celltypes = {'ONSBT','ONSBS', 'ONSST', 'ONSSS', 'ONMBT', 'ONMBS', 'ONMST', 'ONMSS', 'ONLBT', 'ONLBS', 'ONLST', 'ONLSS',...
    'OFFSBT', 'OFFSBS', 'OFFSST', 'OFFSSS', 'OFFMBT', 'OFFMBS', 'OFFMST' ,'OFFMSS', 'OFFLBT', 'OFFLBS', 'OFFLST', 'OFFLSS'};

%% take averages
% average_properties = zeros(24,5);

% % % % keyboard

average_properties = {}; % for each celltype, calculate the average properties of those neurons
for i = 1 : length(celltypes)
    
    
    average_properties{i,1} = celltypes{i};
    
    if isempty(indices{i})
        continue
    end
    
    degree_transiency = [];
    response_latency = [];
    receptive_field = [];
    %     autocorrelation = {};
    %     time_course = {};
    autocorrelation = [];
    time_course = [];
    
    flash_response = [];
    flash_bin_size = [];
    non_lin_index = [];
    
    
    for j = 1 : length(indices{i}) % temporary list of properties of each neuron of that celltype
        degree_transiency = [degree_transiency data{[indices{i}(j)],3}]; % DOT
        response_latency = [response_latency data{[indices{i}(j)],4}]; % RL
        receptive_field = [receptive_field data{[indices{i}(j)],5}]; % RF
        %         autocorrelation{j} = data{[indices{i}(j)],6};
        %         time_course{j} = data{[indices{i}(j)],7};
        autocorrelation = [autocorrelation; data{[indices{i}(j)],6}];
        time_course = [time_course; data{[indices{i}(j)],7}];

        if other_stim == 1
            flash_response = [flash_response; data{[indices{i}(j)],8}];
            flash_bin_size = [flash_bin_size; data{[indices{i}(j)],9}];
%             non_lin_index = [non_lin_index; data{[indices{i}(j)],10}];
        end
    end
    
% % %     keyboard
    
% % % % % % % % % % % % %     figure
% % % % % % % % % % % % %     plot(time_course', 'b')
% % % % % % % % % % % % %     hold on
% % % % % % % % % % % % %     plot(mean(time_course), 'r')
% % % % % % % % % % % % %     
% % % % % % % % % % % % %     figure
% % % % % % % % % % % % %     hist(degree_transiency, 'b')
% % % % % % % % % % % % %     hist(mean(degree_transiencey,'r'))
% % % % % % % % % % % % %     hold on
% % % % % % % % % % % % %     
    
    
    %     keyboard
    
    if strcmp(celltypes{i}(2),'N') % Average DOT (either positive or negative according to ON or OFF classification)
        average_properties{i,2} = mean(degree_transiency); % times polarity
    else
        average_properties{i,2} = -mean(degree_transiency);
    end
    
    average_properties{i,3} = mean(response_latency);
    average_properties{i,4} = mean(receptive_field);
    average_properties{i,5} = mean(autocorrelation);
    average_properties{i,6} = mean(time_course);
    
    average_properties{i,7} = mean(average_properties{i,5}(101:200));
    
    average_properties{i,11} = median(receptive_field);
    
    if other_stim == 1
        average_properties{i,8} = mean(flash_response);
        average_properties{i,9} = mean(flash_bin_size);
        average_properties{i,10} = mean(non_lin_index);
    end
end
% %    keyboard

% Normalize to OFFLBT (21) - these are now all recalculated in
% class_properties_graph.m
norm_ave_properties = {};
for i = 1: length(average_properties)
    
            norm_ave_properties{i,1} = celltypes{i};
    
    if ~isempty(indices(21)) && ~isempty(average_properties{i,2});
        
        
        norm_ave_properties{i,2} = -average_properties{i,2}./average_properties{21,2}; % DOTxpolarity
        norm_ave_properties{i,3} = -average_properties{i,3}./average_properties{21,3}; % RL
        norm_ave_properties{i,4} = average_properties{i,4}./average_properties{21,4}; % RF
        
        % These need to be normalized differently
        norm_ave_properties{i,5} = average_properties{i,5}./average_properties{21,5}; % ACF
        norm_ave_properties{i,6} = average_properties{i,6}./average_properties{21,6}; % TC
        
        norm_ave_properties{i,7} = average_properties{i,7}./average_properties{21,7}; % ACF steadystate
    end
end


end
