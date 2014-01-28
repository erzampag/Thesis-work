function [idList, AllCells, DSI, DSI_error, ratio, noiseratio,...
    DScellsindex, noise, zerofreq, fwhh, reduced_chi_square, MaxSpikes, EIx, EIy,...
    Fone_Fzero, Ftwo_Fzero, OS_x1, OS_y1, OS_x2, OS_y2, Max_Spikes_absolute, spike_rate_all_orientation]=...
    MultiBarsAnalysis(date, datafile, stimfile2)

 tic
% This function determines receptive field properties of neurons stimulated
% with a drifting square-wave stimulus presented in 16 different directions
% at five trials at each orientation.
% By determining the avererage spike rate at each orientation, a preferred orientation is found
%(direction closest to the vector sum of the average spike rates at each orientation).
% From this direction-selectivity index (DSI), width of tuning curve (fwhh),
% and EI position are determined.

%Erin Zampaglione, Summer 2012
%Adapted from drifting_squarewave_shorter_stimulus2, by Willie Tobin circa early 2010 and modified by Jason Triplett
%10-12-2010, Richard Smith in Fall 2011, and Erin Zampaglione in Winter
%2012

%% get neuron info

[idListPreD, classes, ~, EIx, EIy, data, ~, neuronFile] = import_neuron_info(...
    '/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
    '/Users/erinzampaglione/Documents/processed_data/', date, datafile);

[idList, ind] = remove_duplicates(idListPreD, classes);

EIx = EIx(ind);
EIy = EIy(ind);
data = data(ind);

keyboard
%% Determine how many figures to print:
% 0 - print nothing but SAVE VARIABLES!,
% 1 - print everything,
% 2 - print PSTH, FFT, tuningcurve, vector for neuron,
% 3 - print tuning curve and vector together,
% 4 - just final whole retina analysis
% 5 - von Mises fit and data
% 6 - tuning curve FFT
% 7 - OS von Mises fit
printfigures = 0;

%% Create directories to save figures for given cell types found

A= strcat('/Users/erinzampaglione/Documents/Lab_Work/DSOSCells/', date, '/', datafile);
mkdir(A);

%% Define variables for use in analysis

AllCells= zeros(length(idList), 2);   %these arrays will contain 2 numbers for each DS or OS cell: direction and magnitude.

DSI = zeros(1,length(idList));

ratio = zeros(1,length(idList));
noise = zeros(1,length(idList));
noiseratio = zeros(1,length(idList));
zerofreq = zeros(1,length(idList));

DSI_error = zeros(1,length(idList));

fwhh = zeros(1,length(idList));
% size_confint = [];
% fwhh_fractional_error = [];
reduced_chi_square = zeros(1,length(idList));
% reduced_chi_square_binom = [];

MaxSpikes = zeros(1,length(idList));

Fone_Fzero = zeros(1,length(idList));
Ftwo_Fzero = zeros(1,length(idList));
OS_x1 = zeros(1,length(idList)); OS_y1 = zeros(1,length(idList));
OS_x2 = zeros(1,length(idList)); OS_y2 = zeros(1,length(idList));
OS_counter = 0;
Max_Spikes_absolute = zeros(1,length(idList));


spike_rate_all_orientation = cell(length(idList),1);
%% Load stimulus information

[all_orientations, all_spatial_periods, all_temporal_periods, orientation, spatial, temporal, num_trials, frames] =...
    load_stim_info(date, stimfile2);

%% Use TTL pulses to determine the start and end times of the stimuli

% [start_times, end_times] = calculate_start_end(neuronFile);
% 
   keyboard
 
[start_times, end_times, fourD_in, num_trials_NEW]  = calculate_start_end2(neuronFile, frames, all_spatial_periods, all_temporal_periods,...
    all_orientations, spatial, temporal, orientation, num_trials);

%% Start of the loop over neurons

% Troubleshooting subsets of neurons using either their NeuronIDs or indices

% TROUBLE = [180;185;188;189;195;204;211;214;216;219;221;222;254;261;263;264;272;277;278;281;282;291;295]; %2012-02-23-PAPER FIGURE IN THIS GROUP!!

% TROUBLE = [20; 116; 128;138;144;153;155;175;178;191;194;195;203;214;227;228;239;244;251;260;267;271;285;287;377;443;444;448;...
%     482;516;528;530;532;]; % OS from FFT 2010-08-25
% TROUBLE = [6;16;19;33;42;47;48;56;57;61;64;65;68;69;71;75;86;87;89;90;94;99;101;110;111;112;129;130;142;151;157;162;163;170;171;...
%     180;185;188;189;195;197;204;210;211;214;216;219;221;229;254;261;262;263;264;272;277;278;281;282;291;295;304;305;308;313;317;318;...
%     319;324;334;335;347;348;349;353;355;357;358;360;368;372;373;376;380;395;404;414;418;420;458;484;494;516;528;536;]; % DS from FFT 2011-02-23
TROUBLE = [7;16;33;47;53;54;55;57;59;65;68;69;71;83;86;87;89;90;97;100;111;130;133;142;151;157;170;179;180;183;188;194;195;197;205;...
    211;213;223;237;254;263;264;272;280;281;282;294;295;304;305;307;317;319;324;334;347;348;349;353;354;355;357;358;368;370;372;373;377;...
    380;395;414;418;420;433;435;436;484;490;494;501;503;515;518;528;529;]; % OS from FFT 2011-02-23
%     TROUBLE = [25;33;66;74;85;94;98;105;109;118;122;124;130;140;147;165;170;177;187;217;219;235;259;273;323;325;331;332;334;338;363;372;378;...
%         386;394;402;412]; % OS from FFT 2011-06-28
%     TROUBLE = [32;45;52;61;62;68;84;88;92;112;168;200;223;224;249;254;260;...
%        277;307;308;324;330;351;390;392;407]; % OS from FFT 2011-03-25(DSCAM-/-)

% TROUBLE = [812;961;980;1277;3121;3393;3933;4187;4621;5088;6001;7667]; % 2011-08-24 neuron IDs of OS cells
% TROUBLE = [54 55];

% TROUBLE = [2, 3];
% TROUBLE = 50;
% TROUBLE = 27;

TROUBLE = [3;15;16;23;30;31;34;37;39;50;67;68;70;79;120;194;213;225;...
    241;246;247;252;264;275;286;303;331;346;356;363;366;369;]; % DS cells from DS analyis RABBIT, 2013-06-10

TROUBLE = [1;257;287;348;351;527;635;677;828;872;931;964;1098;1142;1443;1548;1908;1940;2011;2044;2116;2194;2388;2418;...
    2674;2747;2792;2793;2838;2993;3124;3182;3202;3482;3528;3586;3740;4384;4413;4520;4741;4847;5192;5257;5267;5342;5462;...
    5582;5642;5761;5783;5837;5881;6078;6406;6618;6663;6681;6694;6811;6857;6860;7082;7130;7156;7234;7323;7671;]; % DS by hand RABBIT, 2013-06-10


fprintf('Analyzing Retina from %s \n', date)
for q=1:length(idList);
% for q = 1: length(TROUBLE);
    
    if mod(q,50)==0,
        fprintf('processing neuron # %d \n',q); 
    end
    
    NeuronID=idList(q); % complete run-thru
    
%     NeuronID = idList(TROUBLE(q)); % if TROUBLE is the indices
     
%     NeuronID = TROUBLE(q); % if TROUBLE is the neuron IDs
    
    temp=neuronFile.getSpikeTimes(NeuronID);
    spikeTimes=double(temp);%converts the int matrix "spikeTimes" to a double matrix
%     spikeTimes = spiketimes{NeuronID};

    if length(spikeTimes) < 100
        fprintf('Neuron ID %d had only %d spikes, and we are gonna skip it\n', NeuronID, length(spikeTimes))
         continue
    end
    
    %% Creation of the raster plots, frequency (PSTH) plots for each orientation, in new figure window
    % Output PSTH for FFT and total spikes for finding rate / dir pref
    % If you want grey screen spike rates, use sort_spikes_by_stimulus_grey     

% % %     keyboard

%     [total_spikes, combined_spikes, segments] = sort_spikes_by_stimulus(orientation, all_orientations, spikeTimes,...
%         start_times, end_times, printfigures);
   
%     
    [total_spikes, combined_spikes, segments] = sort_spikes_by_stimulus2(orientation, temporal, spatial, frames, fourD_in, ...
        num_trials, num_trials_NEW, spikeTimes, start_times, end_times, printfigures);

    total_spikes = squeeze(total_spikes);
    combined_spikes = squeeze(combined_spikes);
% % %      keyboard
    %% Average and Normalize Spikes rates
    
    
    [true_spike_rate, norm_true_spike_rate, Max_Spikes_absolute(q), max_index] =...
        ave_and_norm_spike_rates(orientation, all_orientations, num_trials, start_times, end_times, total_spikes, q);

% 
%     [true_spike_rate, norm_true_spike_rate, Max_Spikes_absolute, max_index] =...
%     ave_and_norm_spike_rates2(spatial, temporal, orientation, num_trials, num_trials_NEW, frames, total_spikes, printfigures);

    true_spike_rate = squeeze(true_spike_rate);
    norm_true_spike_rate = squeeze(norm_true_spike_rate);
    %%%neeed to make maxspikesabsolute to be a cell array
    spike_rate_all_orientation{q} = squeeze(true_spike_rate); % saving all true spike rates for later analysis

    
    %% Calculation of R and theta
    %
    %     if spikeSumx<0
    %         AllCells(q,1)=atan(norm_spikeSumy/norm_spikeSumx)+pi;
    %         %the theta value for the total preferred direction for each cell
    %
    %     else
    %         AllCells(q,1)=atan(norm_spikeSumy/norm_spikeSumx);  %(necessary to get correct sign direction)
    %
    %     end
    
    %     spikeSumx=0;  %to be used in adding up the vectors for different direction-selectivities
    %     spikeSumy=0;
    norm_spikeSumx = 0;
    norm_spikeSumy = 0;
    
    for j = 1:length(orientation);
        norm_spikeSumx = norm_spikeSumx + norm_true_spike_rate(j)*(cos(orientation(j)*(2*pi/360)));
        norm_spikeSumy = norm_spikeSumy + norm_true_spike_rate(j)*(sin(orientation(j)*(2*pi/360)));
    end
    % % %     keyboard
    
    
    % R value for the total preferred direction for each cell
    AllCells(q,2)=((norm_spikeSumx^2)+(norm_spikeSumy^2))^(.5);
    
    % Theta value for total preferred direction for each cell
    AllCells(q,1) = atan2(norm_spikeSumy,norm_spikeSumx);
    
    
    %Convert this theta to degrees temporarily
    temp_theta_degrees = [];
    
    if AllCells(q,1) < 0
        temp_theta_degrees = (AllCells(q,1) + 2*pi) * 180/pi;
    else
        temp_theta_degrees = AllCells(q,1) * 180/pi;
    end
    
    %% New Determination of preferred orientation!
    difference_vector = zeros(1,length(orientation));
    for i = 1 : length(orientation);
        difference_vector(i) = abs(orientation(i) - temp_theta_degrees);
    end
    
    min_diff_vect = min(difference_vector);
    index_of_min = find(difference_vector <= min_diff_vect);
    pref_O = orientation(index_of_min); % the pref_O is based off the Elstrott paper
    
   
    try
    pref_O = pref_O(1);
    catch
        q
    end
    
    index_of_opp = mod((index_of_min+7), 16) + 1;
    
    
    orth_O=pref_O+90;
    if orth_O > 337.5;
        orth_O = orth_O-360;
    end
    
    orth_O_two = pref_O+270; % for ASI
    if orth_O_two > 337.5
        orth_O_two = orth_O_two - 360;
    end
    
    opp_O=pref_O+180;
    if opp_O > 337.5;
        opp_O = opp_O-360;
    end
    idpref_O=find(orientation==pref_O);
    pref_Orate=true_spike_rate(idpref_O);
    idorth_O=find(orientation==orth_O);
    orth_Orate=true_spike_rate(idorth_O);
    idopp_O=find(orientation==opp_O);
    opp_Orate=true_spike_rate(idopp_O);
    
    idorth_O_two=find(orientation==orth_O_two);
    orth_Orate_two=true_spike_rate(idorth_O_two);
    
    %Calculate orientation selectivity and direction selectivity
    %(basically just the percent error between 2 values for each one)
    OSI = (pref_Orate - orth_Orate)/(pref_Orate + orth_Orate);
    OSI = round(OSI/0.01)*.01;
    DSI(q) = (pref_Orate - opp_Orate)/(pref_Orate + opp_Orate);
    %     DSI = round(DSI/0.01)*.01;
    pref_Arate = pref_Orate + opp_Orate;
    orth_Arate = orth_Orate + orth_Orate_two;
    ASI(q) = (pref_Arate - orth_Arate)/(pref_Arate + orth_Arate);
    
    MaxSpikes(q) = pref_Orate*50;
    %% Fit tuning curve, print figure if ==5
    
    [fwhh(q), reduced_chi_square(q)] = fit_tuning_curve(pref_Orate, pref_O, temp_theta_degrees, true_spike_rate,...
        orientation, NeuronID, q, printfigures, AllCells);
    
    % % %     keyboard
    %% Choosing the best stimulus to determine if it is ON/OFF or just ON or OFF
    
    if DSI(q) < 0.5
        thebeststimulus = max_index; % anything that is not a DS cell
    else
        thebeststimulus = index_of_min; % for DS cells, we want to look at the direction closest to the vector sum
    end

%     [ratio(q), noise(q), zerofreq(q), noiseratio(q)] = fourier_transform_PSTH(segments, thebeststimulus, combined_spikes, q);
%     keyboard
    
    [ratio(q), noise(q), zerofreq(q), noiseratio(q)] = fourier_transform_PSTH2([segments{1}], frames,...
                thebeststimulus, orientation, combined_spikes, temporal, spatial, printfigures);
        
    

    %% Calculation of DSI error
    % I did propagation of errors on the DSI by assuming a gaussian and poisson
    % distribution of the total number of spikes (summed)
    
    %     p_spikes = sum(combined_spikes{index_of_min}.*5);
    %     n_spikes = sum(combined_spikes{index_of_opp}.*5);
    
    %     p_spikes = unave_true_spike_rate(index_of_min)*10;
    %     n_spikes = unave_true_spike_rate(index_of_opp)*10;
    
    p_spikes = true_spike_rate(index_of_min)*50;
    n_spikes = true_spike_rate(index_of_opp)*50;
    
    DSI_error(q) = (2/(p_spikes+n_spikes)^2)*(n_spikes*p_spikes*(n_spikes+p_spikes))^(1/2);
    
    %% calculation of DS and OS by FFT
    Y=fft(true_spike_rate);
    n=length(Y);
    % % if n>1;
    % %     Y(1)=[];
    % % end
    % % power = abs(Y(1:floor(n/2))).^2
    power_tuning = abs(Y(1:ceil(n/2)+1)).^2;
    
    amp_tuning = abs(Y(1:ceil(n/2)+1));
    % nyquist = 1/2;
    % % freq = (1:n/2)/(n/2)*nyquist*(length(segments)/10);
    %
    
    % % %     Fone_Fzero(q) = power_tuning(2)/power_tuning(1);
    % % %     Ftwo_Fzero(q) = power_tuning(3)/power_tuning(1);
    
    Fone_Fzero(q) = amp_tuning(2)/amp_tuning(1);
    Ftwo_Fzero(q) = amp_tuning(3)/amp_tuning(1);

    
    % Only use this if statment if you want to fit the tuning curve to only
    % OS cells as defined by the statment
    %     if Ftwo_Fzero(q) > 0.03 && Ftwo_Fzero(q) > Fone_Fzero(q) && true_spike_rate(max_index)*50 > 100
    %     if Ftwo_Fzero(q) > sqrt(0.03) && Ftwo_Fzero(q) > Fone_Fzero(q) && Max_Spikes_absolute > 100
    
%     OS_counter = OS_counter +1;
    
    if Ftwo_Fzero(q) > Fone_Fzero(q)
        % fits a double Von Mises curve and gives angle of largest Rmax
        mu = fit_tuning_curve2(pref_Orate, pref_O, temp_theta_degrees, true_spike_rate, orientation, NeuronID, q,...
            printfigures, AllCells, Ftwo_Fzero, Fone_Fzero);
        
        mu = mod(mu,2*pi);
        
        % finds stim orientation closest to mu
        pref_OS_dir = orientation(find(abs(orientation.*pi/180-mu) == min(abs(orientation.*pi/180-mu))));
        
        % Define the two sides
        OS_side_one = mod(pref_OS_dir + 90, 360);
        OS_side_two = mod(pref_OS_dir - 90, 360);
        
        if OS_side_one > OS_side_two
            OS_heavy = OS_side_one;
            OS_light = OS_side_two;
        else
            OS_heavy = OS_side_two;
            OS_light = OS_side_one;
        end
        
        OS_x1temp = 0; OS_y1temp = 0; OS_x2temp = 0; OS_y2temp = 0;
        
        
        for j = 1 : 16
            if orientation(j) >= OS_light && orientation(j) < OS_heavy
                OS_x1temp = OS_x1temp + norm_true_spike_rate(j)*(cos(orientation(j)*(pi/180)));
                OS_y1temp = OS_y1temp + norm_true_spike_rate(j)*(sin(orientation(j)*(pi/180)));
                
            else
                OS_x2temp = OS_x2temp + norm_true_spike_rate(j)*(cos(orientation(j)*(pi/180)));
                OS_y2temp = OS_y2temp + norm_true_spike_rate(j)*(sin(orientation(j)*(pi/180)));
            end
        end
        
        OS_x1(q) = OS_x1temp;
        OS_x2(q) = OS_x2temp;
        OS_y1(q) = OS_y1temp;
        OS_y2(q) = OS_y2temp;
        
    end

    % draw red line
    if printfigures == 7 && Ftwo_Fzero(q) > sqrt(0.03) && Ftwo_Fzero(q) > Fone_Fzero(q) && max(true_spike_rate)*50 > 100
        gcf;
        hold on;
        h = polar([pref_OS_dir*pi/180-pi/2;pref_OS_dir*pi/180+pi/2],[Max_Spikes_absolute(q);Max_Spikes_absolute(q)], 'r-');
        set(h,'linewidth',4);
        
        % draw vector
        subplot(1,2,2);
        
        x_fake=[0 1 0 -1];
        y_fake=[1 0 -1 0];
        
        h_fake=compass(x_fake,y_fake, '--');
        hold on;
        set(h_fake,'Visible','off');

        ch = compass(OS_x1(q), OS_y1(q), '-k');
        xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
        set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
        set(ch,'linewidth',4);
        ch = compass(OS_x2(q), OS_y2(q), '-k');
        xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
        set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
        set(ch,'linewidth',4);
    end
  
%     
%     max_fourier = Y(2);
%     
%     r_TC = abs(Y(2));
%     theta_TC = angle(Y(2));
%     
%     figure
%     x_fake=[0 1 0 -1];
%     y_fake=[1 0 -1 0];
%     
%     h_fake=compass(x_fake,y_fake, '--');
%     hold on;
%     ch = compass((r_TC*cos(theta_TC)), (r_TC*sin(theta_TC)), '-k');
%     xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
%     set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
%     % set(ch,'linewidth',4);
%     set(ch,'LineWidth',2)
%     %          polar2(AllCells(i,1), AllCells(i,2),[0 1], '-k'); % this is in the Matlab folder
%     hold on
%     title(strcat('Direction Selective Vector for neuron-', num2str(NeuronID)));
%     %         end
%     %     end
%     set(h_fake,'Visible','off')
%     keyboard
   
    
    %% Printing Figures!
    % FFT, single polar plot, etc - based on printed_figures switch
    if printfigures ~=0
        stand_alone_printed_figures(printfigures, combined_spikes, segments,...
            orientation, true_spike_rate, NeuronID, AllCells, index_of_min, q, Fone_Fzero, Ftwo_Fzero, power_tuning, thebeststimulus);
    end
% %        keyboard
end
 toc % end of run through all neurons in retina
% keyboard
%% START OF WHOLE RETINA ANALYSIS
if printfigures == 4
    full_retina_printed_figures(Fone_Fzero, Ftwo_Fzero, AllCells, idList, DSI, ratio, noiseratio,...
    MaxSpikes, Max_Spikes_absolute, OS_x1, OS_y1, OS_x2, OS_y2);
end


% DScellsindex = [];
DScellsindex = find(DSI(:) > 0.5 & noiseratio(:) >5 & MaxSpikes(:) > 100);
DScellsindex_ONorOFF = find(DSI(:) > 0.5 & noiseratio(:) >5 & MaxSpikes(:) > 100 & ratio(:) > 1);
    
% keyboard
DS_classtextfile(idList, DScellsindex, ['/Users/erinzampaglione/Documents/processed_data/', date, '/',...
    datafile, '/DS_classification_matlab', datestr(now,30), '.txt']);

%% 
keyboard
        figure
        x_fake=[0 1 0 -1];
        y_fake=[1 0 -1 0];
        
        h_fake=compass(x_fake,y_fake, '--');
        hold on;
        
        for j = 1 : length(DScellsindex)
                      
                ch = compass((AllCells(DScellsindex(j),2).*cos(AllCells(DScellsindex(j),1))),...
                    (AllCells(DScellsindex(j),2).*sin(AllCells(DScellsindex(j),1))),'-k');
 
            xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
            set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
            
            
            %Removing the label
            % set(findall(gcf, 'String', '30','String','60', 'String', '180') ,'String', ' ')
                        set(ch,'LineWidth',1.5)
            hold on
        end
       
        
        for j = 1 : length(DScellsindex_ONorOFF)
            
            ch = compass((AllCells(DScellsindex_ONorOFF(j),2).*cos(AllCells(DScellsindex_ONorOFF(j),1))),...
                (AllCells(DScellsindex_ONorOFF(j),2).*sin(AllCells(DScellsindex_ONorOFF(j),1))), '-b'); % ON or OFF cells
            
            xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
            set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
            set(ch,'LineWidth',1.5)
                        
            
        end
        
        
        ch = title(strcat('Pref Dir/Mag for DS Cells in Retina-',...
            date),'FontSize', 25);
        
%         delete(findall(gcf, 'type', 'text'));
        
        set(h_fake,'Visible','off')
        
        set(gcf, 'PaperPosition', [0,0,12,12]);
print('-depsc', [A '/' 'DS_PolarPlot'])

% each tuning curve
        for q = DScellsindex'
            figure
            
            p = polar([orientation(:)*pi/180'; orientation(1)],...
                [spike_rate_all_orientation{q}; spike_rate_all_orientation{q}(1)]);
            
            title(['Spike Rates for Neuron ' num2str(idList(q)) ', Index ' num2str(q)],'FontSize', 25)
            
            if ratio(q) <1
                set(p, 'Color', 'r')
            end
            
            set(gcf, 'PaperPosition', [0,0,12,12]);
            print('-depsc', [A '/' 'Tuning' num2str(idList(q))])
            
            close(gcf)
        end
        %   

        
      %%  

keyboard
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % keyboard
if printfigures == 0;
    %     filename = strcat(A, '_', datestr(now, 30)); % maybe you don't want this datestr for now
    filename = A;
    save(filename, 'date', 'idList', 'AllCells', 'DSI', 'OSI','DSI_error', 'ratio',...
        'noiseratio', 'DScellsindex', 'noise', 'zerofreq', 'fwhh',...
        'reduced_chi_square', 'MaxSpikes', 'EIx', 'EIy',...
        'Fone_Fzero','Ftwo_Fzero', 'OS_x1', 'OS_y1', 'OS_x2', 'OS_y2', 'Max_Spikes_absolute', 'spike_rate_all_orientation');
end
   keyboard
