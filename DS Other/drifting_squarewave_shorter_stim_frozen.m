function [idList, AllCells, DSI, DSI_error, ratio, noiseratio,...
    DScellsindex, noise, zerofreq, fwhh, reduced_chi_square, MaxSpikes]=...
    drifting_squarewave_shorter_stimulus2(type, genotype, date, datafile, stimfile2)


%This function determines receptive field properties of neurons stimulated
%with a drifting square-wave stimulus presented in 16 different directions
% at five trials at each orientation.
%By determining the avererage spike rate at each orientation, a preferred orientation is found
%(direction which elicits max spike rate). From this orientation selectivity index (OSI) and
%direction-selectivity index (DSI) are determined.

%Written by Willie Tobin circa early 2010 and modified by Jason Triplett
%10-12-2010, Richard Smith in Fall 2011, and Erin Zampaglione in Winter
%2012

javaaddpath('C:\Users\Arash\workspace\vision8\Vision.jar');

% %define path to raw data file
% full_path=strcat('/Volumes/TEMP_MOUSE/', date, '/', datafile, '/' , datafile, '.bin');
% rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

% spikefile=strcat('/Volumes/TEMP_MOUSE/', date, '/', datafile, '/' , datafile, '.spikes');
% spikeFile=edu.ucsc.neurobiology.vision.io.SpikeFile(spikefile);

neuronfile=strcat('I:\', type, '\', genotype, '\', date, '\', datafile, '\' , datafile, '.neurons');
neuronFile=edu.ucsc.neurobiology.vision.io.NeuronFile(neuronfile);

param_path = strcat('I:\', type, '\', genotype, '\', date, '\', datafile, '\' , datafile,'.params');
pfile = edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);

% *****PRINTFIGURES:
% 0 - print nothing but SAVE VARIABLES!,
% 1 - print everything,
% 2 - print PSTH, FFT, tuningcurve, vector for neuron,
% 3 - print tuning curve and vector together,
% 4 - just final whole retina analysis
% 5 - von Mises fit and data
printfigures = 4;

%Create directories to save figures for given cell types found
A= strcat('C:\Users\Arash\My Documents\DSAnalysis\', type, '\', genotype, '\', date, '\', datafile);
mkdir(A);
% B= strcat('/Users/erinzampaglione/Documents/Second term 1styear/DSCells', date, '/', datafile, '/OS');
% mkdir(B);
% C= strcat('/Users/erinzampaglione/Documents/Second term 1styear/DSCells',date, '/', datafile, '/DS');
% mkdir(C);
% D= strcat('/Users/erinzampaglione/Documents/Second term 1styear/DSCells', date, '/', datafile, '/NonOS');
% mkdir(D);

OSCell=0;
DSCell=0;
NonOSCell=0;
AllCells=[];   %these arrays will contain 2 numbers for each DS or OS cell: direction and magnitude.
DSCells=[];
OSCells=[];
DSCellCounter=0;
OSCellCounter=0;
peakcounter=0;
DSI = [];

fwhh = [];
size_confint = [];
fwhh_fractional_error = [];
reduced_chi_square = [];
reduced_chi_square_binom = [];

MaxSpikes = [];

%Import list of neurons IDs from the file that is defined above as
%'neuronFile';
idListbefore = [];
idListbefore=neuronFile.getIDList();

%Re make the vector of neuronIDs so it doesn't have weird negative numbers
idListPreD = [];

for i = 1:length(idListbefore)
    if idListbefore(i) <=0,
        continue
    else
        idListPreD = [idListPreD, idListbefore(i)];
    end
end
idListPreD = idListPreD';

% STUFF OUT OF THE PARAMETERS FILE (CLASSES ETC)
%%%% Get Classes
classInfo = pfile.getClassIDs(); % returns a JavaHashmap Object

classes = cell(length(idListPreD), 2); % each row is [ID ClassInfo]
for i = 1:length(idListPreD)
    classes{i,1} = idListPreD(i); % because it's a cell it needs a curly bracket
    classes{i,2} = classInfo.get(int32(idListPreD(i)));
end

%% Get Data

data = zeros(length(idListPreD), 5);
% 
% for i = 1:length(idListPreD)
%     data(i,1) = idListPreD(i);
%     data(i,2) = pfile.getDoubleCell(idListPreD(i), 'xOffDS');
%     data(i,3) = pfile.getDoubleCell(idListPreD(i), 'yOffDS');
%     data(i,4) = pfile.getDoubleCell(idListPreD(i), 'magOffDS');
%     data(i,5) = pfile.getDoubleCell(idListPreD(i), 'angOffDS');
% end


%% Make IDlist have only removed duplicates
Dup = [];
for i = 1:length(idListPreD)
    if strcmp(classes{i,2}, 'All/Duplicates')
        Dup = [Dup i]; % making a dupilicates index
    end
end

%Different way of making mag index -
% mag4ind = find(data(:,4) >=0.4); % making an index based on a minimum magnitude

ind = setdiff(1:length(data), Dup); % indices of neurons after duplication removal
% ind4 = setdiff(mag4ind, Dup); % indicies of min. mag. neurons after dup removal


% Use indices of neurons to make a duplicate-removed list of neurons

idList = [];
idList = idListPreD(ind);  % Duplicates only
% idList  = idListPreD(ind4); % Duplicates AND Vision-determined DS cells


%% %load stimuli information from a list of orientations displayed for a
% %particular experiment in proper order from a .txt file in Matlab directory

stim_file = fopen(strcat('I:\', type, '\', genotype, '\', date, '\stimuli\',stimfile2)); % for runthru
% stim_file = fopen(strcat('/Volumes/TEMP_MOUSE/', '2010-08-27-0', '/stimuli/', 's07.txt'));
stimFile = textscan(stim_file, '%s%d%s%d%s%f%s', 'headerlines', 1);
stim = stimFile{6};
% spatial_periods = stimFile{2};
% temporal_periods = stimFile{4};


%% create a matrix of start times and end times in sample #s corresponding
% to when stimuli were being shown based on TTL pulses
% add
temp2=neuronFile.getTTLTimes();% create matrix of TTL pulse times
TTL=double(temp2);%make TTL pulse matrix double

k=0;
stimulus_number=0;
End_times=[];
start_times=[];
TTLdiff = TTL(3)-TTL(2); % was TTL(2) -TTL(1)

for i = 2:length(TTL)
    TTLcheck = TTL(i)-TTL(i-1);
    
    if le(TTLcheck,(TTLdiff + 3000)); % if the time between 2 TTL pulses is "normal", carry on
        %(what does the "k" do?)
        k=k+1;
    else                          %otherwise, must mean that a stimulus has just ended
        if stimulus_number == 0;  % if it's the FIRST stimulus that's just ending....
            stimulus_number = stimulus_number+1;
            start_times(stimulus_number) = 0;
            End_times(stimulus_number) =  200000;  %was 21337;
            %=TTL(i-1)+7967; (for old stimulus length of 10 sec)
            %was +33334,this yielded a wrong stimulus duration of 11.26
            start_times(stimulus_number+1) = TTL(i);
            
        else
            if stimulus_number < 78;   %if it's not the last stimulus.. (shouldn't this be 5*16-1=79???)
                stimulus_number = stimulus_number+1;
                End_times(stimulus_number) =  TTL(i-10)+200000; % was TTL(i-1)+21337
                %=TTL(i-1)+7967;(for old stimulus length of 10 sec)
                %was +33334, yielded a wrong stimulus duration of 11.26;
                start_times(stimulus_number+1) = TTL(i);
            else
                stimulus_number = stimulus_number+1;   %if it IS the last stimulus.....
                End_times(stimulus_number) = TTL(i-10)+200000; %was TTL(i-2)+3*(20000)
                %=TTL(i-1)+7967;(for old length of 10 sec)
                %was +33334, which yielded wrong stim duration of 11.26;
                start_times(stimulus_number+1) = TTL(i);
                End_times(stimulus_number+1) =TTL(length(TTL)-10)+200000;  % was  TTL(i-1)+21337
                % TTL(length(TTL)-1)+3*(20000);  %=TTL(length(TTL))+7967;(for old stim length of 10 sec)
                %=TTL(length(TTL))+7967;(for old stim length of 10 sec)
                %was +33334, this yielded wrong stim duration of 11.26;
            end
        end
    end
end
start_times=start_times';
End_times=End_times';

%% Start of the loop through neurons

%for given neuron, import the spikes times, just doing one right now so
%there's not too many plots popping up.

% TROUBLE = [78;215;]; % DSI > 0.5 2011-08-25, data002
% TROUBLE = [68;391;]; % DSI > 0.5 2011-02-08, data002
% TROUBLE = [126;243;345;]; % DS neurons for DSI > 0.5 in 2011-03-25, data003
%  TROUBLE = [6;16;19;47;64;65;68;71;75;86;89;94;99;101;110;112;127;129;130;151;157;163;171];
% 2011-02-23-0 data002
% TROUBLE = [40;116;124;127;181;187;232;250;271;347;]; % DS cells for DSCAM 2011-08-24, data002
% TROUBLE = [9;75;107;125;144;168;171;180;188;204;226;228;238;250;256;283;285;298;348;454;468;474;514;];
%DS cells for 2011-06-28-0 data002 (weird chi2)
% TROUBLE = [21;22;29;36;37;47;137;140;146;176;187;196;207;216;229;258;265;271;307;308;315;319;321;326;338;339;]; 
%DS 2012-04-26 WT (weird chi-square)
% TROUBLE = [86;90;143;228;337;398;]; % 2012-04-30 WT data000 (weird chi-square)

% TROUBLE = [40;124;127;187;232;250;271;347;]; % 2011-08-24 DSCAM stingent cut
% TROUBLE = [126;243;344;345;]; %2011-03-25 DSCAM stringent cut data003
% TROUBLE = [68;391;]; % 2011-02-08 DSCAM strindgent cut data002
% TROUBLE = [78;215;]; % 2011-08-25 DSCAM stringent cut data002
TROUBLE = [80;144;148;154;272;]; % 2011-03-28 DSCAM data 003 stringent cut
% TROUBLE = [6;16;19;47;64;65;68;71;75;86;89;94;99;101;110;112;127;129;130;151;157;163;171;...
%     180;185;188;189;195;204;211;214;216;219;221;222;254;261;263;264;272;277;278;281;282;291;295;...
%     304;305;313;317;318;319;334;347;348;349;355;357;358;360;368;371;372;376;380;404;418;443;458;463;484;494;516;]; % 2011-02-23 WT data002 stringent
% TROUBLE = [5;20;23;37;44;80;87;94;108;112;123;151;155;162;169;172;180;191;198;275;285;294;296;301;325;331;350;361];...
%     373;380;404;410;418;449;452;455;476;486;487;496;498;502;504;512;519;527;539;553;563;567;579;581;583;585;604;609;...
%     649;659;661;678;683;697;699;702;714;721;734;742;743;755;776;790;792;793;821;824;837;844;848;875;];

fprintf('Analyzing Retina from %s \n', date)
for q=1:length(idList);
% for q = 1: length(TROUBLE);
    
    if mod(q,50)==0,
        fprintf('processing neuron # %d \n',q);
        
    end
    %     if idList(q) <= 0,
    %         continue % skips current value, goes into next loop
    %     end
    
    
    
        NeuronID=idList(q);  %% i just chose a DS looking cell!
%     
%     NeuronID = idList(TROUBLE(q)); % if TROUBLE is the indices
    
    %    NeuronID = TROUBLE(q); % if TROUBLE is the neuron IDs
    
    
    temp=neuronFile.getSpikeTimes(NeuronID);
    spikeTimes=double(temp);%converts the int matrix "spikeTimes" to a double matrix
    
    
    % plot spike times as hash marks on figure subplots, and form an array of
    % arrays of arrays of spike times for each orientation, trial.
    %     fig1 = figure(1);
    Direction=0:22.5:337.5;
    orientation=[0:22.5:337.5];
    allspikes={};
    combinedspikes={};
    extractedspikes={};
    
    %% creation of the raster plot
    if printfigures == 1
        figure
    end
    %     if printfigures == 2
    % figure
    %     end
    for i=1:16;
        %         % These graph-things were turned off
        
        if printfigures == 1
            
            axes('position', [(mod((i-1),4)/4)+.05, ((4-ceil((i*4)/16))/4)+.05, .15, .15 ]);
            axis([-3,13,0,5]);
            title(Direction(i),'FontSize',9,'FontWeight','bold');
            
            grid ON;
        end
        trial_count=0;
        
        for j=1:length(stim);
            if stim(j)==orientation(i);
                trial_count = trial_count+1;
                spikes=[];
                counter=0;
                counter2=0;
                for l=1:length(spikeTimes);
                    if spikeTimes(l)>start_times(j)-60000 && spikeTimes(l)<End_times(j)+60000;
                        %spikeTimes(l)>start_times(j) && spikeTimes(l)<End_times(j);
                        xcoord=((spikeTimes(l)-(start_times(j)))/20000);
                        counter = counter+1;
                        x=[xcoord,xcoord];
                        y=[trial_count,trial_count-1];
                        spikes(counter)=xcoord;
                        
                        if printfigures == 1
                            line(x,y); % turn off to not graph
                        end
                        
                        
                        if spikeTimes(l)>start_times(j) && spikeTimes(l)<End_times(j);  %so that even
                            %though the hash marks display before and after the stimulus,...
                            %the fourier transform will be just of the data
                            counter2=counter2+1;            % that happened during the stimulus
                            allspikes{i, trial_count}(counter2)=spikes(counter);
                            
                        end
                        
                    end
                end
                if counter2==0;
                    allspikes{i, trial_count}=[];
                    
                end
            end
        end
        if printfigures == 1
            hold on;
        end
        %         if printfigures == 2
        %             hold on;
        %         end
        
    end
    if printfigures == 1
        hold off;
    end
    
    
    
    %% Use the allspike arrays to make frequency plots for each orientation,
    %in a new figure window - the PSTH
    
    segments=[.05:.1:9.95];
    %segments = 100;
    minpeakdist=.04*segments;   %to ensure that the peaks aren't right next to each other
    if printfigures == 1
        figure
    end
    
    for i=1:16;
        
        %                 combinedspikes{i}=(1/2)*(hist(allspikes{i,1},segments) + hist(allspikes{i,2},segments));
        %%%THIS IS FOR ANOMALOUS RUN!
        
        %         combinedspikes{i}=(1/3)*(hist(allspikes{i,3},segments) + hist(allspikes{i,4},segments)...
        %             + hist(allspikes{i,5},segments));   %%%THIS IS FOR ANOMALOUS RUN!
        
        combinedspikes{i}=(1/5)*(hist(allspikes{i,1},segments) + hist(allspikes{i,2},segments)...
            + hist(allspikes{i,3},segments) + hist(allspikes{i,4},segments)...
            + hist(allspikes{i,5},segments));   %%%%%THIS IS WHAT YOU NORMALLY USE
        
        if printfigures == 1
            axes('position', [(mod((i-1),4)/4)+.05, ((4-ceil((i*4)/16))/4)+.05, .15, .15 ]);
            
            plot(segments, combinedspikes{i}) % for vector segments
            % plot(combinedspikes{i})
            title(Direction(i),'FontSize',9,'FontWeight','bold');
            xlabel('combined spikes')
            ylabel('')
            hold on
        end
    end
    if printfigures == 1
        hold off;
    end
    
    %% FFT
    if printfigures == 1
        figure
    end
    for i=1:16;
        if printfigures == 1
            axes('position', [(mod((i-1),4)/4)+.05, ((4-ceil((i*4)/16))/4)+.05, .15, .15 ]);
        end
        
        
        Y=fft(combinedspikes{i});
        n=length(Y);
        if n>1;
            Y(1)=[];
        end
        power = abs(Y(1:floor(n/2))).^2;
        nyquist = 1/2;
        freq = (1:n/2)/(n/2)*nyquist*(length(segments)/10); % when segments = [0.1:0.1:10]
        %             freq = (1:n/2)/(n/2)*nyquist*(segments/10); % when
        %             segments = 100
        %          %the segments/10 is to transform the
        %         units from "spikes per segment" to hertz, since the spiketimes go from 0 to 10 seconds
        
        if printfigures == 1
            plot(freq,power);
            title(Direction(i),'FontSize',9,'FontWeight','bold');
            %
            %
            % %loglog(freq,power)   % should i turn this one?!
            xlabel('frequency (Hz)'); set(gca,'XTick',0:1:5)
            ylabel('amplitude')
            
            % axis([0, 50, 0, 2000]);
            hold on;
        end
    end
    if printfigures == 1
        hold off;
    end
    
    %
    %% Average spike rates
    % calculate average spike rate of cell during the display of each of 16 orientations of a moving bar
    trial_number=0;
    orientation=[0:22.5:337.5];
    total_spikes=zeros(16,5); %was (16,10)
    trial_period=zeros(16,5); %was (16,10)
    
    spike_rate=zeros(16,5); %was (16,10)
    unave_true_spike_rate = zeros(16,1);
    true_spike_rate=zeros(16,1); % this was here first!
    norm_true_spike_rate=zeros(16,1); %
    
    x=zeros(16,1);
    y=zeros(16,1);
    
    spikeSumx=0;  %to be used in adding up the vectors for different direction-selectivities
    spikeSumy=0;
    norm_spikeSumx = 0;
    norm_spikeSumy = 0;
    
    error = zeros(16,1);
    unave_error = zeros(16,1);
    norm_error = zeros(16,1);
    
    for i=1:16;
        for j=1:length(stim);
            if stim(j)==orientation(i);
                trial_number = trial_number+1;
                for l=1:length(spikeTimes);
                    if spikeTimes(l)>=start_times(j) && spikeTimes(l)<=End_times(j);
                        total_spikes(i,trial_number)=total_spikes(i,trial_number)+1;
                    end
                end
                trial_period(i,trial_number)=(trial_period(i,trial_number)...
                    +End_times(j)-start_times(j))/20000;
                spike_rate(i,trial_number)=total_spikes(i,trial_number)/trial_period(i,trial_number);
                % spikes per second for a given trial for a given
                % orientation
                
            end
            
        end
        % Here we account for the real number of stimulus display epochs during
        % the recording. If there were some display epochs missed we will still
        % calculate a true average.
        for l = 1:trial_number;
            %               l = 1:2
            %             l = 3:trial_number;  %%%%%THIS IS FOR THE ANOMALOUS RUN!!!
            %             l = 1:trial_number;   %%%%% THIS IS WHAT YOU NORMALLY USE
            
            unave_true_spike_rate(i)=unave_true_spike_rate(i)+spike_rate(i,l);
            
            % probably total spikes could go here?
            
        end
        %%%%%% FOR NORMAL ANALYSIS
        true_spike_rate(i)=unave_true_spike_rate(i)/trial_number; % ave spike rate for a given orientation
        
        
        %%%%% FOR WRONG ANALYSIS
        %            true_spike_rate(i)=unave_true_spike_rate(i)/2; % average spike rate for a given orientation
        %           true_spike_rate(i)=unave_true_spike_rate(i)/3; % average spike rate for a given orientation
        
        
        %create matrix recording # of display epochs during an experiment for
        %each orientation
        x(i)=x(i)+trial_number;
        
        %calculate  std error (on the mean!!)
        unave_error(i)=std(spike_rate(i,1:trial_number)');
        error(i)=unave_error(i)/sqrt(trial_number);
        
        trial_number=0;
    end
%     
%     calculate the average spike rate of cell during grey screens
%     e=0;
%     t=0;
%     u=0;
%     n=0;
%     
%     for j=1:(length(stim)-1);
%         t=t+(start(j+1)-E(j));
%         n=n+1;
%         for l=1:length(spikeTimes);
%             if spikeTimes(l)<start(j+1) && spikeTimes(l)>E(j);
%                 e=e+1;
%             end
%         end
%     end
%     
%     u=e/(t/20000);
    
    
    
    %% plot polar selectivity
    max_rate=max(true_spike_rate); %%%%%% I have to change this because it
    %is no longer the average -  change it to the normalized spike
    %rate
    if printfigures == 1
        figure
        %         subplot ('position',[.25 .5 .45 .45])
        polar(0,max_rate,'-k')
        title('Spike rate (spikes/sec) at each orientation');
        hold on;
        h = polar(orientation(:)*(2*pi/360), true_spike_rate(:), '-k');
        set(h,'LineWidth',2.5);
        h = polar([orientation(16)*(2*pi/360), orientation(1)*(2*pi/360)],...
            [true_spike_rate(16), true_spike_rate(1)], '-k');
        set(h,'LineWidth',2.5);
        hold off
    end
    
    
    %%
    
    for j = 1: 16;
        
        %         if printfigures == 1
        %             polar(orientation(j)*(2*pi/360), true_spike_rate(j), '-.or')
        %             hold on;
        %         end
        spikeSumx=spikeSumx+true_spike_rate(j)*(cos(orientation(j)*(2*pi/360)));
        spikeSumy=spikeSumy+true_spike_rate(j)*(sin(orientation(j)*(2*pi/360)));
        
        if printfigures == 1
            hold off;
        end
    end
    
    %% Getting the normalized spike
    %For each orientation, divide by the sum of the rates of all
    % orientations, and divide errorbars by the same(?)
    
    sum_of_rates = 0;
    for i = 1:16;
        sum_of_rates = sum_of_rates + unave_true_spike_rate(i);
    end
    
    
    for i = 1:16;
        norm_true_spike_rate(i) = unave_true_spike_rate(i) / sum_of_rates;
        norm_error(i) = unave_error(i) / sum_of_rates;
    end
    
    % Then I want normalized x and y vectors to get the coordinates for DS
    % vectors
    
    for j = 1: 16;
        norm_spikeSumx = norm_spikeSumx + norm_true_spike_rate(j)*(cos(orientation(j)*(2*pi/360)));
        norm_spikeSumy = norm_spikeSumy + norm_true_spike_rate(j)*(sin(orientation(j)*(2*pi/360)));
    end
    
%     % This should make a polar plot like before, but with the normalized
%     % spike rate?
%     norm_max_rate=max(norm_true_spike_rate); %%%%%% I have to change this because it
%     %is no longer the average -  change it to the normalized spike
%     %rate
%     if printfigures == 1
%         figure
%         subplot ('position',[.25 .5 .45 .45])
%         polar(0,norm_max_rate,'-k')
%         title('Spike rate (spikes/sec) at each orientation');
%         hold on;
%     end
%     for j = 1 : 16;
%         if printfigures == 1
%             
%             polar(orientation(j)*2*pi/360, norm_true_spike_rate(j), '-.or')
%             title('Normalized Spike Rates for a Given Neuron')
%             hold on
%         end
%     end
%     if printfigures == 1
%         hold off
%     end
    
    
    %% plot neurons tuning curve %% this was commented out
%     if printfigures == 1
%         subplot('position',[0.12 0.1 0.75 0.3])
%     end
%     max_rate=round(max_rate);
%     if printfigures == 1
%         %         errorbar(orientation,true_spike_rate,error,'or');
%         errorbar(orientation,norm_true_spike_rate,norm_error,'or');
%         hold on;
%         %         axis([0,360,0,max_rate+2]);
%         %         y=[0:max_rate+2];
%         
%         axis([0,360,0,1]);
%         y=[0:1];
%         
%         x=[0:30:360];
%         grid off; legend off;
%         set(gca,'XTick',x,'FontSize',7);
%         xlabel('Sq. Wave Orientation (Deg)', 'FontSize', 9);
%         set(gca,'YTick',y,'FontSize',7);
%         %         ylabel('Spike Rate (Hz)','FontSize', 9);
%         ylabel(' Normalized Spike Rate','FontSize', 9);
%         hold on;
%     end
%     
%     
%     peakratio(DSCellCounter)=ratio(q);%the ratio of the two peaks we're interested in
%     
%     if spikeSumx<0
%         AllCells(q,1)=atan(norm_spikeSumy/norm_spikeSumx)+pi;
%         %the theta value for the total preferred direction for each cell
%         
%     else
%         AllCells(q,1)=atan(norm_spikeSumy/norm_spikeSumx);  %(necessary to get correct sign direction)
%         
%     end
    
    
    % Theta value for total preferred direction for each cell
    AllCells(q,1) = atan2(norm_spikeSumy,norm_spikeSumx);
    
    
    %Convert this theta to degrees temporarily
    temp_theta_degrees = [];
    
    if AllCells(q,1) < 0
        temp_theta_degrees = (AllCells(q,1) + 2*pi) * 180/pi;
    else
        temp_theta_degrees = AllCells(q,1) * 180/pi;
    end
    
    % R value for the total preferred direction for each cell
    AllCells(q,2)=((norm_spikeSumx^2)+(norm_spikeSumy^2))^(.5);
    
    
    % % if angle < 0
    % %   angle += 2 pi
    % % angle = angle *180/pi
    % % (only angle += 2 pi is contingent of angle < 0)
    % %
    % % to convert from -pi to pi, you add 2pi to all the angles between -pi and 0, i.e. all the angles < 0
    % % then you have 0 to 2pi, and you can just convert to 0 to 360 by multiplying by 180/pi
    
    
    
    
    %what is all this stuff?
    %     AllCells(q,2)=((spikeSumx^2)+(spikeSumy^2))^(.5); %non-normalized
    
    %% New Determination of preferred orientation!
    difference_vector = [];
    for i = 1 : length(orientation);
        difference_vector(i) = abs(orientation(i) - temp_theta_degrees);
    end
    
    min_diff_vect = min(difference_vector);
    index_of_min = find(difference_vector <= min_diff_vect);
    pref_O = orientation(index_of_min); % the pref_O is based off the Elstrott paper
    
    index_of_opp = mod((index_of_min+7), 16) + 1;
    
    
    %     %Determine preferred orientation, orthogonal orientation and opposite
    %     %orientation
    %     [max_rate,i]=max(true_spike_rate);
    %     pref_O=orientation(i);  % T
    
    
    orth_O=pref_O+90;
    if orth_O > 337.5;
        orth_O = orth_O-360;
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
    
    %Calculate orientation selectivity and direction selectivity
    %(basically just the percent error between 2 values for each one)
    OSI = (pref_Orate - orth_Orate)/(pref_Orate + orth_Orate);
    OSI = round(OSI/0.01)*.01;
    DSI(q) = (pref_Orate - opp_Orate)/(pref_Orate + opp_Orate);
    %     DSI = round(DSI/0.01)*.01;
    
    
    
    %% Tuning Curve Fitting
%     
%     global Rmax mu
%     global mu
%     
%     Rmax = pref_Orate*50;
%     mu = pref_O*pi/180; let mu be the direction closest to the vector sum
%     
%     mu = temp_theta_degrees*pi/180; % let mu be the vector sum    
%     
%     f = fittype('vonMises(k,x)'); % function in /Sherlab
%     f = fittype('vonMises(Rmax,k,x)'); % function in /Sherlab
    f = fittype('vonMises(Rmax,mu,k,x)'); % function in /Sherlab
    
    
    options = fitoptions('Method','NonlinearLeastSquares');
    options.StartPoint = [pref_Orate*50  1 temp_theta_degrees*pi/180];
    options.Lower = [0 0 0];
    
    f = setoptions(f, options);
    
    [fo, gof] = fit((orientation.*pi/180)',true_spike_rate.*50,f); %,'Weights', (1/(error.^2))
    
    %     [f, MSGID] = lastwarn();
    %     warning('off', MSGID)
    
    goodness(:,q) = struct2cell(gof);
    
    MyCoeffs = coeffvalues(fo);
    
    Rmax = MyCoeffs(1);
    mu = MyCoeffs(3);
    MyCoeffs = MyCoeffs(2);
    
%     fprintf('Rmax= %f, mu = %f, k = %f \n', Rmax, mu, MyCoeffs)
    
    %     test = temp_theta_degrees - mu*180/pi
    
    
    chi_square = 0;
    binomial_error = [];
    predicted_spike_rate = [];
    norm_predicted_spike_rate = [];
    
    ci = [];
    ci = confint(fo);
    
    
    
    for i = 1:16
        predicted_spike_rate(i) = (Rmax*exp(MyCoeffs*cos((orientation(i)*(pi/180))-mu)))/exp(MyCoeffs);
    end
    
    
    % for binomial error, SEM = sqrt(p*(1-p)), estimate p as normalized predicted spike rate,
    %
    
    for i = 1: 16
        norm_predicted_spike_rate(i) = predicted_spike_rate(i) / sum(predicted_spike_rate);
    end
    
    for i = 1:16
        binomial_error(i) = sqrt(norm_predicted_spike_rate(i)*(1-norm_predicted_spike_rate(i))/5);
    end
    
    for i = 1:16
        chi_square = chi_square + ...
            ((predicted_spike_rate(i) - true_spike_rate(i)*50)^2)/predicted_spike_rate(i);
    end
    
    for i = 1:16
        chi_square_binom = chi_square + ...
            ((predicted_spike_rate(i) - true_spike_rate(i))^2)/(binomial_error(i)^2);
    end
    
    %      for i = 1:16
    %         chi_square = chi_square + ...
    %             ((norm_predicted_spike_rate(i) - norm_true_spike_rate(i))^2)/(norm_predicted_spike_rate(i));
    %      end
    
    reduced_chi_square(q) = chi_square/15;
%      reduced_chi_square(q) = chi_square/13;

    
%     reduced_chi_square_binom(q) = chi_square_binom/15;
    
    %      test= 0;
    %     test = sum((predicted_spike_rate' - true_spike_rate).^2);
    %     test
    
    
    fwhh(q) = 2*acos(log((1/2)*exp(MyCoeffs)+(1/2)*exp(-MyCoeffs))/MyCoeffs)  ; % in radians
    
    lowerconf = 2*acos(log((1/2)*exp(ci(1))+(1/2)*exp(-ci(1)))/ci(1));
    upperconf = 2*acos(log((1/2)*exp(ci(2))+(1/2)*exp(-ci(2)))/ci(2));
    
    size_confint(q) = lowerconf-upperconf;
    
    fwhh_fractional_error(q) = size_confint(q)/fwhh(q); % smaller numbers are good
    
    if printfigures == 5
        
        figure %von Mises fits
        
        h = polar(orientation(:)*(pi/180), true_spike_rate(:)*50, 'ok'); % data points
        set(h, 'MarkerFaceColor', 'r');
        
        hold on
        
        h = polar(0:0.01:2*pi,(Rmax*exp(MyCoeffs*cos((0:0.01:2*pi)-mu)))/exp(MyCoeffs), '-'); % fitted curve
        set(h,'linewidth',4)
        
        
        h = polar([AllCells(q,1), AllCells(q,1)], [0, 100],'k-'); % DS vector arbitrary length
        set(h, 'linewidth', 2)
        
        h = polar(pref_O*pi/180, pref_Orate*50, 'ok'); % determined maximum point
        set(h, 'MarkerFaceColor', 'y');
        
        title({strcat('ID-', num2str(NeuronID), ', X^2=', num2str(reduced_chi_square(q))),...
            strcat('Pref #Spikes=', num2str(pref_Orate*50), ', CI/FWHH=', num2str(fwhh_fractional_error(q)))},...
            'FontSize', 25);
        
    end
    %     error
    %
    
    MaxSpikes(q) = pref_Orate*50;
    
    
    %%
    
    %make a histogram of the ratios of the two biggest peaks, in the
    %direction which elicited the best response for each particular cell
    %
    %     %     strongestresponsedirection=(AllCells(q,1)*(360/(2*pi)));  %what
    %     it used to be!!
    %     strongestresponsedirection = temp_theta_degrees; % what i thought
    %     would work?
    %
    %     thebeststimulus=round(strongestresponsedirection/22.5); %to make an integer between 1 and 16
    
    
    %
    %     if thebeststimulus<0
    %         thebeststimulus=17-abs(thebeststimulus);
    %     else
    %         thebeststimulus=thebeststimulus+1;
    %
    %     end
    %
    
    
    Y=fft(combinedspikes{index_of_min});  %% all these index_of_mins were thebesttimulus
    n=length(Y);
    if n>1;
        Y(1)=[];
    end
    power = abs(Y(1:floor(n/2))).^2;
    nyquist = 1/2;
    freq = (1:n/2)/(n/2)*nyquist*(length(segments)/10);
    %freq = (1:n/2)/(n/2)*nyquist*(segments/17); %the segments/17 is to transform the units
    %from "spikes per segment" to hertz, since the spiketimes go from -3 seconds to 14 seconds
    %(but the length of each trial actually might not be exactly 17 seconds.....needs to be addressed)
    
    
    %% Calculation of DSI error
    % I did propagation of errors on the DSI by assuming a gaussian and poisson
    % distribution of the total number of spikes (summed)
    
    %     p_spikes = sum(combinedspikes{index_of_min}.*5);
    %     n_spikes = sum(combinedspikes{index_of_opp}.*5);
    
    %     p_spikes = unave_true_spike_rate(index_of_min)*10;
    %     n_spikes = unave_true_spike_rate(index_of_opp)*10;
    
    p_spikes = true_spike_rate(index_of_min)*50;
    n_spikes = true_spike_rate(index_of_opp)*50;
    
    DSI_error(q) = (2/(p_spikes+n_spikes)^2)*(n_spikes*p_spikes*(n_spikes+p_spikes))^(1/2);
    
    
    %%
    if printfigures == 2
        figure
        subplot(2,2,1)
        plot(segments, combinedspikes{index_of_min}) % for vector segments
        % plot(combinedspikes{i})
        title(Direction(index_of_min),'FontSize',9,'FontWeight','bold');
        xlabel('combined spikes')
        ylabel('')
        hold on
        
        subplot(2,2,3)
        plot(freq,power);
        xlabel('Frequency (Hz)'); ylabel('power');
        
        subplot(2,2,2)        %         subplot ('position',[.25 .5 .45 .45])
        polar(0,max_rate,'-k')
        title('Spike rate (spikes/sec) at each orientation');
        hold on;
        h = polar(orientation(:)*(2*pi/360), true_spike_rate(:), '-k');
        set(h,'LineWidth',2.5);
        h = polar([orientation(16)*(2*pi/360), orientation(1)*(2*pi/360)],...
            [true_spike_rate(16), true_spike_rate(1)], '-k');
        set(h,'LineWidth',2.5);
        hold off
        
        subplot(2,2,4)
        
        x_fake=[0 1 0 -1];
        y_fake=[1 0 -1 0];
        
        h_fake=compass(x_fake,y_fake, '--');
        hold on;
        
        %     for i = 1 : length(AllCells)
        %     if AllCells(i,2) >0.4 %normalized magnitude vector
        %         if DSI(i) > 0.5 && ratio(i) > 1 % Direction selectivity index and f2/f1
        ch = compass((AllCells(q,2)*cos(AllCells(q,1))), (AllCells(q,2)*sin(AllCells(q,1))), '-k');
        xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
        set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
        % set(ch,'linewidth',4);
        set(ch,'LineWidth',2)
        %polar2(AllCells(i,1), AllCells(i,2),[0 1], '-k'); % this is in the Matlab folder - wtf polar
        hold on
        title(strcat('Direction Selective Vector for neuron-', num2str(NeuronID)));
        %         end
        %     end
        set(h_fake,'Visible','off')
        
        
    end
    
    if printfigures == 3
        figure
        subplot(1,2,1)
        polar(0,max_rate,'-k')
        title(strcat('Spike rate (spikes/sec) at each orientation for neuron-', num2str(NeuronID)));
        hold on;
        h = polar(orientation(:)*(2*pi/360), true_spike_rate(:), '-k');
        set(h,'LineWidth',2.5);
        h = polar([orientation(16)*(2*pi/360), orientation(1)*(2*pi/360)],...
            [true_spike_rate(16), true_spike_rate(1)], '-k');
        set(h,'LineWidth',2.5);
        hold off
        
        subplot(1,2,2)
        
        x_fake=[0 1 0 -1];
        y_fake=[1 0 -1 0];
        
        h_fake=compass(x_fake,y_fake, '--');
        hold on;
        ch = compass((AllCells(q,2)*cos(AllCells(q,1))), (AllCells(q,2)*sin(AllCells(q,1))), '-k');
        xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
        set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
        % set(ch,'linewidth',4);
        set(ch,'LineWidth',2)
        %          polar2(AllCells(i,1), AllCells(i,2),[0 1], '-k'); % this is in the Matlab folder
        hold on
        title(strcat('Direction Selective Vector for neuron-', num2str(NeuronID)));
        %         end
        %     end
        set(h_fake,'Visible','off')
        
    end
    
    
    %         h = polar(orientation.*pi/180,(Rmax*exp(MyCoeffs*cos((orientation.*pi/180)-mu)))/exp(MyCoeffs));
    %         set(h,'linewidth',4);
    %         hold on
    %         h = polar([orientation(16)*pi/180, orientation(1)*pi/180], ...
    %             [(Rmax*exp(MyCoeffs*cos((orientation(16)*pi/180)-mu)))/exp(MyCoeffs),...
    %             (Rmax*exp(MyCoeffs*cos((orientation(1).*pi/180)-mu)))/exp(MyCoeffs)]);
    %         set(h,'linewidth',4);
    
    %         title(strcat('von Mises fit for neuron-', num2str(NeuronID)));
    
    
    
    %     peaknumberone=power(13:19); %finds the peak in the range where we'd expect the on OR off peak
    %     peaknumbertwo=power(23:29);  %finds the peak in the range where we'd expect the on/off peak
    peaknumberone=power(7:12); %finds the peak in the range where we'd expect the on OR off peak
    peaknumbertwo=power(17:22);  %finds the peak in the range where we'd expect the on/off peak
    
    %[maximums, locs]=findpeaks(power,'sortstr', 'descend','minpeakdistance',minpeakdist);
    [~,I]=max(peaknumberone); %finds the peak in the range where we'd expect the on/off peak
    [~,H]=max(peaknumbertwo); %finds the peak in the range where we'd expect the on OR off peak
    
    
    
    %     ratio(q)=(power(12+I))/(power(22+H));  %the ratio of the two peaks we're interested in
    ratio(q)=(power(16+H))/(power(6+I));  %the ratio of the two peaks we're interested in
    
    noise(q) = mean(horzcat(power(2:6), power(13:16), power(23:26), power(33:36), power(43:46)));
    zerofreq(q) = power(1); % how does this compare to noise?  it should correspond to ambient spikerate?
    
    % noise ratio for either f1 or f2 (depending on which is dominant)
    if (power(16+H)) > (power(6+I));
        noiseratio(q) = (power(16+H))/noise(q);
    elseif (power(16+H)) <= (power(6+I))
        noiseratio(q) = (power(6+I))/noise(q);
    end
    
    
    %     if ((spikeSumx^2)+(spikeSumy^2))^(.5) >40;
    if DSI(q) >0.6
        
        %in case we only want the ratios
        %of the peaks of the highly DS Cells
        
        DSCellCounter = DSCellCounter+1;
        DSCells(DSCellCounter,1)=AllCells(q,1);
        DSCells(DSCellCounter,2)=AllCells(q,2);
        peakratio(DSCellCounter)=ratio(q);
        
        
    end
    
    
    
%     if DSI >= 0.5;
%         f = fittype('gauss1');
%         options = fitoptions('gauss1');
%         [gfit,gof] = fit(orientation',true_spike_rate,f);
%         plot(gfit,'r');
%         legend off;
%     else
%         if DSI < 0.5;
%             f1 = fittype('gauss2');
%             options = fitoptions('gauss2');
%             [gfit,gof] = fit(orientation',true_spike_rate,f1);
%             plot(gfit,'r');
%             legend off;
%         end
%     end
%     hold on;
%     
%     show OSI & DSI on plot
%     OSIstr=num2str(OSI);
%     str4=strcat('OSI:   ', OSIstr);
%     text(360,2,str4,'FontSize',8);
%     DSIstr=num2str(DSI);
%     str5=strcat('DSI:   ', DSIstr);
%     text(360,1,str5,'FontSize',8)
%     ,
%     %plot the cells average spike waveform
%     e=neuronFile.getElectrode(NeuronID);
%     n=length(spikeTimes);
%     x=0:50:100;
%     y=-200:100:100;
%     subplot('position',[0.7 0.75 0.25 0.2])
%     m=0;
%     S=zeros(n,70);
%     for j=1:n
%         t=spikeTimes(j);
%         d1=rawFile.getData(e,t-20,70);
%         S(j,:)=d1;
%     end
%     a=mean(S);
%     a=a/1.8;
%     plot(a)
%     xlabel('samples(20KHz)','FontSize',9)
%     ylabel('microvolts','FontSize',9)
%     set(gca,'YTick',y,'FontSize',7);
%     set(gca,'XTick',x,'FontSize',7);
%     
%     plot an amplitude histogram
%     First, extract spike times and amplitudes
%     
%     amptime=spikeFile.getSpikeTimesAmplitudes(e);
%     
%     
%     
%     %spike times are the first half of the list
%     times=amptime(1:(length(amptime)/2));
%     
%     %Spike amps in ADC counts are the second half of list
%     amps=amptime((length(amptime)/2)+1:length(amptime));
%     
%     %Verify that the # of spike times in this list is same as found before
%     %(lines 47 & 48)
%     neuronAmps=zeros(length(spikeTimes),1);
%     for j=1:length(spikeTimes)
%         for k=1:length(times)
%             if spikeTimes(j)==times(k)
%                 neuronAmps(j)=amps(k);
%             end
%         end
%     end
%     %Convert amps from ADC counts to microvolts
%     neuronAmps=abs(neuronAmps);
%     neuronAmps=neuronAmps/1.8;
%     
%     %plot an amplitude histogram
%     subplot('position',[0.7 0.4 0.25 0.2]);
%     histx=[1:1:max(neuronAmps)];
%     [N,x]=hist(neuronAmps,histx);
%     bar(N);
%     hold on
%     
%     % Fit a gaussian
%     f = fittype('gauss1');
%     options = fitoptions('gauss1');
%     options.Lower = [0 -Inf 0 0 -Inf 0]; % sets the boundaries
%     [gfit,gof]=fit(x',N',f);
%     plot(gfit);
%     legend off;
%     rsqrd2=gof.adjrsquare;
%     rsqrd2=round(rsqrd2/.001)*.001;
%     rsqrd2str=num2str(rsqrd2);
%     rsqrd2str=['R^2: ' rsqrd2str];
%     text(5,max(N)-20,rsqrd2str,'FontSize',7);
%     xlabel('Spike amplitude (microvolts)','FontSize',8);
%     ylabel('Counts','FontSize',8);
%     axis ([0,max(neuronAmps),0,max(N)+.05*max(N)]);
%     hold off
%     
%     %make a scatter plot of spike amplitudes against experiment duration
%     subplot('position', [0.7 0.1 0.25 0.2]);
%     U=zeros(1,length(neuronAmps));
%     V=zeros(1,length(neuronAmps));
%     for j=1:length(neuronAmps)
%         U(j)=neuronAmps(j);
%         V(j)=spikeTimes(j)/20000;
%     end
%     
%     plot(V,U,'.')
%     x=0:1000:length(spikeTimes)/20000;
%     y=0:50:200;
%     xlabel('Time (sec)','FontSize',9);
%     ylabel('SpikeAmp (microvolts)','FontSize', 9);
%     set(gca,'YTick',y,'FontSize',7);
%     set(gca,'XTick',x,'FontSize',7);
%     
%     NeuronID=num2str(NeuronID);
%     
%     if OSI >=.5 && DSI < 0.5;
%         cd(B);
%         saveas(gcf,NeuronID);
%         OSCell=OSCell+1;
%     else
%         if OSI >=.5 && DSI >= 0.5;
%             cd(C);
%             saveas(gcf,NeuronID);
%             DSCell= DSCell+1;
%         else
%             if OSI < 0.5;
%                 cd(D);
%                 saveas(gcf,NeuronID);
%                 NonOSCell=NonOSCell+1;
%             end
%         end
%     end
%     
%     close all
%     keyboard
end

%% START OF WHOLE RETINA ANALYSIS
if printfigures == 4
    
%     hold on;
%     
%     figure
%     for i = 1: length(q);
%         if ratio(i) >1
%             plot(DSI(i), AllCells(i,2))
%             xlabel('DSI'); ylabel('Magnitude of DS vector')
%         end
%     end
%     
%     
%     This stuff is useful when dealing with LOTS of cells in one data set:
%     the first figure makes a polar plot of the most highly-DS cells
%     the second figure makes a histogram of the ratios of the amplitudes of
%     the highest and second-highes peaks in the fourier transform.
%     
%     figure
%     grid OFF
%     
%     for i = 1 : length(AllCells)
%         %     if AllCells(i,2) >0.4 %normalized magnitude vector
%         if DSI(i) > 0.5 % Direction selectivity index
%             %     polar(AllCells(i,1), AllCells(i,2), 'ko');
%             
%             polar2(AllCells(i,1), AllCells(i,2),[0 1], 'ko'); % this is in the Matlab folder - wtf polar
%             hold on
%         end
%     end
    %% Compass or Polar Plot for the whole Retina
    figure
    x_fake=[0 1 0 -1];
    y_fake=[1 0 -1 0];
    
    h_fake=compass(x_fake,y_fake, '--');
    hold on;
    
    for i = 1 : length(AllCells)
        %     if AllCells(i,2) >0.4 %normalized magnitude vector
        if DSI(i) > 0.5 && ratio(i) > 1 && noiseratio(i) > 5% Direction selectivity index and f2/f1
            ch = compass((AllCells(i,2)*cos(AllCells(i,1))), (AllCells(i,2)*sin(AllCells(i,1))), '-k');
            xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
            set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
            % set(ch,'linewidth',4);
            
            %polar2(AllCells(i,1), AllCells(i,2),[0 1], '-k'); % this is in the Matlab folder - wtf polar
            hold on
        end
    end
    set(h_fake,'Visible','off')
    
    %% DSI histogram
    figure
    
    rgreat1ind = find(ratio(:) >1 & noiseratio(:) > 5);
    hist(DSI(rgreat1ind), (-0.975:0.05:0.975))
    
    figure
    hist(AllCells((rgreat1ind),2), (0.0125:0.025:0.9875))
    
    xlabel('Magnitude of DS Vector', 'FontSize', 15);
    ylabel('# of Healthy ON/OFF Neurons', 'FontSize', 15);
    
    
    
%     %
%     figure
%     subplot('position', [.05, .05, .8, .8])
%     polar(0,60,'-k')
%     hold on;
%     for j=1:DSCellCounter;
%         
%         polar(DSCells(j,1), DSCells(j,2),'.r')
%         
%         hold on;
%         
%     end
%     
%     
%     hold off;
%     figure
%     
%     hist(ratio(:),4000)
%     axis([0,5,0,25])

    
    
end
DScellsindex = [];
DScellsindex = find(ratio(:) >1 & DSI(:) > 0.5 & noiseratio(:) >5);
% DScells index was originally  1, 0.5 , 5

% keyboard
if printfigures == 0;
    %     filename = strcat(A, '_', datestr(now, 30)); % maybe you don't want this datestr for now
    filename = A;
    save(filename, 'idList', 'AllCells', 'DSI','DSI_error', 'ratio',...
        'noiseratio', 'DScellsindex', 'noise', 'zerofreq', 'fwhh',...
        'reduced_chi_square', 'MaxSpikes');
end
%   keyboard
