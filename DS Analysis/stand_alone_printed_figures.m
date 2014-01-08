function [] = stand_alone_printed_figures(printfigures, combinedspikes, segments,...
    orientation, true_spike_rate, NeuronID, AllCells, index_of_min, q, Fone_Fzero, Ftwo_Fzero, power_tuning, thebeststimulus)

%% FFT
max_rate=max(true_spike_rate); %%%%%% I have to change this because it
%is no longer the average -  change it to the normalized spike
%rate
if printfigures == 1
    figure
    for i=1:16;
        axes('position', [(mod((i-1),4)/4)+.05, ((4-ceil((i*4)/16))/4)+.05, .15, .15 ]);
        
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
        
        plot(freq,power);
        title(orientation(i),'FontSize',9,'FontWeight','bold');
        %
        %
        % %loglog(freq,power)   % should i turn this one?!
        xlabel('frequency (Hz)'); set(gca,'XTick',0:1:5)
        ylabel('amplitude')
        
        % axis([0, 50, 0, 2000]);
        hold on;
        
    end
    hold off;
end

%% Polar Plot
if printfigures == 1
    figure
    %         subplot ('position',[.25 .5 .45 .45])
    polar(0,max_rate,'-k')
    title(strcat('Spike rate (spikes/sec) at each orientation for neuron-', num2str(NeuronID)));
    hold on;
    h = polar(orientation(:)*(2*pi/360), true_spike_rate(:), '-k');
    set(h,'LineWidth',2.5);
    h = polar([orientation(16)*(2*pi/360), orientation(1)*(2*pi/360)],...
        [true_spike_rate(16), true_spike_rate(1)], '-k');
    set(h,'LineWidth',2.5);
    
    h = polar([AllCells(q,1), AllCells(q,1)], [0, 5],'k-'); % DS vector arbitrary length
    set(h, 'linewidth', 2)
    
% %     set(findall(gcf, 'String', '0', '-or', 'String', '30','-or', 'String', '60','-or', 'String', '90', '-or',...
% %         'String', '120', '-or', 'String', '150','-or', 'String', '180','-or', 'String', '210', '-or',...
% %         'String', '240', '-or', 'String', '270','-or', 'String', '300','-or', 'String', '330','-or',...
% %         'String', '2', '-or', 'String', '4','-or', 'String', '6','-or', 'String', '8'), 'String', '')

% TURN ON FOR CLEAN GRAPH
% delete(findall(gcf, 'type', 'text'));    


    hold off
end
%%
if printfigures == 2
    % repeated for purposes of printing figure
     Y=fft(combinedspikes{thebeststimulus});  %% all these index_of_mins were thebeststimulus
    n=length(Y);
    if n>1;
        Y(1)=[];
    end
    power = abs(Y(1:floor(n/2))).^2;
    nyquist = 1/2;
    freq = (1:n/2)/(n/2)*nyquist*(length(segments)/10);
    
    figure
    subplot(2,2,1)
    plot(segments, combinedspikes{thebeststimulus}) % for vector segments
    % plot(combinedspikes{i})
    title(orientation(thebeststimulus),'FontSize',9,'FontWeight','bold');
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
%%
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
%%
% % keyboard
if printfigures == 6
    max_index = find(true_spike_rate == max(true_spike_rate));
    if length(max_index) >1
        max_index = max_index(1);
    end
    
%     if Ftwo_Fzero(q) > 0.03 && Ftwo_Fzero(q) > Fone_Fzero(q) && true_spike_rate(max_index)*50 > 100
    if Ftwo_Fzero(q) > sqrt(0.03) && Ftwo_Fzero(q) > Fone_Fzero(q) && true_spike_rate(max_index)*50 > 100

%     if Ftwo_Fzero(q) > 0.03 && Fone_Fzero(q) > 0.45

        figure
        subplot(2,2,4)
        plot(power_tuning(2:end))
        title('FFT on tuning curve'); xlabel('components'); ylabel('power')
        
        subplot(2,2,2)
        polar(0,max_rate,'-k')
        title(strcat('Spike rate (spikes/sec) at each orientation for neuron-', num2str(NeuronID)));
        hold on;
        h = polar(orientation(:)*(2*pi/360), true_spike_rate(:), '-k');
        set(h,'LineWidth',2.5);
        h = polar([orientation(16)*(2*pi/360), orientation(1)*(2*pi/360)],...
            [true_spike_rate(16), true_spike_rate(1)], '-k');
        set(h,'LineWidth',2.5);
        hold off
        
%         max_index = find(true_spike_rate == max(true_spike_rate));
            Y=fft(combinedspikes{max_index});  %% all these index_of_mins were thebesttimulus
    n=length(Y);
    if n>1;
        Y(1)=[];
    end
    power = abs(Y(1:floor(n/2))).^2;
    nyquist = 1/2;
    freq = (1:n/2)/(n/2)*nyquist*(length(segments)/10);
    
    
    subplot(2,2,1)
    plot(segments, combinedspikes{max_index}) % for vector segments
    % plot(combinedspikes{i})
    title(orientation(max_index),'FontSize',9,'FontWeight','bold');
    xlabel('combined spikes')
    ylabel('')
    hold on
    
    subplot(2,2,3)
    plot(freq,power);
    xlabel('Frequency (Hz)'); ylabel('power');
        
        
    end
end



end