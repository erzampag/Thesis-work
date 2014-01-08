%% UNUSED CODE IN drifting_squarewave_shorter_stimulus2
%% normalized polar plot
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


%% fits to gaussians, cell average spike waveforms, amplitude histograms    
    
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
%     
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

%% Old Whole Retina Plots
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

 
