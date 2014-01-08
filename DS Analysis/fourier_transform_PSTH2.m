function [ratio, noise, zerofreq, noiseratio] = fourier_transform_PSTH2(segments, frames, thebeststimulus, orientation,...
    combinedspikes, this_temporal, this_spatial, printfigures)


Y=fft(combinedspikes{thebeststimulus});  %% all these index_of_mins were thebesttimulus
n=length(Y);
if n>1;
    zero_freq_mode = Y(1);
    Y(1)=[];
end
power = abs(Y(1:floor(n/2))).^2;
nyquist = 1/2;

% %     freq = (1:n/2)/(n/2)*nyquist*(length(segments)/10); % when segments = [0.1:0.1:10]

freq = (1:n/2)/(n/2)*nyquist*(length(segments)/(frames/120));




%%%%% Need to re-define peak 1 and peak 2 to be correct for any temporal period

peak1 = 1/(this_temporal/120);
peak2 = 2*peak1;

closestPeak1 = find(abs(freq - peak1) == min(abs(freq - peak1))); % this is the closest point to peak 1
closestPeak2 = find(abs(freq - peak2) == min(abs(freq - peak2))); % this is the closest point to peak 2

%     range1 =  floor(closestPeak1 - length(freq)/8) : ceil(closestPeak1 + length(freq)/8);
%     range2 =  floor(closestPeak2 - length(freq)/8) : ceil(closestPeak2 + length(freq)/8);


range1 =  floor(closestPeak1 - 1) : ceil(closestPeak1 + 1);
range2 =  floor(closestPeak2 - 1) : ceil(closestPeak2 + 1);

[maxpower1,I]=max(power(range1)); %finds the peak in the range where we'd expect the on/off peak
[maxpower2,H]=max(power(range2)); %finds the peak in the range where we'd expect the on OR off peak

ratio = maxpower2/maxpower1;

% temporary noise stand in:
%%%%% Need to re-define the noise to be the non-harmonic sections of
%%%%% the data
noise_range = horzcat(1 : range1(1)-1,  range1(end)+1 : range2(1) - 1);

noise = mean(power(noise_range));


zerofreq = zero_freq_mode;

% noise ratio for either f1 or f2 (depending on which is dominant)

if maxpower2 > maxpower1; %f2, or ON/OFF
    noiseratio = maxpower2/noise;
else %f1 or ON or OFF
    noiseratio = maxpower1/noise;
end


%%%%% PRINT
if printfigures == 9
    
    figure; plot(freq, power)
    
    xlabel('freq Hz'); ylabel('power');
    title(['sp: ', num2str(this_spatial), ' temp: ', num2str(this_temporal), ' orientation: ', num2str(orientation(thebeststimulus))]);
    hold on;
    
    % want to plot interval for first and second peaks, mark where
    % those two peaks are, and plot the intervals used for the noise to
    % ensure they are not overlapping with harmonic contributions.
    
    plot(freq(range1), ones(length(range1),1)*60, '-r');
    plot(freq(range2), ones(length(range2),1)*60, '-b');
    
    plot(freq(power == maxpower1), power(power == maxpower1), '.r');
    plot(freq(power == maxpower2), power(power == maxpower2), '.r');
    
    plot(freq(noise_range), ones(length(noise_range),1)*60, 'x');
    
    % % % %         index = find(power==max(power));
    % % % %         plot(freq(index), power(index), '.r');
    
end




if printfigures == 1 % graph every FFT (only one is used above)
    

            figure
            for i = 1 : length(orientation)
                
                axes('position', [(mod((i-1),4)/4)+.05, ((4-ceil((i*4)/16))/4)+.05, .15, .15 ]); % a subplot for each orientation
                
                Y=fft(combinedspikes{i});
                n=length(Y);
                if n>1;
                    Y(1)=[];
                end
                power = abs(Y(1:floor(n/2))).^2;
                nyquist = 1/2;
%                 freq = (1:n/2)/(n/2)*nyquist*(length(segments)/10); % when segments = [0.1:0.1:10]
                
                freq = (1:n/2)/(n/2)*nyquist*(length(segments)/(frames/120));

                %             freq = (1:n/2)/(n/2)*nyquist*(segments/10); % when
                %             segments = 100
                %          %the segments/10 is to transform the
                %         units from "spikes per segment" to hertz, since the spiketimes go from 0 to 10 seconds
                
                plot(freq,power);
                title(orientation(i),'FontSize',9,'FontWeight','bold');
                %
                %
                % %loglog(freq,power)   % should i turn this one?!
                xlabel('frequency (Hz)'); %set(gca,'XTick',0:1:5)
                ylabel('amplitude')
                
                % axis([0, 50, 0, 2000]);
                hold on;
                
            end
                        mtit(['Spatial: ', num2str(this_spatial), '  Temporal: ' num2str(this_temporal)]) % after the figure was made

            hold off;

    
end


end


%%% OLD HARDCODED PEAK FINDING

%freq = (1:n/2)/(n/2)*nyquist*(segments/17); %the segments/17 is to transform the units
%from "spikes per segment" to hertz, since the spiketimes go from -3 seconds to 14 seconds
%(but the length of each trial actually might not be exactly 17 seconds.....needs to be addressed)


% % % % % % %     plot(freq,power)


%     peaknumberone=power(13:19); %finds the peak in the range where we'd expect the on OR off peak
%     peaknumbertwo=power(23:29);  %finds the peak in the range where we'd expect the on/off peak
%     peaknumberone=power(7:12); %finds the peak in the range where we'd expect the on OR off peak
%     peaknumbertwo=power(17:22);  %finds the peak in the range where we'd expect the on/off peak

%[maximums, locs]=findpeaks(power,'sortstr', 'descend','minpeakdistance',minpeakdist);
%     [~,I]=max(peaknumberone); %finds the peak in the range where we'd expect the on/off peak
%     [~,H]=max(peaknumbertwo); %finds the peak in the range where we'd expect the on OR off peak

%     ratio(q)=(power(12+I))/(power(22+H));  %the ratio of the two peaks we're interested in
%     ratio=(power(16+H))/(power(6+I));  %the ratio of the two peaks we're interested in

%     noise = mean(horzcat(power(2:6), power(13:16), power(23:26), power(33:36), power(43:46)));

% % % %     zerofreq = power(1); % how does this compare to noise?  it should correspond to ambient spikerate?
%^- would it correspond to offset but that value is removed?  %%%%NOTE
%i do not think that this is the actual zerofrequency - that is taken
%out at the beginning

%     if (power(16+H)) > (power(6+I));
%         noiseratio = (power(16+H))/noise;
%     elseif (power(16+H)) <= (power(6+I))
%         noiseratio = (power(6+I))/noise;
%     end




