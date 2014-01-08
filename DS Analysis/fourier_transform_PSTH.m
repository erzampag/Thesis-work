function [ratio, noise, zerofreq, noiseratio] = fourier_transform_PSTH(segments, thebeststimulus, combinedspikes, q)


Y=fft(combinedspikes{thebeststimulus});  %% all these index_of_mins were thebesttimulus
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
    
    
    % % % % % % %     plot(freq,power)
    
    
    %     peaknumberone=power(13:19); %finds the peak in the range where we'd expect the on OR off peak
    %     peaknumbertwo=power(23:29);  %finds the peak in the range where we'd expect the on/off peak
    
    peaknumberone=power(7:12); %finds the peak in the range where we'd expect the on OR off peak
    peaknumbertwo=power(17:22);  %finds the peak in the range where we'd expect the on/off peak
    
    %[maximums, locs]=findpeaks(power,'sortstr', 'descend','minpeakdistance',minpeakdist);
    [~,I]=max(peaknumberone); %finds the peak in the range where we'd expect the on/off peak
    [~,H]=max(peaknumbertwo); %finds the peak in the range where we'd expect the on OR off peak
    
    
    
    %     ratio(q)=(power(12+I))/(power(22+H));  %the ratio of the two peaks we're interested in
    ratio=(power(16+H))/(power(6+I));  %the ratio of the two peaks we're interested in
    
    noise = mean(horzcat(power(2:6), power(13:16), power(23:26), power(33:36), power(43:46)));
    zerofreq = power(1); % how does this compare to noise?  it should correspond to ambient spikerate?
    %^- would it correspond to offset but that value is removed?
    
    % noise ratio for either f1 or f2 (depending on which is dominant)
    if (power(16+H)) > (power(6+I));
        noiseratio = (power(16+H))/noise;
    elseif (power(16+H)) <= (power(6+I))
        noiseratio = (power(6+I))/noise;
    end
    
end
    