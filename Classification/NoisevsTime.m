% Want to look at noise at various points in the raw data


%% import noise data from textfiles labled "noise0000.txt", where number is "time to skip"
% directory_path = ('/Users/erinzampaglione/Desktop/2011-02-23-raw');
% 
directory_path = ('/Users/erinzampaglione/Desktop/2011-06-28-raw');

files = dir(directory_path);


% figure(111)
all_electrode_noise = zeros(513, 9);
time_at_noise_measurement = zeros(9,1);
conformity_ratio = zeros(9,1);

counter = 0;

for i = 1 : length(files)
    if strcmp(files(i).name(1:1), 'n')
        
        counter = counter + 1;
        
        file_path = strcat(directory_path, '/', files(i).name);
        file = fopen(file_path);
        
        text1 = textscan(file, '%f');
        fclose(file);
        
        all_electrode_noise(:,counter) = text1{1};
        
        time_at_noise_measurement(counter) = str2double(files(i).name(6:9));
        
        
        % calculate conformity ratio
        conformity_ratio(counter) = mean(all_electrode_noise(2:end,counter)) / std(all_electrode_noise(2:end,counter));
        
    end
end
% keyboard
%% plot histograms of all the electrodes RMS

for i = 1 : length(time_at_noise_measurement)
    
    figure

    hist(all_electrode_noise(2:end,i), 2.5:5:55.5) % bins chosen for 2011-02-23
    %  hist(all_electrode_noise, 2.5:5:55.5) % bins chosen for 2011-02-23
    
    title(['Distribution of noise for electrodes at time ' num2str(time_at_noise_measurement(i))])
    
    
    text(10,350, ['mean = ' num2str(mean(all_electrode_noise(2:end,i))) ', std = ' num2str(std(all_electrode_noise(2:end,i)))...
        ', conformity ratio = ' num2str(conformity_ratio(i))])
    
    ylim([0 360])

    saveas(gcf, ['NoiseDistAt' num2str(time_at_noise_measurement(i))], 'eps')
end


% keyboard

%% plot
figure

plot(time_at_noise_measurement, mean(all_electrode_noise,1), '-');
hold on
errorbar(time_at_noise_measurement, mean(all_electrode_noise,1), std(all_electrode_noise) / sqrt(length(all_electrode_noise)))

title(['Noise as a Function of Time for Retina ' directory_path(32:41)])

ylim([12 30]); xlim([-50 2500])


saveas(gcf, ['NoiseofTfor' directory_path(32:41)], 'eps')


