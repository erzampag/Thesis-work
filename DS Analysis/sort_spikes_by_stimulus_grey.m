function [total_spikes, total_spikes_grey, combinedspikes, segments] = sort_spikes_by_stimulus_grey(orientation, stim, spikeTimes,...
    start_times, End_times, printfigures)

% Creation of the raster plot      

% plot spike times as hash marks on figure subplots, and form an array of
    % arrays of arrays of spike times for each orientation, trial.
    %     fig1 = figure(1);

    allspikes={};
    combinedspikes={};
    extractedspikes={};
    total_spikes = zeros(16,5);
    total_spikes_grey = zeros(16,5);


    if printfigures == 1
        figure
    end
   %% allspikes and total_spikes
    for i=1:16;
        
        if printfigures == 1
            
            axes('position', [(mod((i-1),4)/4)+.05, ((4-ceil((i*4)/16))/4)+.05, .15, .15 ]);
            axis([-3,13,0,5]);
            title(orientation(i),'FontSize',9,'FontWeight','bold');
            
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
                            total_spikes(i,trial_count)=total_spikes(i,trial_count)+1;
                            
                        elseif spikeTimes(l) > End_times(j) && spikeTimes(l) < End_times(j)+60000; %correlate with graphs, not start/end (longer than 3s)
                            total_spikes_grey(i,trial_count)=total_spikes_grey(i,trial_count)+1;   
                        end
                        
                    end
                end
                if counter2==0;
                    allspikes{i, trial_count}=[];
                    
                end
            end
        end
        
        
% % %          keyboard
        
        if printfigures == 1
            hold on;
        end
        
    end
    if printfigures == 1
        hold off;
    end
    
    
    segments= .05:.1:9.95;
    %segments = 100;
    minpeakdist=.04*segments;   %to ensure that the peaks aren't right next to each other
    if printfigures == 1
        figure
    end
    
    %% Combined Spikes
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
            title(orientation(i),'FontSize',9,'FontWeight','bold');
            xlabel('combined spikes')
            ylabel('')
            hold on
        end
    end
    if printfigures == 1
        hold off;
    end
    
end
    