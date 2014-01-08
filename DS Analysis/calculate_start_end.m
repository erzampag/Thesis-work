function [start_times, End_times] = calculate_start_end(neuronFile)
% create a matrix of start times and end times in sample #s corresponding
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

end