function [Fone_Fzero, Ftwo_Fzero, OS_x1, OS_y1, OS_x2, OS_y2] = ...
    determine_prefO_by_FFT(true_spike_rate, norm_true_spike_rate,  AllCells, pref_Orate, pref_O,...
    orientation, printfigures, q, NeuronID, Max_Spikes_absolute)


    temp_theta_degrees = [];
    
    if AllCells(q,1) < 0
        temp_theta_degrees = (AllCells(q,1) + 2*pi) * 180/pi;
    else
        temp_theta_degrees = AllCells(q,1) * 180/pi;
    end
    
%%% calculation of DS and OS by FFT
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
    
    Fone_Fzero = amp_tuning(2)/amp_tuning(1);
    Ftwo_Fzero = amp_tuning(3)/amp_tuning(1);

    
    % Only use this if statment if you want to fit the tuning curve to only
    % OS cells as defined by the statment
    %     if Ftwo_Fzero(q) > 0.03 && Ftwo_Fzero(q) > Fone_Fzero(q) && true_spike_rate(max_index)*50 > 100
    %     if Ftwo_Fzero(q) > sqrt(0.03) && Ftwo_Fzero(q) > Fone_Fzero(q) && Max_Spikes_absolute > 100
    
%     OS_counter = OS_counter +1;
    
    if Ftwo_Fzero > Fone_Fzero
        % fits a double Von Mises curve and gives angle of largest Rmax
        mu = fit_tuning_curve2(pref_Orate, pref_O, temp_theta_degrees, true_spike_rate, orientation, NeuronID, q,...
            printfigures, AllCells, Ftwo_Fzero, Fone_Fzero);
        
        mu = mod(mu,2*pi);
        
        % finds stim orientation closest to mu
        pref_OS_dir = orientation(abs(orientation.*pi/180-mu) == min(abs(orientation.*pi/180-mu)));
        
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
        
        OS_x1 = OS_x1temp;
        OS_x2 = OS_x2temp;
        OS_y1 = OS_y1temp;
        OS_y2 = OS_y2temp;
        
    end

    % draw red line
% % %     if printfigures == 7 && Ftwo_Fzero > sqrt(0.03) && Ftwo_Fzero > Fone_Fzero && max(true_spike_rate)*50 > 100
    if printfigures == 7  
        gcf;
        hold on;
        h = polar([pref_OS_dir*pi/180-pi/2;pref_OS_dir*pi/180+pi/2],[Max_Spikes_absolute;Max_Spikes_absolute], 'r-');
        set(h,'linewidth',4);
        
        % draw vector
        subplot(1,2,2);
        
        x_fake=[0 1 0 -1];
        y_fake=[1 0 -1 0];
        
        h_fake=compass(x_fake,y_fake, '--');
        hold on;
        set(h_fake,'Visible','off');

        ch = compass(OS_x1, OS_y1, '-k');
        xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
        set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
        set(ch,'linewidth',4);
        ch = compass(OS_x2, OS_y2, '-k');
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