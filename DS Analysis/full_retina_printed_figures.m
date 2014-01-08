function full_retina_printed_figures(Fone_Fzero, Ftwo_Fzero, AllCells, idList, DSI, ratio, noiseratio,...
    MaxSpikes, Max_Spikes_absolute, OS_x1, OS_y1, OS_x2, OS_y2)
    %% Histograms and scatter plots of amplitudes (or power = amp^2) of FFT of Tuning Curves for DS (high F1) vs OS (high F2) cells
    
    figure
    hist(Fone_Fzero);
    xlabel('F1/F0'); ylabel('# of neurons');
    figure
    hist(Ftwo_Fzero);
    xlabel('F2/F0'); ylabel('# of neurons');
    figure
    scatter(Fone_Fzero, Ftwo_Fzero);
    xlabel('Fone_Fzero'); ylabel('Ftwo_Fzero');
    figure
    scatter(Fone_Fzero.^2, Ftwo_Fzero.^2);
    xlabel('Fone_Fzero squared'); ylabel('Ftwo_Fzero sqared');
    
    %% Compass or Polar Plot for the whole Retina for DS
    figure
    x_fake=[0 1 0 -1];
    y_fake=[1 0 -1 0];
    
    h_fake=compass(x_fake,y_fake, '--');
    hold on;
    
    for i = 1 : length(AllCells)
        %     if AllCells(i,2) >0.4 %normalized magnitude vector
        if DSI(i) > 0.5 && ratio(i) > 1 && noiseratio(i) > 5 && MaxSpikes(i) > 100 % Direction selectivity index and f2/f1
%         if DSI(i) > 0.5 && noiseratio(i) > 5 && MaxSpikes(i) > 100 % WITHOUT ON/OFF requirement

            ch = compass((AllCells(i,2)*cos(AllCells(i,1))), (AllCells(i,2)*sin(AllCells(i,1))), '-k');
            xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
            set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
            % set(ch,'linewidth',4);
            
            %polar2(AllCells(i,1), AllCells(i,2),[0 1], '-k'); % this is in the Matlab folder - wtf polar
            hold on
            set(ch,'LineWidth',1.5)
            
% % %             if ratio(i) > 1   %%% for ON or OFF vs ON/OFF
% % %             set(ch, 'Color', 'k');
% % %             else
% % %                set(ch, 'Color', 'b');  
% % %             end
% % %             
            
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
    
    %% OS
    figure
    x_fake=[0 1 0 -1];
    y_fake=[1 0 -1 0];
    
    h_fake=compass(x_fake,y_fake, '--');
    hold on;
    
    for i = 1: length(idList)
        if Ftwo_Fzero(i) > sqrt(0.03) && Ftwo_Fzero(i) > Fone_Fzero(i) && Max_Spikes_absolute(i) > 100
%         if Ftwo_Fzero(i) > sqrt(0.03) && Ftwo_Fzero(i) > Fone_Fzero(i) %&& max(true_spike_rate)*50 > 100
        ch = compass(OS_x1(i), OS_y1(i), '-k');
        xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
        set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
        set(ch,'linewidth',4);
        
        ch = compass(OS_x2(i), OS_y2(i), '-k');
        xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
        set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
        set(ch,'linewidth',4);
        hold on

        set(ch,'LineWidth',1.5)
        end
    end
    set(h_fake,'Visible','off')
end
