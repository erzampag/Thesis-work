function [R, theta] = calculate_R_theta(norm_true_spike_rate, orientation)

%%% Calculation of R and theta
    %
    %     if spikeSumx<0
    %         AllCells(q,1)=atan(norm_spikeSumy/norm_spikeSumx)+pi;
    %         %the theta value for the total preferred direction for each cell
    %
    %     else
    %         AllCells(q,1)=atan(norm_spikeSumy/norm_spikeSumx);  %(necessary to get correct sign direction)
    %
    %     end
    
    %     spikeSumx=0;  %to be used in adding up the vectors for different direction-selectivities
    %     spikeSumy=0;
    norm_spikeSumx = 0;
    norm_spikeSumy = 0;
    
    for j = 1:length(orientation);
        norm_spikeSumx = norm_spikeSumx + norm_true_spike_rate(j)*(cos(orientation(j)*(pi/180)));
        norm_spikeSumy = norm_spikeSumy + norm_true_spike_rate(j)*(sin(orientation(j)*(pi/180)));
    end
    % % %     keyboard
    
    
    % R value for the total preferred direction for each cell
    R = ((norm_spikeSumx^2)+(norm_spikeSumy^2))^(.5);
    
    % Theta value for total preferred direction for each cell
    theta = atan2(norm_spikeSumy,norm_spikeSumx); %%%% -pi <= ATAN2(Y,X) <= pi.
    
end