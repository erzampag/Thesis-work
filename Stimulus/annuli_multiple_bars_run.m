function [new_VBL_timestamps] = annuli_multiple_bars_run(num_textures, TTL_interval, parameters, annulus_positions, annulus_radius)
    
    global win;
    global screen_width;
    global screen_height;
    global handle;
    
% % % sca 
% % % keyboard

    keyIsDown = 0;
    texture_cnt = 0;
    gratings_deviation_red = parameters(1);
    gratings_deviation_green = parameters(2);
    gratings_deviation_blue = parameters(3);
    gratings_temporal_period = parameters(4);
    gratings_spatial_period = parameters(5);
    gratings_angle = parameters(6);
    
    
    % Computing speed in pixels per frame
        gratings_speed = gratings_spatial_period / gratings_temporal_period;


    
    % Compute texture dimensions based on gratings_angle and screen_width, screen_height. Reduces angles to first quadrant if needed.
    [texture_width, texture_height] = texture_dimensions(screen_width, screen_height, gratings_angle);
    
    if annulus_radius ~= 0
        %%%%%%% Annulus Code %%%%%%%%%%%%%%%%%%%
        annulus = 255*(~Circle(annulus_radius)); % points inside annulus are 0, outside are 255
        annulus_texture = repmat(255, screen_height, screen_width); % creates the background
         for a = 1:100 % in Bars and Gratings, it runs multi_annulus_positions_generator.m instead of Daniele's
            annulus_texture(annulus_positions(a,1):annulus_positions(a,2),annulus_positions(a,3):annulus_positions(a,4)) = annulus; % puts the circle in the background
         end
        annulus_mask = ones(screen_height, screen_width, 2) * 128; % first layer, gray as the gratings background
        annulus_mask(:,:,2) = annulus_texture; % second layer, 0 where the circle is, 255 where it isn't

        mask_texture = Screen('MakeTexture', win, annulus_mask);
        %%%%%%% Annulus Code %%%%%%%%%%%%%%%%%%%
    end
    
    % Vector that holds one line of the gratings. The other lines will be drawn by the GPU
    x = meshgrid(-gratings_spatial_period/2:gratings_spatial_period/2,1);
    x(gratings_spatial_period/2 + 1) = []; % remove the zero at the center of the vector
    single_period = sign(x); % compute the period
    num_periods = ceil(texture_width / gratings_spatial_period); % number of periods needed to cover the whole texture
    sign_function = single_period;
    for a=1:(num_periods - 1) % repeats the period over the whole texture width
        sign_function = horzcat(sign_function, single_period);
    end
    
    % 3d matrix that holds the full color gratings
    gratings = zeros(1, length(sign_function), 3);
    % Add colors
    gratings(1,:,1) = 128 + gratings_deviation_red * sign_function;
    gratings(1,:,2) = 128 + gratings_deviation_green * sign_function;
    gratings(1,:,3) = 128 + gratings_deviation_blue * sign_function;
    
    new_VBL_timestamps = zeros(1,num_textures);
    srcRect = [0 0 texture_width texture_height];
    
    while texture_cnt < num_textures && keyIsDown == 0
        
        % Checks if any key has been pressed
% %         keyIsDown = KbCheck(1);
                keyIsDown = KbCheck;

        
        grating_texture = Screen('MakeTexture', win, gratings);
        % Draw grating texture, rotated by "angle":
        Screen('DrawTexture', win, grating_texture, srcRect, [], gratings_angle);
        
        if annulus_radius ~= 0
            % Draw mask texture, overlaps the gratings texture
            Screen('DrawTexture', win, mask_texture, [], []);
        end
                
        % Define Timestamps:
        VBLTimestamp = Screen('Flip', win, 0, 0, 0);
        new_VBL_timestamps(texture_cnt + 1) = VBLTimestamp;
        
        % sends out a pulse on the RTS serial line
        
        TTL_pulse(texture_cnt, handle, TTL_interval);
        texture_cnt = texture_cnt + 1;
        
        
        % gratings = circshift(gratings, [0 gratings_speed]);
        
        
        if gratings_speed >= 1
            gratings = circshift(gratings, [0 gratings_speed]);
        elseif gratings_speed < 1 && mod(texture_cnt, 1/gratings_speed) == 0 % if mod ~= 0, gratings = gratings
            gratings = circshift(gratings, [0 1]); % if texcont = 2 and speed = 1/2, shift by 1? 
        end



Screen('Close', grating_texture);
    end
    
    % grey background
    Screen('Flip', win, 0, 0, 0);
    % % % % % %     WaitSecs(3);
    % waits for all keys to be released, avoids stick keys
    while KbCheck; end
    
    % remove trailing zeros, if any
    new_VBL_timestamps = new_VBL_timestamps(new_VBL_timestamps~=0);
end