try
    HideCursor;
    [win, winRect] = Screen('OpenWindow', 0, 128);
    
    
    
    screen_width = winRect(3);
    screen_height = winRect(4);
    
    black = BlackIndex(win);
    white = WhiteIndex(win);
    
    % From Daneile's Noise.m
    objRect = SetRect(0,0, 40, 20);
    numRects = 1; % Number of rectangles to fit in the window
    dstRect = ArrangeRects(numRects, objRect, winRect);
    [xc, yc] = RectCenter(dstRect);
    dstRect = CenterRectOnPoint(objRect*16, xc, yc);
    
    left = dstRect(1);
    top = dstRect(2);
    right = dstRect(1)+16;
    bottom = dstRect(2)+16;
    newRect = [];
    
    for i = 0:16:304 % moving vertically
        
        
        top = dstRect(2)+i;
        bottom = dstRect(2)+16+i;
        
        for j = 0:16:624 % moving horizontally
            newRect = [newRect; [left+j top right+j bottom]];
            
            
        end
        
        
    end
    
    
    
    
    gratings_spatial_period = 16;
    gratings_temporal_period = 64;
    
    gratings_deviation_red = 127;
    gratings_deviation_blue = 127;
    gratings_deviation_green = 127;
    
    
    gratings_speed = gratings_spatial_period / gratings_temporal_period;
    
    
    
    % % %     noiseimg = round(255.*rand(20, 40));
    % % %     noiseimg = ceil(16.*rand(20,40));
    
    for i = 1:3
        noiseimg = ceil(16.*rand(800,1));
        
        gratings_angle = noiseimg.*(360/16);
        
        
        % Compute texture dimensions based on gratings_angle and screen_width, screen_height. Reduces angles to first quadrant if needed.
        
        % % %      [texture_width, texture_height] = texture_dimensions(screen_width, screen_height, gratings_angle);
        % % %
        % % %    srcRect = [0 0 texture_width texture_height];
        
        
        
        texture_cnt = 0;
        
        % From Daniele's multiple bar run
        x = meshgrid(-gratings_spatial_period/2:gratings_spatial_period/2,1);
        x(gratings_spatial_period/2 + 1) = []; % remove the zero at the center of the vector
        single_period = sign(x); % compute the period
        num_periods = ceil(screen_width / gratings_spatial_period); % number of periods needed to cover the whole texture
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
        
        % % % %     new_VBL_timestamps = zeros(1,num_textures);
        srcRect = [0 0 screen_width screen_height];
        
        while texture_cnt < 1000
            
            %%%% && keyIsDown == 0
            
            % Checks if any key has been pressed
            % %         keyIsDown = KbCheck(1);
            keyIsDown = KbCheck;
            
            
            grating_texture = Screen('MakeTexture', win, gratings);
            % Draw grating texture, rotated by "angle":
            Screen('DrawTextures', win, grating_texture, newRect', newRect', gratings_angle', [], [], [], [], kPsychUseTextureMatrixForRotation); %%randomly trying to see if this will work? newRect(1,:)'
            Screen('Flip', win, 0, 0, 1);
            
            % % % % %         if annulus_radius ~= 0
            % % % % %             % Draw mask texture, overlaps the gratings texture
            % % % % %             Screen('DrawTexture', win, mask_texture, [], []);
            % % % % %         end
            
            % % % % %         % Define Timestamps:
            % % % % %         VBLTimestamp = Screen('Flip', win, 0, 0, 0);
            % % % % %         new_VBL_timestamps(texture_cnt + 1) = VBLTimestamp;
            % % % % %
            % % % % %         % sends out a pulse on the RTS serial line
            % % % % %
            % % % % %         TTL_pulse(texture_cnt, handle, TTL_interval);
            % % % % %         texture_cnt = texture_cnt + 1;
            
            
            % gratings = circshift(gratings, [0 gratings_speed]);
            
            texture_cnt = texture_cnt +1;
            if gratings_speed >= 1
                gratings = circshift(gratings, [0 gratings_speed]);
            elseif gratings_speed < 1 && mod(texture_cnt, 1/gratings_speed) == 0 % if mod ~= 0, gratings = gratings
                gratings = circshift(gratings, [0 1]); % if texcont = 2 and speed = 1/2, shift by 1?
            end
            
            
            
            Screen('Close', grating_texture);
        end
    end
    
    
    
    % %     STuff that just draws color
    % % % % % % % % % % % % % % % % % % %  texture = Screen('MakeTexture', win, noiseimg);
    % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % %     Screen('DrawTexture', win, texture, [], dstRect, [], 0)
    % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % %     Screen('Flip', win, 0, 0, 1);
    
    % % %     sca
    % % %     keyboard
    
    pause
    WaitSecs(2);
    ShowCursor;
    Screen('CloseAll')
    
catch
    ShowCursor;
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end