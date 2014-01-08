
% Look at DriftDemo6

try
    HideCursor;
    [win, winRect] = Screen('OpenWindow', 0, 128);
    
    %     global screen_width
    %     global screen_height
    
    screen_width = winRect(3);
    screen_height = winRect(4);
    
    black = BlackIndex(win);
    white = WhiteIndex(win);
    
    
    
    
    gratings_spatial_period = 100;
    gratings_temporal_period = 64;
    
    gratings_deviation_red = 127;
    gratings_deviation_blue = 127;
    gratings_deviation_green = 127;
    
    
    gratings_speed = gratings_spatial_period / gratings_temporal_period;
    
    
    
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
    
    
    % %     sca
    % %     keyboard
    
    %     Screen('FillOval',win, 0,[0,0,100,100])
    
    
    
    % % % % % % % % %     annulus_positions = multi_annulus_positions_generator(100, 3, '/Users/erinzampaglione/Documents/Stimulus_output/');
    % % % % % % % % %     annulus_radius = 100;
    % % % % % % % % %
    % % % % % % % % %     %%%%%%% Annulus Code %%%%%%%%%%%%%%%%%%%
    % % % % % % % % %     annulus = 255*(~Circle(annulus_radius)); % points inside annulus are 0, outside are 255
    % % % % % % % % %     annulus_texture = repmat(255, screen_height, screen_width); % creates the background
    % % % % % % % % %     for a = 1:3 % in Bars and Gratings, it runs multi_annulus_positions_generator.m instead of Daniele's
    % % % % % % % % %         annulus_texture(annulus_positions(a,1):annulus_positions(a,2),annulus_positions(a,3):annulus_positions(a,4)) = annulus; % puts the circle in the background
    % % % % % % % % %     end
    % % % % % % % % %     annulus_mask = ones(screen_height, screen_width, 2) * 128; % first layer, gray as the gratings background
    % % % % % % % % %     annulus_mask(:,:,2) = annulus_texture; % second layer, 0 where the circle is, 255 where it isn't
    % % % % % % % % %
    % % % % % % % % %     mask_texture = Screen('MakeTexture', win, annulus_mask);
    % % % % % % % % %     %%%%%%% Annulus Code %%%%%%%%%%%%%%%%%%%
    
    % % % % % % % % %         mask_texture = Screen('MakeTexture', win, annulus_mask);
    % % % % % % % % %         Screen('DrawTexture', win, mask_texture, [], []);
    % % % % % % % % %
    %%
    
    
    % Make sure this GPU supports shading at all:
    AssertGLSL;
    
    % Enable alpha blending for typical drawing of masked textures:
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Create a special texture drawing shader for masked texture drawing:
    glsl = MakeTextureDrawShader(win, 'SeparateAlphaChannel');
    
    
    % Create circular aperture for the alpha-channel:
    texsize =300;
    [x,y]=meshgrid(-texsize:texsize, -texsize:texsize);
    circle = white * (x.^2 + y.^2 <= (texsize)^2);
    
    % Set 2nd channel (the alpha channel) of 'grating' to the aperture
    % defined in 'circle':
    gratings2 = gratings;
    gratings2(:,:,2) = 0;
    gratings2(1:2*texsize+1, 1:2*texsize+1, 2) = circle;
    
    
    
    
    
    
    %%
    texture_cnt = 0;
    
    while texture_cnt < 1000
        
        %%%% && keyIsDown == 0
        
        % Checks if any key has been pressed
        % %         keyIsDown = KbCheck(1);
        keyIsDown = KbCheck;
        
        %%
            % Set 2nd channel (the alpha channel) of 'grating' to the aperture
    % defined in 'circle':
    gratings2 = gratings;
    gratings2(:,:,2) = 0;
    gratings2(1:2*texsize+1, 1:2*texsize+1, 2) = circle;
        %%
        
        grating_texture = Screen('MakeTexture', win, gratings);
        
        
        % Draw grating texture, rotated by "angle":
        Screen('DrawTexture', win, grating_texture, srcRect)
        
        %%
        % Store alpha-masked grating in texture and attach the special 'glsl'
        % texture shader to it:
        gratingtex1 = Screen('MakeTexture', win, gratings2, [], [], [], [], glsl);
        %         Screen('DrawTexture', win, gratingtex1, srcRect, [], 90, [], [], [], [], [], [0, yoffset, 0, 0]);
        Screen('DrawTexture', win, gratingtex1, srcRect,[], 90);
        
        
        %%
        
        
        
        %         Screen('FillOval',win, gratings ,[0,0,300,300])
        
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
        
        if mod(texture_cnt, 30) == 0
            temp_rand = (1-2*rand(1));
            if temp_rand < 0
                gratings_speed = -1;
            else
                gratings_speed = 1;
            end
            
            gratings = circshift(gratings, [0 gratings_speed]);
            
        end
        
        
        Screen('Close', grating_texture);
        Screen('Close', gratingtex1);
    end
    
    
    % %     pause
    % %     WaitSecs(2);
    ShowCursor;
    Screen('CloseAll')
    
catch
    ShowCursor;
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end