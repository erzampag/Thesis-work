function []= color_test(screen_id, red, green, blue, secs)

% screen_id

try
    [win, winRect] = Screen('OpenWindow', screen_id, 128);
    
    screen_width = winRect(3);
    screen_height = winRect(4);
    
    frame = zeros(screen_height, screen_width, 3);
    
    frame(:, :, 1) = red*255*ones(screen_height, screen_width); % red
    frame(:, :, 2) = green*255*ones(screen_height, screen_width); % green
    frame(:, :, 3) = blue*255*ones(screen_height, screen_width); %blue
    
    
    color_texture = Screen('MakeTexture', win, frame);
    
    
    Screen('DrawTexture',win, color_texture);
    Screen('Flip', win);
    
%     WaitSecs(secs);
    keyisDown =0;
    while keyisDown == 0 
        keyisDown = KbCheck(-3);
        
    end
    
    ShowCursor;
    Screen('CloseAll');
    
catch
    ShowCursor;
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
