try
    [win, winRect] = Screen('OpenWindow', 0, 128);
    
    screen_width = winRect(3);
    screen_height = winRect(4);
    
    black = BlackIndex(win);
    white = WhiteIndex(win);
    
    dark = black*ones(screen_height, screen_width);
    light = white*ones(screen_height, screen_width);
    
    dark_texture = Screen('MakeTexture', win, dark);
    light_texture = Screen('MakeTexture', win, light);
    
    tic
    
    for i = 1:3;
        
        Screen('DrawTexture',win, dark_texture);
        Screen('Flip', win);
        
        WaitSecs(1);
        
        Screen('Flip', win);
        
        WaitSecs(1);
        
        Screen('DrawTexture',win, light_texture);
        Screen('Flip', win);
        
        WaitSecs(1);
        
        Screen('Flip', win);
        
        WaitSecs(1);
    end
    
    toc
    
    ShowCursor;
    Screen('CloseAll');
    
catch
    ShowCursor;
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
    
    