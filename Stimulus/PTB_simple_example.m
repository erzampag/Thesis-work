try
    % Open Window
    HideCursor;
    [win, winRect] = Screen('OpenWindow', 0, 128);
    
    % Generate an image from a file
    [I map alpha] = imread('/Users/erinzampaglione/Desktop/owl.png');
    I(:,:,4) = alpha(:,:);
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    owl = Screen('MakeTexture', win, I);
    Screen('DrawTexture', win, owl)
    
    % OR
    
    % Generate a circle of defined location, radius
%     Start = [winRect(3)/2 winRect(4)/2];
%     End = Start+ 40;
%     Screen('FillArc', win, 0, [Start End], 0, 360);
%     
    % Flip Screen to display image drawn in buffer
    Screen('Flip', win);
    
    % Wait, then close screen
    WaitSecs(3);
    
    ShowCursor;
    Screen('CloseAll');
    
catch
    % Often PTB code has a try catch block like this so that the PTB window
    % closes if there is an error while running the code
    ShowCursor;
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end