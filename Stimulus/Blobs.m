
try
    HideCursor;
% % %     [win, winRect] = Screen('OpenWindow', 0, 128);
    [win, winRect] = Screen('OpenWindow', 0, 128);

    
% % %     priorityLevel=MaxPriority(win);
% % % 
    
% % % % %     Screen('SetVideoCaptureParameter',grabber, ':CodecType=xvidenc')
% % % % % % % % %     movie = Screen('CreateMovie', win, '/Users/erinzampaglione/Desktop/test2', appsrc name=ptbvideoappsrc do-timestamp=0 stream-type=0 max-bytes=0 block=1 is-live=0 emit-signals=0 ! capsfilter caps="video/x-raw-rgb, bpp=(int)32, depth=(int)32, endianess=(int)4321, red_mask=(int)16711680, green_mask=(int)65280, blue_mask=(int)255, width=(int)1920, height=(int)1200, framerate=30/1 "5! videorate ! ffmpegcolorspace !  x264enc key-int-max=30 speed-preset=1  ! qtmux faststart=1 movie-timescale=1000 ! filesink name=ptbfilesink async=0 location=MyTestMovie.mov);
% % % % % movie = Screen('CreateMovie', win,  '/Users/erinzampaglione/Desktop/test2');
% % % % % %%%':CodecType=VideoCodec=xvidenc',
% % % % appsrc name=ptbvideoappsrc do-timestamp=0 stream-type=0 max-bytes=0 block=1 is-live=0 emit-signals=0 ! capsfilter caps="video/x-raw-rgb, bpp=(int)32, depth=(int)32, endianess=(int)4321, red_mask=(int)16711680, green_mask=(int)65280, blue_mask=(int)255, width=(int)1920, height=(int)1200, framerate=30/1 "5! videorate ! ffmpegcolorspace !  x264enc key-int-max=30 speed-preset=1  ! qtmux faststart=1 movie-timescale=1000 ! filesink name=ptbfilesink async=0 location=MyTestMovie.mov
% % %      movie = Screen('CreateMovie', win, '/Users/erinzampaglione/Desktop/test2.mov', 'profile', []);



    black = BlackIndex(win);
    StartPos = round(rand(1,2)*500);
    NewDirection = [];
    
    
    imageArray = [];
    
    % % % % % %        Example 1: Generate values from the uniform distribution on the
    % % % % % %        interval [a, b].
    % % % % % %           r = a + (b-a).*rand(100,1);
    
    numblob = 100;
    blobsize = 40;
    StartPos = [round(winRect(3).*(rand(numblob,1))) round(winRect(4).*(rand(numblob,1)))]; % Take up the entire screen
    
    
    
    for i = 1: 3 % number of direction changes
        
        % % %     StartPos = round(rand(1,2)*40);
        % % %     NewDirection = StartPos + round(-1 + 2.*rand(2,1))';
        
        % % %     StartPos = round(rand(40,2));
        % % %     NewDirection = StartPos + round(-1 + 2.*rand(40,2));
        
        % % %      NewDirection = round(-1 + 2.*rand(40,2));
        NewDirection = [];
        while length(NewDirection) < length(StartPos)
            rand_direction = round(-1 + 2.*rand(1,2));
            if rand_direction(1) == 0 && rand_direction(2) == 0;
                continue
            else
                NewDirection = [NewDirection; rand_direction]; % one of 8 directions
            end
        end
        
        
        for k = 1:250 % number of pixels before direction change
            if mod(k, 1) ~=0 % slows down speed from 1 pixel/frame
                StartPos = StartPos;
            else
                StartPos = StartPos+ NewDirection;
                for j = 1 : length(StartPos) % is there an easier way to wrap?
                    if StartPos(j,1) > winRect(3) - blobsize/2
                        StartPos(j,1) = StartPos(j,1) - winRect(3);
                    end
                    if StartPos(j,1) < winRect(1)- blobsize/2
                        StartPos(j,1) = StartPos(j,1) + winRect(3);
                    end
                    if StartPos(j,2) > winRect(4) - blobsize/2
                        StartPos(j,2) = StartPos(j,2) - winRect(4);
                    end
                    if StartPos(j,2) < winRect(2) - blobsize/2
                        StartPos(j,2) = StartPos(j,2) + winRect(4);
                    end
                end
            end
            
            
            EndPos = StartPos + blobsize;
            for j = 1: length(StartPos) % Draw blobs
                Start = StartPos(j, :);
                End = EndPos(j, :);
                Screen('FillArc', win, 0, [Start End], 0, 360);
            end
            
            Screen('Flip', win);
            
            
% % % % %                          Screen('AddFrameToMovie', win)
% % % % % % %             if mod(k,10) ==0
% % % % % % %                      imageArray = [imageArray; {Screen('GetImage', win)}];
% % % % % % %             end
        end
    end
    
    ShowCursor;
    Screen('CloseAll');
    
% % % % % % % %         Screen('FinalizeMovie', movie);
    
    
% % % % % % % %     	%Creates the .gif
% % % % % % % %     	delayTime = 1/60; %Screen refresh rate of 60Hz
% % % % % % % %     	for i=1:length(imageArray)
% % % % % % % %     		%Gifs can't take RBG matrices: they have to be specified with the pixels as indices into a colormap
% % % % % % % %     		%See the help for imwrite for more details
% % % % % % % %     		[y, newmap] = cmunique(imageArray{i});
% % % % % % % %     
% % % % % % % %     		%Creates a .gif animation - makes first frame, then appends the rest
% % % % % % % %     		if i==1
% % % % % % % %     			imwrite(y, newmap, '/Users/erinzampaglione/Desktop/test3.gif');
% % % % % % % %     		else
% % % % % % % %     			imwrite(y, newmap, '/Users/erinzampaglione/Desktop/test3.gif', 'DelayTime', delayTime, 'WriteMode', 'append');
% % % % % % % %     		end
% % % % % % % %     	end
catch
    ShowCursor;
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
