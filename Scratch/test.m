 
  		% make movie
  		window=Screen('OpenWindow', 0, 0);
  		rect=[0 0 200 200];
  		for i=1:100
  			movie(i)=Screen('OpenOffscreenWindow', window, 0, rect);
  			Screen('FillOval', movie(i), 255, [0 0 2 2]*(i-1));
  		end;
 
  		% show movie
  		for i=[1:100 100:-1:1] % forwards and backwards
  			Screen('CopyWindow',movie(i),window,rect,rect);
  			Screen('Flip', window);
  		end;
  		Screen('CloseAll');
 