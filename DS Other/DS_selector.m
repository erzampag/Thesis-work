

figure
        x_fake=[0 1 0 -1];
        y_fake=[1 0 -1 0];
        
        h_fake=compass(x_fake,y_fake, '--');
        hold on;
        
        for j = 1 : length(DScellsindex_new)
            %     for i = 1 : length(AllCells)
            %     if AllCells(i,2) >0.4 %normalized magnitude vector
            %         if DSI(i) > 0.5 && ratio(i) > 1 % Direction selectivity index and f2/f1
            %  if strncmp(textdata{i,4}, 'W', 1)   
                
                ch = compass((AllCells(DScellsindex_new(j),2).*cos(AllCells(DScellsindex_new(j),1))),...
                    (AllCells(DScellsindex_new(j),2).*sin(AllCells(DScellsindex_new(j),1))), '-k');

            %Removing the label
            % set(findall(gcf, 'String', '30','String','60', 'String', '180') ,'String', ' ')
            
            %polar2(AllCells(i,1), AllCells(i,2),[0 1], '-k');% this is in the Matlab folder
            set(ch,'LineWidth',1.5)
            hold on

        end
        
        
        [x,y] = ginput(2);
    
    width = x(2) - x(1);
    height = y(2) - y(1);
        
        
        
        
        
        
        
        
        