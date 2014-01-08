

ID_list = [1;257;287;348;351;527;635;677;828;872;931;964;1098;1142;1443;1548;1908;1940;2011;2044;2116;2194;2388;2418;...
    2674;2747;2792;2793;2838;2993;3124;3182;3202;3482;3528;3586;3740;4384;4413;4520;4741;4847;5192;5257;5267;5342;5462;...
    5582;5642;5761;5783;5837;5881;6078;6406;6618;6663;6681;6694;6811;6857;6860;7082;7130;7156;7234;7323;7671;];

index = [];
counter = 1;
for i = 1:length(idList)
    if idList(i) == ID_list(counter)
        counter = counter +1;
        index = [index i];
    end
end



DS_INDEX = find(idList(:) == ID_list);

DS_INDEX = find(DSI(:) > 0.5 & ratio(:) > 1 & noiseratio(:) > 5 & MaxSpikes(:) > 100); % with f2/f1
DS_INDEX = find(DSI(:) > 0.5 & noiseratio(:) > 5 & MaxSpikes(:) > 100); % WITHOUT f2/f1

DS_INDEX = find(DSI(:) > 0.5 & ratio(:) <= 1 & noiseratio(:) > 5 & MaxSpikes(:) > 100); % with f2/f1




figure
x_fake=[0 1 0 -1];
y_fake=[1 0 -1 0];

h_fake=compass(x_fake,y_fake, '--');
hold on;


counter = 1;

for i = 1 : length(AllCells)
    %     if AllCells(i,2) >0.4 %normalized magnitude vector
%     if DSI(i) > 0.5 && ratio(i) > 1 && noiseratio(i) > 5 && MaxSpikes(i) > 100 % Direction selectivity index and f2/f1
%     if DSI(i) > 0.5 && noiseratio(i) > 5 && MaxSpikes(i) > 100 % no f2/f1 > 1 requirement
    if idList(i) == ID_list(counter)
        counter = counter +1; 
        

        
        ch = compass((AllCells(i,2)*cos(AllCells(i,1))), (AllCells(i,2)*sin(AllCells(i,1))), 'ok');
        xd=get(ch,'xdata'); %TURN THESE TWO LINES OFF
        set(ch,'xdata',[xd(1:2),nan,nan,nan]); % TO ADD ARROWHEADS
        % set(ch,'linewidth',4);
        
        %polar2(AllCells(i,1), AllCells(i,2),[0 1], '-k'); % this is in the Matlab folder - wtf polar
        hold on
        set(ch,'LineWidth',1.5)
        
        if ratio(i) > 1
            set(ch, 'Color', 'k');
        else
            set(ch, 'Color', 'b');
        end
        
        
    end
end
set(h_fake,'Visible','off')





[neuronIDs, classes, ~, ~, ~, ~, ~, ~] = ...
    import_neuron_info('/Users/erinzampaglione/Documents/workspace/vision8/Vision.jar',...
    '/Users/erinzampaglione/Documents/processed_data/', '2013-06-10-0', 'data004_mapped_data007');

n_index = [];
s_index = [];
e_index = [];
w_index = [];

for i = 1:length(classes)
    
    if length(classes{i,2}) < 22
        continue
    end
    
    if strcmp(classes{i,2}(20:22), 'nor')
        n_index = [n_index i];
    elseif strcmp(classes{i,2}(20:22), 'sou')
        s_index = [s_index i];
    elseif strcmp(classes{i,2}(20:22), 'eas')
        e_index = [e_index i];
    elseif strcmp(classes{i,2}(20:22), 'wes')
        w_index = [w_index i];
    end
end

DS_INDEX = w_index;
DS_INDEX = index;
figure

hold on

p = polar(AllCells(DS_INDEX,1), AllCells(DS_INDEX,2), 'og');


set(p,'LineWidth',1.5)


