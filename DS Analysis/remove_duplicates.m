function [idList, ind] = remove_duplicates(idListPreD, classes)

% Make IDlist have only removed duplicates
Dup = [];
for i = 1:length(idListPreD)
    if strcmp(classes{i,2}, 'All/Duplicates')
        Dup = [Dup i]; % making a dupilicates index
    end
end

%Different way of making mag index -
% mag4ind = find(data(:,4) >=0.4); % making an index based on a minimum magnitude

ind = setdiff(1:length(idListPreD), Dup); % indices of neurons after duplication removal
% ind4 = setdiff(mag4ind, Dup); % indicies of min. mag. neurons after dup removal


% Use indices of neurons to make a duplicate-removed list of neurons

idList = [];
idList = idListPreD(ind);  % Duplicates only
% idList  = idListPreD(ind4); % Duplicates AND Vision-determined DS cells


end
