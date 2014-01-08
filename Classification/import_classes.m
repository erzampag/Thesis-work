function classes = import_classes(pfile, idListPreD)

% Get Classes
% ENZ, Fall 2012
classInfo = pfile.getClassIDs(); % returns a JavaHashmap Object

classes = cell(length(idListPreD), 2); % each row is [ID ClassInfo]
for i = 1:length(idListPreD)
    classes{i,1} = idListPreD(i); % because it's a cell it needs a curly bracket
    classes{i,2} = classInfo.get(int32(idListPreD(i)));
end

end