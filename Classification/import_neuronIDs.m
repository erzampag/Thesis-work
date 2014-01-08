function idListPreD =  import_neuronIDs(neuronFile)

% Import list of neurons IDs from 'neuronFile';
% ENZ, Fall 2012

idListbefore = [];
idListbefore=neuronFile.getIDList();

%Re make the vector of neuronIDs so it doesn't have weird negative numbers
idListPreD = [];

for i = 1:length(idListbefore)
    if idListbefore(i) <=0,
        continue
    else
        idListPreD = [idListPreD, idListbefore(i)];
    end
end
idListPreD = idListPreD';

end