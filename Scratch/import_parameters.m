function [data, EIx, EIy] = import_parameters(pfile, idListPreD)
% STUFF OUT OF THE PARAMETERS FILE

% Get Data

data = zeros(length(idListPreD), 5);
%
% for i = 1:length(idListPreD)
%     data(i,1) = idListPreD(i);
%     data(i,2) = pfile.getDoubleCell(idListPreD(i), 'xOffDS');
%     data(i,3) = pfile.getDoubleCell(idListPreD(i), 'yOffDS');
%     data(i,4) = pfile.getDoubleCell(idListPreD(i), 'magOffDS');
%     data(i,5) = pfile.getDoubleCell(idListPreD(i), 'angOffDS');
% end





% keyboard
% % % % % % % % % % % % 
EIx = zeros(length(idListPreD),1);
EIy = zeros(length(idListPreD),1);
% % % % % % % % % % % % % % EI_all = zeros(length(idListPreD), 3);
% % % % % % % % % % % % 
for i = 1:length(idListPreD)
    
%         EI_all(i,1) = idListPreD(i);
%         EI_all(i,2) = pfile.getDoubleCell(idListPreD(i), 'EIx0');
%         EI_all(i,3) = pfile.getDoubleCell(idListPreD(i), 'EIy0');
    EIx(i,1) = pfile.getDoubleCell(idListPreD(i), 'EIx0');
    EIy(i,1) = pfile.getDoubleCell(idListPreD(i), 'EIy0');
end
end
