function cumul=cumulative(data,line, col)
%function cumul=cumulative(data,line, col)
% creates cumulative frequency histograms
% line: line in which data starts
% col: column in the data matrix
%
% Marianne Renner oct 09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x=data(line:size(data,1),col);
datacol=sortrows(x(:,1))';	

if ~isempty(datacol)
       proba = linspace(0,1,length(datacol));
       cumul = [datacol', proba'];
end 

clear data datacol proba 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%