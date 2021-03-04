function [RowPosCell,ColPosCell,CellSelectedIndex] = InitialMotilitySelect(ntrack, domain, domain_x, Xmax)

Index = reshape(1:numel(domain_x), size(domain_x));
Index(domain ~= 1) = NaN; %identify positions which are not red

Index = Index(~isnan(Index));
ValidX = domain_x(Index) > 0 & domain_x(Index) < Xmax;
Index = [Index(ValidX),domain_x(Index(ValidX)), domain_x(Index(ValidX)) < Xmax/2];
IndexLHS = sortrows(Index(Index(:,3) == 1,:), 2, 'descend');
IndexRHS = sortrows(Index(Index(:,3) == 0,:), 2, 'ascend');

CellSelectedIndex = [IndexLHS(1:ceil(ntrack/2),1);IndexRHS(1:ceil(ntrack/2),1)];
CellSelectedIndex = CellSelectedIndex(1:ntrack); %ensure only ntrack cells are tracked
RowPosCell = mod(CellSelectedIndex-1,size(domain_x,1))+1;
ColPosCell = ceil(CellSelectedIndex/size(domain_x,1));
        
end

