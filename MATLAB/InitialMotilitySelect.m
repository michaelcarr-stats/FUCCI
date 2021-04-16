function [RowPosCell,ColPosCell,CellSelectedIndex] = InitialMotilitySelect(ntrack, domain, domain_x, Xmax)

Index = reshape(1:numel(domain_x), size(domain_x)); %domain indicies
Index(domain ~= 1) = NaN; %identify positions which are not red

Index = Index(~isnan(Index)); %remove indicies which are not red
ValidX = domain_x(Index) > 0 & domain_x(Index) < Xmax; %determine which cells are in the domain
Index = [Index(ValidX),domain_x(Index(ValidX)), domain_x(Index(ValidX)) < Xmax/2]; % matrix with columns: Indicies in domain, x position of indicies, if on LHS of scratched region
IndexLHS = sortrows(Index(Index(:,3) == 1,:), 2, 'descend'); %sort indicies on LHS of domain by x-position in descending order
IndexRHS = sortrows(Index(Index(:,3) == 0,:), 2, 'ascend'); %sort indicies on LHS of domain by x-position in ascending  order

CellSelectedIndex = [IndexLHS(1:ceil(ntrack/2),1);IndexRHS(1:ceil(ntrack/2),1)]; % select cells on leading edge
CellSelectedIndex = CellSelectedIndex(1:ntrack); %ensure only ntrack cells are tracked
RowPosCell = mod(CellSelectedIndex-1,size(domain_x,1))+1; %compute row index from indecies
ColPosCell = ceil(CellSelectedIndex/size(domain_x,1)); %compute column index from indecies
        
end

