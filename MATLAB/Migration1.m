function [MigRow, MigColumn] = Migration1(rowIndex,columnIndex,rowNum,columnNum,BC)
MigRow = rowIndex; %Preassign the rwo index as the current row
MigColumn = columnIndex; %Preassign the column index as the current column
P = rand; %Genearte a random number between 0 and 1

%{
Target Sites arround lattice site L.
   ___  ___
  / 3 \/ 5 \
  \___/\___/__
/ 1 \/ L \/ 2 \
\___/\___/\___/ 
  / 4 \/ 6 \
  \___/\___/

%}

if BC==1
    if P<1/6 %Target site 1
        if rowIndex==1
            MigRow = rowNum;
        else
            MigRow = MigRow - 1;
        end
    elseif P>=1/6 && P<1/3 %Target site 2
        if rowIndex==1
            MigRow = MigRow - 1 + rowNum;
        else
            MigRow = MigRow - 1;
        end
        if columnIndex==columnNum-1
            MigColumn = 1;
        elseif columnIndex<columnNum-1
            MigColumn = MigColumn + 1;
        end
    elseif P>=1/3 && P<1/2 %Target site 3
        if columnIndex==columnNum-1
            MigColumn = 1;
        elseif columnIndex<columnNum-1
            MigColumn = MigColumn + 1;
        end
    elseif P>=1/2 && P<2/3 %Target site 4
        if rowIndex==rowNum
            MigRow = 1;
        else
            MigRow = MigRow + 1;
        end
        if columnIndex==columnNum-1
            MigColumn = 1;
        elseif columnIndex<columnNum-1
            MigColumn = MigColumn + 1;
        end
    elseif P>=2/3 && P<5/6 %Target site 5
        if rowIndex==rowNum
            MigRow = 1;
        else
            MigRow = MigRow + 1;
        end
    elseif P>=5/6 %Target site 6
        if columnIndex==1
            MigColumn = columnNum - 1;
        else
            MigColumn = MigColumn - 1;
        end
    end
elseif BC==2
    if P<1/6 %Target site 1
        if rowIndex==1
            MigRow = rowIndex;
        else
            MigRow = MigRow - 1;
        end
    elseif P>=1/6 && P<1/3 %Target site 2
        if rowIndex==1
            MigRow = rowIndex;
        else
            MigRow = MigRow - 1;
        end
        if columnIndex==columnNum-1
            MigColumn = columnIndex;
        elseif columnIndex<columnNum-1
            MigColumn = MigColumn + 1;
        end
    elseif P>=1/3 && P<1/2 %Target site 3
        if columnIndex==columnNum-1
            MigColumn = columnIndex;
        elseif columnIndex<columnNum-1
            MigColumn = MigColumn + 1;
        end
    elseif P>=1/2 && P<2/3 %Target site 4
        if rowIndex==rowNum
            MigRow = rowIndex;
        else
            MigRow = MigRow + 1;
        end
        if columnIndex==columnNum-1
            MigColumn = columnIndex;
        elseif columnIndex<columnNum-1
            MigColumn = MigColumn + 1;
        end
    elseif P>=2/3 && P<5/6 %Target site 5
        if rowIndex==rowNum
            MigRow = rowIndex;
        else
            MigRow = MigRow + 1;
        end
    elseif P>=5/6 %Target site 6
        if columnIndex==1
            MigColumn = columnIndex;
        else
            MigColumn = MigColumn - 1;
        end
    end    
end
