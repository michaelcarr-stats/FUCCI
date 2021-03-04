function [MigRow, MigColumn] = Migration2(rowIndex,columnIndex,rowNum,columnNum)
MigRow = rowIndex; %Preassign the rwo index as the current row
MigColumn = columnIndex; %Preassign the column index as the current column
P = rand; %Genearte a random number between 0 and 1
BC = 2; %Boundary condition: 1.Periodic; 2.No flux
if BC==1
    if P<1/6 %Target site 1
        if rowIndex==1
            MigRow = rowNum;
        else
            MigRow = MigRow - 1;
        end
        if columnIndex==1
            MigColumn = columnNum - 1;
        else
            MigColumn = MigColumn - 1;
        end
    elseif P>=1/6 && P<1/3 %Target site 2
        if rowIndex==1
            MigRow = rowNum;
        else
            MigRow = MigRow - 1;
        end
        if columnIndex==columnNum
            MigColumn = 1;
        end
    elseif P>=1/3 && P<1/2 %Target site 3
        if columnIndex==columnNum
            MigColumn = 1;
        else
            MigColumn = MigColumn + 1;
        end
    elseif P>=1/2 && P<2/3 %Target site 4
        if rowIndex==rowNum
            MigRow = 1;
        else
            MigRow = MigRow + 1;
        end
        if columnIndex==columnNum
            MigColumn = 1;
        end
    elseif P>=2/3 && P<5/6 %Target site 5
        if rowIndex==rowNum
            MigRow = 1;
        else
            MigRow = MigRow + 1;
        end
        if columnIndex==1
            MigColumn = columnNum - 1;
        else
            MigColumn = MigColumn - 1;
        end
    elseif P>=5/6 %Target site 6
        if columnIndex==1
            MigColumn = columnNum;
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
        if columnIndex==1
            MigColumn = columnIndex;
        else
            MigColumn = MigColumn - 1;
        end
    elseif P>=1/6 && P<1/3 %Target site 2
        if rowIndex==1
            MigRow = rowIndex;
        else
            MigRow = MigRow - 1;
        end
        if columnIndex==columnNum
            MigColumn = columnIndex;
        end
    elseif P>=1/3 && P<1/2 %Target site 3
        if columnIndex==columnNum
            MigColumn = columnIndex;
        else
            MigColumn = MigColumn + 1;
        end
    elseif P>=1/2 && P<2/3 %Target site 4
        if rowIndex==rowNum
            MigRow = rowIndex;
        else
            MigRow = MigRow + 1;
        end
        if columnIndex==columnNum
            MigColumn = columnIndex;
        end
    elseif P>=2/3 && P<5/6 %Target site 5
        if rowIndex==rowNum
            MigRow = rowIndex;
        else
            MigRow = MigRow + 1;
        end
        if columnIndex==1
            MigColumn = columnIndex;
        else
            MigColumn = MigColumn - 1;
        end
    elseif P>=5/6 %Target site 6
        if columnIndex==1
            MigColumn = columnIndex;
        else
            MigColumn = MigColumn - 1;
        end
    end
end