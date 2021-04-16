function [Nred, Nyellow, Ngreen] = InitialCellCount(domain, domain_x, domain_y, rowNum, columnNum, Xmax, Ymax)
    % this function counts the initial number of cells in the domain 
    
    Nred = 0;
    Nyellow = 0; 
    Ngreen = 0;

    for i = 1:rowNum
        for j = 1:columnNum
            if (domain_x(i,j) < Xmax && domain_y(i,j) < Ymax) %if cell is in domain
                if domain(i,j) == 1 %if cell is red
                    Nred = Nred + 1;
                elseif domain(i,j) == 2 %if cell is yellow
                    Nyellow = Nyellow + 1;
                elseif domain(i,j) == 3 %if cell is green
                    Ngreen = Ngreen + 1;
                end
            end
        end
    end
end
