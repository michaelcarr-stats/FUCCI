function domain = initialise_domain(domain, domain_x, rowNum, columnNum, initialDensity, Xmax, initialL, percentR, percentY, Nred_0, Nyellow_0, Ngreen_0, SetCells)

    if (SetCells == false)
        for i=1:rowNum % row loop
            for j=1:columnNum % col loop
                % place initial condition at the left and right end of domain    
                if ((domain_x(i,j) >= 0 && domain_x(i,j) <= initialL/2) || (domain_x(i,j) <= Xmax && domain_x(i,j) >= Xmax - initialL/2)) % check if inside domain
                    if rand <= initialDensity
                        h = rand;
                        if h <= percentR
                            domain(i,j) = 1;
                        elseif h > percentR && h <= (percentR + percentY)
                            domain(i,j) = 2;
                        else
                            domain(i,j) = 3;
                        end
                    end
                end
            end  % end for j
        end % end for i
    else
        Index = randsample(rowNum*columnNum,rowNum*columnNum,false); %randomly sample domain index
        count = 0;
        for i = 1:length(Index)
            row = mod(Index(i) - 1,rowNum) + 1; %compute row index
            col = floor((Index(i) - 1)/rowNum) + 1; %compute column index
            if ((domain_x(row, col) >= 0 && domain_x(row, col) <= initialL/2) || (domain_x(row, col) >= Xmax - initialL/2 && domain_x(row, col) <= Xmax)) % check if inside domain
                count = count + 1;
                if (count <= Nred_0)
                    domain(row,col) = 1; %assign red cell
                elseif (count <= Nred_0 + Nyellow_0) 
                    domain(row,col) = 2; %assign yellow cell
                elseif (count <= Nred_0 + Nyellow_0 + Ngreen_0)
                    domain(row,col) = 3; %assign green cell
                end
            end
        end
    end
    
end

