function domain = initialise_domain(domain, domain_x, rowNum, columnNum, initialDensity, Xmax, initialL, percentR, percentY, Nred_0, Nyellow_0, Ngreen_0, SetCells)
% FUNCTION INITIALISE DOMAIN
% This function initialises the domain. The initial
% seeding of the domain is based on the percentages of 
% red, yellow and green that are specified in the Main 
% script. The populated domain, and start populations 
% of red, yellow and green agents are returned. 

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
        Index = randsample(rowNum*columnNum,rowNum*columnNum,false);
        count = 0;
        for i = 1:length(Index)
            row = mod(Index(i) - 1,rowNum) + 1;
            col = floor((Index(i) - 1)/rowNum) + 1;
            if ((domain_x(row, col) >= 0 && domain_x(row, col) <= initialL/2) || (domain_x(row, col) >= Xmax - initialL/2 && domain_x(row, col) <= Xmax)) % check if inside domain
                count = count + 1;
                if (count <= Nred_0)
                    domain(row,col) = 1;
                elseif (count <= Nred_0 + Nyellow_0)
                    domain(row,col) = 2;
                elseif (count <= Nred_0 + Nyellow_0 + Ngreen_0)
                    domain(row,col) = 3;
                end
            end
        end
    end
    
end

