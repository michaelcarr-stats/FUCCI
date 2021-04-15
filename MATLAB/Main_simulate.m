 function [SummaryStatData, ExitSimStatus] = Main_simulate(theta_prop, s, T_record, CellTracking)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                 Setup                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ntrack = s.ntrack; %number of cells to track
    simuNum = s.simuNum; %The total number of realisations
    delta = s.delta; %Lattice spacing [mum]
    t = s.t; %Time clock starting with 0
    tstop = s.tstop; % Total iterations at each realisation
    Xmax = s.Xmax; %width of domain
    Ymax = s.Ymax; %heigth of domain
    columnNum = s.columnNum; %Corresponding total column nodes
    rowNum = s.rowNum; %Corresponding total row nodes
    BC = s.BC; %Boundary condition: 1.Periodic; 2.No flux
    domain_x = s.domain_x; %x coordinates for all the nodes
    domain_y = s.domain_y; %y coordinates for all the nodes
    domain = s.domain;
    NstartRed = s.NstartRed;
    NstartYellow = s.NstartYellow;
    NstartGreen = s.NstartGreen;
    
    %model parameters
    Pr = theta_prop(4); % migration rate of red per time unit
    Py = theta_prop(5); % migration rate of yellow per time unit
    Pg = theta_prop(6); % migration rate of green per time unit
    Kry = theta_prop(1); % transition rate of red to yellow per time unit
    Kyg = theta_prop(2); % transition rate of yellow to green per time  unit
    Kgr = theta_prop(3); % transition rate of green to red per time unit
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             Initialisation                              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    simuIndex = 1;
    
    % update population counts
    Nred = NstartRed;
    Nyellow = NstartYellow;
    Ngreen = NstartGreen;
    
    SummaryStatData = NaN;
    
    %initialise output to handle exit simulation event
    ExitSimStatus = true; %assigned false if simulation is completed 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             Main framework                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    while simuIndex <= simuNum  %If the current simulation index is not larger than the total realisation number
        
        transID = 0; %reset transition ID for next iteration
        migFailed = true; %determine if migration was aborted
         
        
        RowPosCell = s.RowPosCell; %Row position of tracked cells
        ColPosCell = s.ColPosCell; %column position of tracked cells
        CellSelected = s.CellSelected;
        CellSelectedStart = CellSelected;
        
        while t < tstop  %Iterations from t = delta t to T
            
            
            % calculate total propensity function
            ar = Pr * Nred; % red movement
            ay = Py * Nyellow; % yellow movement
            ag = Pg * Ngreen; % green movement
            tr = Kry * Nred; % transition red to yellow
            ty = Kyg * Nyellow; % transition yellow to green
            tg = Kgr * Ngreen; % transition green to red
            a0 = ar + ay + ag + tr + ty + tg; % total propensity function

            % calculate new time step
            tau = (1/a0)*log(1/rand);
            R = rand;

            % Determine what event occurs by splitting the interval up based on
            % the propensity functions          
                        
            % Red migration
            if R <= ar/a0 
                [rowIndexs,columnIndexs] = find(domain == 1); %find locations of red cells
                Index = randi([1 length(rowIndexs)]); %randomly choose location
                rowIndex = rowIndexs(Index); %assign rowIndex to corresponding location
                columnIndex = columnIndexs(Index); %assign columnIndex to corresponding location 
                if rem(rowNum-rowIndex,2)==0 %If the selected agent is at bottom/(bottom-2)/(bottom-4)/... rows
                    [migPosition(1), migPosition(2)] = Migration1(rowIndex,columnIndex,rowNum,columnNum,BC);
                    if domain(migPosition(1),migPosition(2))==0 && ...%Check (1)if the target site is vacant; (2)if the target site is within the actual domain
                            domain_x(migPosition(1),migPosition(2))<=Xmax && ...
                            domain_y(migPosition(1),migPosition(2))<=Ymax
                        domain(rowIndex,columnIndex) = 0; %Remove agent at the prevous site
                        domain(migPosition(1),migPosition(2)) = 1; %Move agent to the new site
                        migFailed = false;
                    end
                else %If the selected agent is at (bottom-1)/(bottom-3)/(bottom-5)/... rows
                    [migPosition(1), migPosition(2)] = Migration2(rowIndex,columnIndex,rowNum,columnNum,BC);
                    if domain(migPosition(1),migPosition(2))==0 && ...%Check (1)if the target site is vacant; (2)if the target site is within the actual domain
                            domain_x(migPosition(1),migPosition(2))<=Xmax && ...
                            domain_y(migPosition(1),migPosition(2))<=Ymax
                        domain(rowIndex,columnIndex) = 0; %Remove agent at prevous site
                        domain(migPosition(1),migPosition(2)) = 1; %Move agent to the new site
                        migFailed = false;
                    end
                end

                % Yellow migration
            elseif R > ar/a0 && R <= (ay+ar)/a0
                [rowIndexs,columnIndexs] = find(domain == 2); %find locations of yellow cells
                Index = randi([1 length(rowIndexs)]); %randomly choose location
                rowIndex = rowIndexs(Index); %assign rowIndex to corresponding location
                columnIndex = columnIndexs(Index); %assign columnIndex to corresponding location 
                if rem(rowNum-rowIndex,2)==0  %If the selected agent is at bottom/(bottom-2)/(bottom-4)/... rows
                    [migPosition(1), migPosition(2)] = Migration1(rowIndex,columnIndex,rowNum,columnNum,BC);
                    if domain(migPosition(1),migPosition(2))==0 && ...%Check (1)if the target site is vacant; (2)if the target site is within the actual domain
                            domain_x(migPosition(1),migPosition(2))<=Xmax && ...
                            domain_y(migPosition(1),migPosition(2))<=Ymax
                        domain(rowIndex,columnIndex) = 0; %Remove agent at the prevous site
                        domain(migPosition(1),migPosition(2)) = 2; %Move agent to the new site
                        migFailed = false;
                    end
                else %If the selected agent is at (bottom-1)/(bottom-3)/(bottom-5)/... rows
                    [migPosition(1), migPosition(2)] = Migration2(rowIndex,columnIndex,rowNum,columnNum,BC);
                    if domain(migPosition(1),migPosition(2))==0 && ...%Check (1)if the target site is vacant; (2)if the target site is within the actual domain
                            domain_x(migPosition(1),migPosition(2))<=Xmax && ...
                            domain_y(migPosition(1),migPosition(2))<=Ymax
                        domain(rowIndex,columnIndex) = 0; %Remove agent at prevous site
                        domain(migPosition(1),migPosition(2)) = 2; %Move agent to the new site
                        migFailed = false;
                    end
                end

                % Green migration
            elseif R > (ay+ar)/a0 && R <=(ar+ay+ag)/a0
                [rowIndexs,columnIndexs] = find(domain == 3); %find locations of green cells
                Index = randi([1 length(rowIndexs)]); %randomly choose location
                rowIndex = rowIndexs(Index); %assign rowIndex to corresponding location
                columnIndex = columnIndexs(Index); %assign columnIndex to corresponding location 
                if rem(rowNum-rowIndex,2)==0  %If the selected agent is at bottom/(bottom-2)/(bottom-4)/... rows
                    [migPosition(1), migPosition(2)] = Migration1(rowIndex,columnIndex,rowNum,columnNum,BC);
                    if domain(migPosition(1),migPosition(2))==0 && ...%Check (1)if the target site is vacant; (2)if the target site is within the actual domain
                            domain_x(migPosition(1),migPosition(2))<=Xmax && ...
                            domain_y(migPosition(1),migPosition(2))<=Ymax
                        domain(rowIndex,columnIndex) = 0; %Remove agent at the prevous site
                        domain(migPosition(1),migPosition(2)) = 3; %Move agent to the new site
                        migFailed = false;
                    end
                else %If the selected agent is at (bottom-1)/(bottom-3)/(bottom-5)/... rows
                    [migPosition(1), migPosition(2)] = Migration2(rowIndex,columnIndex,rowNum,columnNum,BC);
                    if domain(migPosition(1),migPosition(2))==0 && ...%Check (1)if the target site is vacant; (2)if the target site is within the actual domain
                            domain_x(migPosition(1),migPosition(2))<=Xmax && ...
                            domain_y(migPosition(1),migPosition(2))<=Ymax
                        domain(rowIndex,columnIndex) = 0; %Remove agent at prevous site
                        domain(migPosition(1),migPosition(2)) = 3; %Move agent to the new site
                        migFailed = false;
                    end
                end
                

                % Transition red to yellow
            elseif R > (ar+ay+ag)/a0 && R <= (ar+ay+ag+tr)/a0
                [rowIndexs,columnIndexs] = find(domain == 1); %find locations of red cells
                Index = randi([1 length(rowIndexs)]); %randomly choose location
                rowIndex = rowIndexs(Index); %assign rowIndex to corresponding location
                columnIndex = columnIndexs(Index); %assign columnIndex to corresponding location 
                transID = 2; %identify if cell has changed color
                domain(rowIndex, columnIndex) = 2; % Change domain entry to be a yellow agent
                Nyellow = Nyellow + 1; % update yellow population
                Nred = Nred - 1; % update red population
                
                
                % Transition yellow to green
            elseif R > (ar+ay+ag+tr)/a0 && R <= (ar+ay+ag+tr+ty)/a0 
                [rowIndexs,columnIndexs] = find(domain == 2); %find locations of yellow cells
                Index = randi([1 length(rowIndexs)]); %randomly choose location
                rowIndex = rowIndexs(Index); %assign rowIndex to corresponding location
                columnIndex = columnIndexs(Index); %assign columnIndex to corresponding location 
                transID = 3; %identify if cell has changed color
                domain(rowIndex, columnIndex) = 3; % Change domain entry to be green agent
                Ngreen = Ngreen + 1; % update green population
                Nyellow = Nyellow - 1; % update yellow population
                
                % Transition green to red and proliferation
            elseif R > (ar+ay+ag+tr+ty)/a0 && R <= (ar+ay+ag+tr+ty+tg)/a0 
                [rowIndexs,columnIndexs] = find(domain == 3); %find locations of green cells
                Index = randi([1 length(rowIndexs)]); %randomly choose location
                rowIndex = rowIndexs(Index); %assign rowIndex to corresponding location
                columnIndex = columnIndexs(Index); %assign columnIndex to corresponding location 
                                
                % move the daughter cell to the new position
                if rem(rowNum-rowIndex,2)==0  %If the selected agent is at bottom/(bottom-2)/(bottom-4)/... rows
                    [migPosition(1), migPosition(2)] = Migration1(rowIndex,columnIndex,rowNum,columnNum,BC);
                    if domain(migPosition(1),migPosition(2))==0 && ...%Check (1)if the target site is vacant; (2)if the target site is within the actual domain
                            domain_x(migPosition(1),migPosition(2))<=Xmax && ...
                            domain_y(migPosition(1),migPosition(2))<=Ymax
                        domain(rowIndex, columnIndex) = 1; % change the original agent to red
                        domain(migPosition(1),migPosition(2)) = 1; % place red daughter agent at the new site
                        Nred = Nred + 2; % increase red cell population
                        Ngreen = Ngreen - 1; % decrease green cell population
                        transID = 1; %identify if cell has changed color
                    end
                else %If the selected agent is at (bottom-1)/(bottom-3)/(bottom-5)/... rows
                    [migPosition(1), migPosition(2)] = Migration2(rowIndex,columnIndex,rowNum,columnNum,BC);
                    if domain(migPosition(1),migPosition(2))==0 && ...%Check (1)if the target site is vacant; (2)if the target site is within the actual domain
                            domain_x(migPosition(1),migPosition(2))<=Xmax && ...
                            domain_y(migPosition(1),migPosition(2))<=Ymax
                        domain(rowIndex, columnIndex) = 1; % change the original agent to red
                        domain(migPosition(1),migPosition(2)) = 1; % place red daughter agent at the new site
                        Nred = Nred + 2; % increase red cell population
                        Ngreen = Ngreen - 1; % decrease green cell population
                        transID = 1; %identify if cell has changed color
                    end
                end
            end
            
            %generate summary statistic data 
            SummaryStatData = GenerateSummaryStatData(CellTracking, SummaryStatData, ntrack, T_record, t, Nred, Nyellow, Ngreen, RowPosCell, ColPosCell, CellSelected, CellSelectedStart, rowIndex, columnIndex, migFailed, transID, delta, simuIndex, simuNum, domain, domain_x); 
            
            if ~migFailed && any((RowPosCell == rowIndex).*(ColPosCell== columnIndex))
                index = (1:ntrack).*(RowPosCell == rowIndex).*(ColPosCell == columnIndex); %determine which tracked cell moved
                index = index(index ~= 0); 
                RowPosCell(1,index) = migPosition(1); %update positions
                ColPosCell(1,index) = migPosition(2); %update positions
            elseif transID ~= 0 && any((RowPosCell == rowIndex).*(ColPosCell== columnIndex))
                index = (1:ntrack).*(RowPosCell == rowIndex).*(ColPosCell == columnIndex); %determine which tracked cell transitioned
                index = index(index ~= 0);
                CellSelected(1,index) = transID; %update phase id
            end
            
            transID = 0; %reset transition ID for next iteration
            migFailed = true; %reset migration abortion indicator
            
            %Update time
            t = t + tau;
            
        end
        
        
        % reset the system if the realisation is not the last one
        if simuIndex ~= simuNum %If the current realisation is not the last one
            
            t = 0; %Reset time     

            % re-initialise the domain
            domain = s.domain;
            NstartRed = s.NstartRed;
            NstartYellow = s.NstartYellow;
            NstartGreen = s.NstartGreen;
            
            % update total population counts
            Nred = NstartRed;
            Nyellow = NstartYellow;
            Ngreen = NstartGreen;
        end
        simuIndex = simuIndex + 1; %Update simulation index for the next realisation
    end
    
   ExitSimStatus = false; %indicate simulation finished successfully
end
