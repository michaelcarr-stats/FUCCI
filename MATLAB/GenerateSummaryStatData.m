function SummaryStatData = GenerateSummaryStatData(CellTracking, SummaryStatData, ntrack, T_record, t, Nred, Nyellow, Ngreen, RowPosCell, ColPosCell, CellSelected, CellSelectedStart, rowIndex, columnIndex, migFailed, transID, delta, simuIndex, simuNum, domain, domain_x)
    
    if CellTracking %if using cell tracking data
        %initialise structure
        if ~isstruct(SummaryStatData)
            %summary statistics
            clear SummaryStatData
            SummaryStatData.Nred_ss = NaN(simuNum,1);
            SummaryStatData.Nyellow_ss = NaN(simuNum,1);
            SummaryStatData.Ngreen_ss = NaN(simuNum,1);
            SummaryStatData.red_distance = zeros(simuNum, ntrack);
            SummaryStatData.yellow_distance = zeros(simuNum, ntrack);
            SummaryStatData.green_distance = zeros(simuNum, ntrack);
            
            %additional helper variables
            SummaryStatData.red_time = zeros(simuNum, ntrack);
            SummaryStatData.yellow_time = zeros(simuNum, ntrack);
            SummaryStatData.green_time = zeros(simuNum, ntrack);
            SummaryStatData.phasestart_time = zeros(simuNum, ntrack);
            SummaryStatData.TrackingCompleted = zeros(simuNum, ntrack, 2);
        end
       
        red_distance = SummaryStatData.red_distance;
        yellow_distance = SummaryStatData.yellow_distance; 
        green_distance = SummaryStatData.green_distance;
        red_time = SummaryStatData.red_time;
        yellow_time = SummaryStatData.yellow_time;
        green_time = SummaryStatData.green_time;
        phasestart_time = SummaryStatData.phasestart_time;
        TrackingCompleted = SummaryStatData.TrackingCompleted;
        
        %record number of cells in simulation
        if (t >= T_record - 1)
            SummaryStatData.Nred_ss(simuIndex) = Nred;
            SummaryStatData.Nyellow_ss(simuIndex) = Nyellow;
            SummaryStatData.Ngreen_ss(simuIndex) = Ngreen;
        end
        
        % record cell trajectory data
        if ~migFailed || transID ~= 0
            if any((RowPosCell == rowIndex).*(ColPosCell== columnIndex))

                %determine tracked cell which moved or transitioned
                 index = (1:ntrack).*(RowPosCell == rowIndex).*(ColPosCell == columnIndex);
                 index = index(index ~= 0);

                if transID ~= 0  %if transition event took place                                        
                    % stop tracking cell if it has transition through
                    % the 3 phases.

                    %record phase duration
                    if TrackingCompleted(simuIndex,index,1) == 0
                        if transID == 1
                            green_time(simuIndex, index) = t - phasestart_time(simuIndex, index); %calculate phase duration
                            if red_time(simuIndex, index) == 0 %ensure not overiding first phase
                                red_time(simuIndex, index) = T_record - t; %assign next phase termination time as 48hrs - to be overriden if completed
                                phasestart_time(simuIndex, index) = t; %reset reference time
                            end
                        elseif transID == 2
                            red_time(simuIndex, index) = t -  phasestart_time(simuIndex, index); %calculate phase duration
                            if yellow_time(simuIndex, index) == 0 %ensure not overiding first phase
                                yellow_time(simuIndex, index) = T_record - t; %assign next phase termination time as 48hrs - to be overriden if completed
                                phasestart_time(simuIndex, index) = t; %reset reference time
                            end
                        elseif transID == 3
                            yellow_time(simuIndex, index) = t -  phasestart_time(simuIndex, index); %calculate phase duration
                            if green_time(simuIndex, index) == 0 %ensure not overiding first phase
                                green_time(simuIndex, index) = T_record - t; %assign next phase termination time as 48hrs - to be overriden if completed
                                phasestart_time(simuIndex, index) = t; %reset reference time
                            end
                        end
                    end


                    if transID == CellSelectedStart(index)% determine if tracked cell has moved through the 3 phases 

                        TrackingCompleted(simuIndex,index,1) = 1;
                        %record time tracking was completed
                        if TrackingCompleted(simuIndex,index,2) == 0                                  
                            TrackingCompleted(simuIndex,index,2) = t;
                        end
                    end  


                elseif ~migFailed  %if movement event took place measure radial distance
                    if TrackingCompleted(simuIndex, index, 1) == 0 % record while cell cycle incomplete
                        if CellSelected(index) == 1 %if cell is red
                            red_distance(simuIndex,index) = red_distance(simuIndex,index) + delta;
                        elseif CellSelected(index) == 2 %if cell is yellow
                            yellow_distance(simuIndex,index) = yellow_distance(simuIndex,index) + delta;
                        elseif CellSelected(index) == 3 %if cell is green
                            green_distance(simuIndex,index) = green_distance(simuIndex,index) + delta;
                        end 
                    end
                end
            end    
        end
        
        % store summary statistics
        SummaryStatData.red_distance = red_distance;
        SummaryStatData.yellow_distance = yellow_distance;
        SummaryStatData.green_distance = green_distance;

        %save additional helper variables
        SummaryStatData.red_time = red_time;
        SummaryStatData.yellow_time = yellow_time;
        SummaryStatData.green_time = green_time;
        SummaryStatData.phasestart_time = phasestart_time;
        SummaryStatData.TrackingCompleted = TrackingCompleted;
        
    else % if cell density was used
        
        if (t >= T_record - 1)
            %initialise structure
            if ~isstruct(SummaryStatData)
                %summary statistics
                clear SummaryStatData
                SummaryStatData.Nred_ss = NaN(simuNum,1);
                SummaryStatData.Nyellow_ss = NaN(simuNum,1);
                SummaryStatData.Ngreen_ss = NaN(simuNum,1);
                SummaryStatData.red_position = cell(simuNum, 1);
                SummaryStatData.yellow_position = cell(simuNum, 1);
                SummaryStatData.green_position = cell(simuNum, 1);               
            end
            
            %compute number of cells in simulation 
            SummaryStatData.Nred_ss(simuIndex) = Nred;
            SummaryStatData.Nyellow_ss(simuIndex) = Nyellow;
            SummaryStatData.Ngreen_ss(simuIndex) = Ngreen;
            
            %compute x position of cells in simulation
            SummaryStatData.red_position(simuIndex) = mat2cell(domain_x(domain == 1), length(domain_x(domain == 1)));
            SummaryStatData.yellow_position(simuIndex) = mat2cell(domain_x(domain == 2), length(domain_x(domain == 2)));
            SummaryStatData.green_position(simuIndex) = mat2cell(domain_x(domain == 3), length(domain_x(domain == 3)));
        end
        
    end
    
    

end
