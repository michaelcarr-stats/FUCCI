function s = SetupStruct(ntrack, Xmax, Ymax, InitPosData, CellTrackingData)

%% Input Values
    simuNum = 3; %number of simulation realisations
    delta = 20; %Lattice spacing [mum]
    t = 0; %Time clock starting with 0
    tstop = 48; % Total iterations at each realisation
    columnNum = ceil(Xmax/delta) + 1; %Corresponding total column nodes
    rowNum = ceil(Ymax/(delta*sqrt(3)/2)) + 1; %Corresponding total row nodes
    initialL = 20*delta; %Length of the initial agent pack
    initialDensity = 1; %Density of the initial pack

    percentR = 0.33; % percentage of the initial domain to be red
    percentY = 0.33; % percentage of the initial domain to be yellow
    percentG = 0.33; % percentage of the initial domain to be green
     
    BC = 2; %Boundary condition, 1.Periodic, 2. No flux
 
 %% Initialisation
    domain = zeros(rowNum,columnNum); %cell phase identity for nodes: 1 if red, 2 if yellow, 3 if green, 0 if unoccupied
    domain_x = zeros(rowNum,columnNum); %x coordinates for all the nodes
    domain_y = zeros(rowNum,columnNum); %y coordinates for all the nodes
    %Construct a coordinate system in the domain (See Jin et al., 2017)
    for i=1:rowNum
        for j=1:columnNum
            if rem((rowNum-i),2) == 0
                domain_x(i,j) = (j-0.5) * delta; %This ensures the first node at the bottom row always is associated with a x coordinate of 0.5
            else
                domain_x(i,j) = (j-1) * delta; %This ensures the first node at the second bottom row always is associated with a x coordinate of 0
            end
            domain_y(i,j) = (rowNum-i) * sqrt(3) / 2 * delta; %Constructing y coordinates
        end
    end
    
    %Construct the initial domain
    if (isempty(InitPosData) && isempty(CellTrackingData))
        domain = initialise_domain(domain, domain_x, rowNum, columnNum, initialDensity, Xmax, initialL, percentR, percentY);
        [RowPosCell,ColPosCell,CellSelectedIndex] = InitialMotilitySelect(ntrack, domain, domain_x, Xmax);
    else
        i = round(rowNum - 2*InitPosData(:,2)/(delta*sqrt(3)));
        j = round(InitPosData(:,1)/delta + 1);
        Indexs = unique([i, j, InitPosData(:,3)], 'rows');
        
        for row = 1:size(Indexs,1)
            domain(Indexs(row,1), Indexs(row,2)) = Indexs(row,3);
        end
        
        if (isempty(CellTrackingData))
            [RowPosCell,ColPosCell,CellSelectedIndex] = InitialMotilitySelect(ntrack, domain, domain_x, Xmax);
        else 
            %initialise motility parameters to track
            RowPosCell = NaN(ntrack,1);
            ColPosCell = NaN(ntrack,1);
            for i = 1:ntrack
                Initx = CellTrackingData(CellTrackingData(:,5) == 1,1);
                Inity = CellTrackingData(CellTrackingData(:,5) == 1,2);

                tmp = [sqrt((InitPosData(:,1) - Initx(i)).^2+(InitPosData(:,2) - Inity(i)).^2),(1:size(InitPosData,1))'];
                index = tmp(tmp(:,1) == min(tmp(:,1)),2); %compute distance between Inital Position data and move tracking data
                RowPosCell(i) = round(rowNum - 2*InitPosData(index,2)/(delta*sqrt(3))); 
                ColPosCell(i) = round(InitPosData(index,1)/delta +1);
            end

            CellSelectedIndex = NaN(ntrack,1);
            for i = 1:ntrack
                CellSelectedIndex(i) = (ColPosCell(i)-1)*size(domain,1) + RowPosCell(i);
            end
        end
    end
    
    [NstartRed, NstartYellow, NstartGreen] = InitialCellCount(domain, domain_x, domain_y, rowNum, columnNum, Xmax, Ymax);
    CellSelected = domain(CellSelectedIndex);
    
    %assign values to structure
    s.ntrack = ntrack;
    s.simuNum = simuNum;
    s.delta = delta; 
    s.t = t;
    s.tstop = tstop;
    s.Xmax = Xmax;
    s.Ymax = Ymax;
    s.columnNum = columnNum; 
    s.rowNum = rowNum;
    s.BC = BC;
    s.domain_x = domain_x;
    s.domain_y = domain_y;
    s.domain = domain;
    s.NstartRed = NstartRed;
    s.NstartYellow = NstartYellow;
    s.NstartGreen = NstartGreen;
    s.RowPosCell = RowPosCell';
    s.ColPosCell = ColPosCell';
    s.CellSelected = CellSelected';
    
    end