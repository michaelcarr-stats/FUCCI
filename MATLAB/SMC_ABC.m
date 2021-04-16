function SMC_ABC
    %% parrallel computing setup
    ncpus = 16; %number of CPUs to use
    clust = parcluster('local');
    clust.JobStorageLocation = tempdir;
    par = parpool(clust,ncpus); 

    rng(1234); % set seed
    
    %% Inputs
    N = 1000; %initial number of simulations
    T_record = 48; %time period to record summary statistics

    tol_target = 0;  % TARGET TOLERANCE.  IF SET TO ZERO, THE STOPPING RULE IS BASED ON MCMC ACCEPTANCE RATE BELOW
    pacc_target = 0.01; % STOP ALGORITHM IF MCMC ACCEPTANCE RATE IS TOO SMALL.  IF SET TO 0, STOPPING RULE IS BASED ON ABOVE TARGET TOLERANCE

    %tuning parameters
    alpha = 0.5;
    c = 0.01;
    
    %initial R_t value
    MCMC_trials = 30;% 
    
    %transition prior - uniform U(t_lb,t_ub)
    t_lb = 0; t_ub = 1;
   
    %Movement prior - uniform U(m_lb,m_ub)
    m_lb = 0; m_ub = 10;

    num_params = 6; %number of parameters to estimate
    
    %% Observed Data
    CellTracking = true; %TRUE - if using cell trajectory data; FALSE - if using cell density data
    
    %load data
    InitPosData = readmatrix("FUCCI_proccessed.xlsx", "sheet","InitPos");
    FinalPosData = readmatrix("FUCCI_proccessed.xlsx", "sheet","FinalPos");
    CellTrackingData = readmatrix("FUCCI_proccessed.xlsx", "sheet","CellTracking");
    ntrack = max(CellTrackingData(:,4));
    
    Xmax = 1309.09; %Length of the domain
    Ymax = 1745.35; %Width of the domain
    
    %summary statistics: 
    Nred = sum(FinalPosData(:,3) == 1);
    Nyellow = sum(FinalPosData(:,3) == 2);
    Ngreen = sum(FinalPosData(:,3) == 3);
    
    if CellTracking %if using cell tracking summary statistics
        
        %compute distance travelled through each cell phase by all tracked cells
        RedDistance = 0;
        YellowDistance = 0;
        GreenDistance = 0;
        j = 2;
        for i = 1:ntrack
           while CellTrackingData(j,4) == i
               if CellTrackingData(j,3) == 1 && CellTrackingData(j-1,3) == 1
                   RedDistance = RedDistance + sqrt((CellTrackingData(j,1)-CellTrackingData(j-1,1))^2+(CellTrackingData(j,2)-CellTrackingData(j-1,2))^2);
               elseif CellTrackingData(j,3) == 2 
                   YellowDistance = YellowDistance + sqrt((CellTrackingData(j,1)-CellTrackingData(j-1,1))^2+(CellTrackingData(j,2)-CellTrackingData(j-1,2))^2);
               elseif CellTrackingData(j,3) == 3 
                   GreenDistance = GreenDistance + sqrt((CellTrackingData(j,1)-CellTrackingData(j-1,1))^2+(CellTrackingData(j,2)-CellTrackingData(j-1,2))^2);
               end
               j = j + 1;
               if (j > size(CellTrackingData,1))
                   break
               end
           end
        end

       TotalDistance = [RedDistance/ntrack, YellowDistance/ntrack, GreenDistance/ntrack];
       
    else %if using cell density summary statistics
        
        %compute median positon of cells on LHS and RHS of scratch
        RedDistanceMed = [median(FinalPosData(FinalPosData(:,3) == 1 & FinalPosData(:,1) <= Xmax/2, 1)), median(FinalPosData(FinalPosData(:,3) == 1 & FinalPosData(:,1) > Xmax/2, 1))];
        YellowDistanceMed = [median(FinalPosData(FinalPosData(:,3) == 2 & FinalPosData(:,1) <= Xmax/2, 1)), median(FinalPosData(FinalPosData(:,3) == 2 & FinalPosData(:,1) > Xmax/2, 1))];
        GreenDistanceMed = [median(FinalPosData(FinalPosData(:,3) == 3 & FinalPosData(:,1) <= Xmax/2, 1)), median(FinalPosData(FinalPosData(:,3) == 3 & FinalPosData(:,1) > Xmax/2, 1))];
        
        %compute IQR of cells on the LHS and RHS of scratch
        RedDistanceIQR = [iqr(FinalPosData(FinalPosData(:,3) == 1 & FinalPosData(:,1) <= Xmax/2, 1)), iqr(FinalPosData(FinalPosData(:,3) == 1 & FinalPosData(:,1) > Xmax/2, 1))];
        YellowDistanceIQR = [iqr(FinalPosData(FinalPosData(:,3) == 2 & FinalPosData(:,1) <= Xmax/2, 1)), iqr(FinalPosData(FinalPosData(:,3) == 2 & FinalPosData(:,1) > Xmax/2, 1))];
        GreenDistanceIQR = [iqr(FinalPosData(FinalPosData(:,3) == 3 & FinalPosData(:,1) <= Xmax/2, 1)), iqr(FinalPosData(FinalPosData(:,3) == 3 & FinalPosData(:,1) > Xmax/2, 1))];
        
        TotalDistance = [RedDistanceMed,YellowDistanceMed,GreenDistanceMed,RedDistanceIQR,YellowDistanceIQR,GreenDistanceIQR]; 
        
    end
    
    sy = [Nred, Nyellow, Ngreen, TotalDistance]; %store observed summary statistics
    
    %% Simulation environment setup
    s = SetupStruct(ntrack, Xmax, Ymax, InitPosData, CellTrackingData);

    %% SMC ABC Algorithm - proliferation and motility
    
    %initialisation
    theta = zeros(N,num_params+1);
    summaries = zeros(N,length(sy));
    
    % creating initial sample
    parfor (i = 1:N,ncpus)
        while ~all(theta(i,:)) %repeat untill 1 proposed value is accepted
                % Draw theta ~ pi(theta)
                theta_prop = [unifrnd(t_lb,t_ub,1,3),unifrnd(m_lb,m_ub,1,3)];

                %simulate x ~ f(x|theta)
                [SummaryStatData, ExitSimStatus] = Main_simulate(theta_prop, s, T_record, CellTracking);
                              
                if ExitSimStatus %check if simulation finished correctly
                    continue;
                end 
                
                sx = GenerateSummaryStatistics(SummaryStatData, CellTracking, Xmax); %compute simulated summary stats
                if any(isnan(sx)) %check integrity of summary statistics
                    continue  
                end 
                
                d = norm(sx - sy, 2); %discrepency function
                %store results
                theta(i,:) = [theta_prop, d];
                summaries(i,:) = sx; %store simulated summary statistics
        end
    end

    %transform parameters to be strictly within the priors and combine with
    %summaries & tolerance
    theta_trans = [TransformationFun(theta(:,1:3),1,false),TransformationFun(theta(:,4:6),10,false),theta(:,num_params + 1), summaries];
    
    %update N and N_a - due to rejecting parameters
    N = size(theta_trans, 1);
    N_a = floor(alpha*N);

    %sort theta by the discrepency function
    theta_trans = sortrows(theta_trans, num_params + 1);
    
    %initialise
    p_accept = 1;
    tolmax = theta_trans(N-N_a,num_params+1);
    
    while p_accept > pacc_target && tolmax > tol_target

        %set tol
        tol = theta_trans(N-N_a,num_params+1);

        %resampling
        r = randsample(1:N-N_a, N_a); %row index to be drawn
        theta_trans((N-N_a + 1):N,:) = theta_trans(r,:);

        %compute tuning parameter
        cov_rw = cov(theta_trans(:,1:num_params));

        %seperate discrepency vector from theta for parrallel computing
        d_vec = theta_trans(:, num_params + 1);
        summaries = theta_trans(:, (num_params + 2):end);
        theta_trans = theta_trans(:, 1:num_params);

        accept = zeros(N,1); % initialise acceptence count matrix 
        
        %trial mcmc step
        parfor (j = (N - N_a +1):N,ncpus)

            theta_trans_j = theta_trans(j,:);

            for k = 1:MCMC_trials 
                
                %Genetate theta ~ q(.|theta_prop)
                theta_trans_prop = mvnrnd(theta_trans_j,cov_rw);
                
                %MH ratio - f(theta_trans) = exp(-theta_trans)./(1+exp(-theta_trans)).^2
                r = exp(sum(log(exp(-theta_trans_prop))) - sum(log((1+exp(-theta_trans_prop)).^2)) - (sum(log(exp(-theta_trans_j))) - sum(log((1+exp(-theta_trans_j)).^2))));
                
                if rand > r
                    continue %reject proposal
                end
                
                %transform back to theta to simulate
                theta_prop = [TransformationFun(theta_trans_prop(:,1:3),1,true),TransformationFun(theta_trans_prop(:,4:6),10,true)];
                
                 %simulate x ~ f(x|theta)
                [SummaryStatData, ExitSimStatus] = Main_simulate(theta_prop, s, T_record, CellTracking);
                              
                if ExitSimStatus %check if simulation finished correctly
                    continue;
                end 
                
                sx = GenerateSummaryStatistics(SummaryStatData, CellTracking, Xmax);
                if any(isnan(sx)) %check integrity of summary statistics
                    continue 
                end                 
                
                d = norm(sx - sy, 2); %discrepency function
                if d > tol
                    continue %reject proposal
                end
                                
                theta_trans_j = theta_trans_prop; %update theta_j
                theta_trans(j,:) = theta_trans_j;%overwrite with best theta_j
                d_vec(j,1) = d; %overwrite with new tolerance
                accept(j,1) = accept(j,1) + 1; %count of MH acceptances
                summaries(j,:) = sx %store simulated summary statistics

            end

        end
        
        
        p_accept = sum(accept((N - N_a +1):N))/(N_a*MCMC_trials);%calculate acceptance rate
        MCMC_iters = ceil(log(c)/log(1-p_accept)); %calculate additional mcmc moves required
               
        % Move Step
        parfor (j = (N - N_a +1):N,ncpus)

            theta_trans_j = theta_trans(j,:);

            for k = (MCMC_trials + 1):(MCMC_iters - MCMC_trials)
                %Genetate theta ~ q(.|theta_prop)
                theta_trans_prop = mvnrnd(theta_trans_j,cov_rw);
                
                %MH ratio - f(theta_trans) = exp(-theta_trans)./(1+exp(-theta_trans)).^2
                r = exp(sum(log(exp(-theta_trans_prop))) - sum(log((1+exp(-theta_trans_prop)).^2)) - (sum(log(exp(-theta_trans_j))) - sum(log((1+exp(-theta_trans_j)).^2))));
                
                if rand > r
                    continue %reject proposal
                end
                
                %transform back to theta to simulate
                theta_prop = [TransformationFun(theta_trans_prop(:,1:3),1,true),TransformationFun(theta_trans_prop(:,4:6),10,true)];
                
                 %simulate x ~ f(x|theta)
                [SummaryStatData, ExitSimStatus] = Main_simulate(theta_prop, s, T_record, CellTracking);
                              
                if ExitSimStatus %check if simulation finished correctly
                    continue;
                end 
                
                sx = GenerateSummaryStatistics(SummaryStatData, CellTracking, Xmax);
                if any(isnan(sx)) %check integrity of summary statistics
                    continue 
                end                 
                
                d = norm(sx - sy, 2); %discrepency function
                if d > tol
                    continue %reject proposal
                end
                                
                theta_trans_j = theta_trans_prop; %update theta_j
                theta_trans(j,:) = theta_trans_j;%overwrite with best theta_j
                d_vec(j,1) = d; %overwrite with new tolerance
                accept(j,1) = accept(j,1) + 1; %count of MH acceptances
                summaries(j,:) = sx %store simulated summary statistics
            end

        end
             
        %recombine theta and discrepency vector for sorting
        theta_trans(:,num_params + 1) = d_vec;
        theta_trans(:, (num_params + 2):(num_params + 1 + size(summaries,2))) = summaries;
        
        %sort theta by the discrepency function
        theta_trans = sortrows(theta_trans, num_params + 1);

        %calculate acceptance ratio
        num_MCMC_iters = max(0,MCMC_iters - MCMC_trials) + MCMC_trials; %calculate number of mcmc iterations
        accept = accept((N - N_a + 1):N,1); %remove unsampled rows
        p_accept = sum(accept)/(length(accept)*num_MCMC_iters); %compute overall acceptance rate
        
        MCMC_iters = ceil(log(c)/log(1-p_accept)); %compute number of mcmc iterations required
        MCMC_trials = ceil(MCMC_iters/2); %compute pilot number of mcmc iterations for next sequence

        tolmax = theta_trans(N-N_a,num_params + 1); %compute max tolerance for next sequence
        p_unique = length(unique(theta_trans, 'rows'))/N; %compute number of unique particles
        
        %display performance stats
        display(p_accept)
        display(tolmax)
        display(p_unique)
        
        save(sprintf('progress%d.mat',ntrack),'theta_trans', 'p_accept','tolmax','p_unique'); %write performance stats to excel
    end
    
    % Correct way to close the parallel pool
    delete(par);
    
   %% Regression Adjustment
    
    %initialise
    theta_trans_regressed = zeros(N,num_params);
    
    %seperate discrepency vector from theta for parrallel computing
    d_vec = theta_trans(:, num_params + 1);
    summaries = theta_trans(:, (num_params + 2):end);
    theta_trans = theta_trans(:, 1:num_params);
   
    %weighting
    the_weights = 0.75*(1-(d_vec./max(d_vec)).^2);
    
    %apply adjustment to one parmater at a time
    for i = 1:num_params
        b = glmfit(summaries,theta_trans(:,i),'normal','weights',the_weights); %calculate coefficients
        theta_trans_regressed(:,i) = theta_trans(:,i)-(summaries-repmat(sy,N,1))*b(2:end); %apply adjustment
    end    
   
   %inverse transform particles
   theta = [TransformationFun(theta_trans(:,1:3),1,true), TransformationFun(theta_trans(:,4:6),10,true)];
   theta_regressed = [TransformationFun(theta_trans_regressed(:,1:3),1,true), TransformationFun(theta_trans_regressed(:,4:6),10,true)];
   
%% save results   

    filename = sprintf('results_SMCABC_track%d.mat',ntrack);
    save(filename,'theta','theta_regressed','sy','d_vec', 'summaries'); %save data

    quit;
end  
