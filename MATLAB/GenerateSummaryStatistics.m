function sx = GenerateSummaryStatistics(SummaryStatData, CellTracking, Xmax)  
    
    Nred = mean(SummaryStatData.Nred_ss);
    Nyellow = mean(SummaryStatData.Nyellow_ss);
    Ngreen = mean(SummaryStatData.Ngreen_ss);


    if CellTracking
        red_distance = SummaryStatData.red_distance;
        yellow_distance = SummaryStatData.yellow_distance;
        green_distance = SummaryStatData.green_distance;
        red_time = SummaryStatData.red_time;
        yellow_time = SummaryStatData.yellow_time;
        green_time = SummaryStatData.green_time;

        %red_distance(red_time == 0) = NaN; %by default red is always recorded
        yellow_distance(yellow_time == 0) = NaN;
        green_distance(green_time == 0) = NaN;

        red_distance_avg = mean(red_distance,'all','omitnan');
        yellow_distance_avg = mean(yellow_distance,'all','omitnan');
        green_distance_avg = mean(green_distance,'all','omitnan');

        sx = [Nred, Nyellow, Ngreen, red_distance_avg, yellow_distance_avg, green_distance_avg];
    else
        red_position = SummaryStatData.red_position;
        yellow_position = SummaryStatData.yellow_position;
        green_position = SummaryStatData.green_position;
        
        red_distance_med = [median(red_position{1}(red_position{1} <= Xmax/2)),median(red_position{1}(red_position{1} > Xmax/2))];
        yellow_distance_med = [median(yellow_position{1}(yellow_position{1} <= Xmax/2)),median(yellow_position{1}(yellow_position{1} > Xmax/2))]; 
        green_distance_med = [median(green_position{1}(green_position{1} <= Xmax/2)),median(green_position{1}(green_position{1} > Xmax/2))]; 
        
        red_distance_iqr = [iqr(red_position{1}(red_position{1} <= Xmax/2)),iqr(red_position{1}(red_position{1} > Xmax/2))];
        yellow_distance_iqr = [iqr(yellow_position{1}(yellow_position{1} <= Xmax/2)),iqr(yellow_position{1}(yellow_position{1} > Xmax/2))]; 
        green_distance_iqr = [iqr(green_position{1}(green_position{1} <= Xmax/2)),iqr(green_position{1}(green_position{1} > Xmax/2))]; 
        
        sx = [Nred, Nyellow, Ngreen, red_distance_med, yellow_distance_med, green_distance_med, red_distance_iqr, yellow_distance_iqr, green_distance_iqr];
        
    end
end