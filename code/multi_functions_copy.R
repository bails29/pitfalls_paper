#function for smoothing flow
ma <- function(x, n){stats::filter(x, rep(1 / n, n), sides = 1)}

#function to calculate rolling percentile. x = day of the year, n = number of days on either side of x, drought_percentile = threshold desired
all_roll <- function(x, n, percentile, data, attribute){ # x is day of year, n is number of days on each side
  tryCatch({
    data <- data.table(data)
    index <- which(data$day == x)
    index <- index[which(index > n)]
    before <- index-n
    after <- index+n
    
    index2 <- integer(0)
    
    for(t in seq_along(before)) {
      index2 <- c(index2, before[t]:after[t])
    }
    
    index <- index2
    
    return(data[index, quantile(get(attribute), percentile, na.rm = TRUE)][[1]])},
    error = function(e) {
      message("An Error Occurred")
      print(paste(e, "for day ", x))
    }
  )}


########################################
# test GPT function
# Function to calculate percentile for each day based on the previous and following 15 days
rolling_value_percentile <- function(x, n, data, attribute, rank_name = "rank", dir = 1){ # x is day of year, n is number of days on each side
  tryCatch({
    data <- data.table(data)
    index <- which(data$day == x)
    index <- index[which(index > n)] #next step would fail if index is <=15
    before <- index-n
    after <- index+n
    
    index2 <- integer(0)
    
    for (t in seq_along(before)) {
      index2 <- c(index2, before[t]:after[t])
    }
    
    window_data <- data[index2, get(attribute)]
    
    # Locate the central observations within the window.
    indices <- match(index, index2)
    
    # Calculate ranks for those observations.
    test <- sapply(indices, rank_ba, window_data, dir = dir) 
    
    data2 <- data[index]
    
    return(data2[,  (rank_name) := test])},
    error = function(e) {
      message("An Error Occurred")
      print(paste(e, "for day ", x))
    }
  )}




# Function for getting percentiles of each value 
# "dir" allows you to flip it for you to flip the distribution for when you are
# estimating the percentiles of the CDPM series vs the regular time series. 1 is 
# for regular time series, 0 is for CDPM time series. 
rank_ba <- function(i, value, dir){ if (dir == 1){rank <-1 - mean(as.numeric(value) > value[i], na.rm = T)
}else if (dir == 0){rank <- mean(as.numeric(value) > value[i], na.rm = T)}
  return(rank)}


### try to speed up the detection function (currently annoyingly slow)#
#Function to estimate drought and non-drought periods. 
# includes options for restricting drought period length (min_length) 
# adjusting smoothing parameters (smoothing_param)
# adjusting temporal aggregation period for estimation of moving threshold (n)
# for merging events which are close together (min_recovery) e.g. if they are less than n distance, they're merged 
# processed_events_path = where to store the indices files 
# rapid_threshold is the max distance between events to call it a "rapid transition"
# and seasonal_threshold is the max gap for considering a transition "seasonal"
# flood_definition is whether you want a fixed threshold or a variable one for floods
# both thresholds will be returned, but only one will be used to determine if an event is a flood.
# min_lag_flood merge flood events closer than x days together to ensure independence

# ********edit so that you can choose the method.********* 
# ff - both thresholds are fixed, vf - drought is variable, flood is fixed, 
# sf - drought is seasonal, flood is fixed, i - indices approach, 
# fv - fixed drought, variable flood, vv both variable, ss = both seasonal

# snow melt method: D = bi directional degree day split, Bb = bi-direcitonal blind (winter, summer) split,
# B = all blind (winter, spring, summer, fall). 
#NOTE: if rolling average temperature is never below 0, D will just be a fixed threshold

#approach is either fixed, raw, or monthly. 

getSeasonal_p <- function(df, Dmethod = "s", Fmethod = "f", approach = "index",
                          drought_percentile = 0.20, flood_percentile,
                          smoothing_param = 30, min_length = 30,
                          n = 14, min_recovery = 30, processed_events_path = "./events_p/", 
                          rapid_threshold = 14, seasonal_threshold = 90, 
                          snowmelt_method = "D", var_name = "streamflow", est_var = "smoothed",  calc_thresholdD = "threshold",
                          calc_thresholdF = "flood_threshold"){ # df <- flow;  
  # df = flow;  Dmethod = "f"; Fmethod = "f";
  # drought_percentile = 0.30; flood_percentile = 0.70;
  # smoothing_param = 28; min_length = 30;
  # n = 14; min_recovery = 15;
  # processed_events_path = "./test/"; var_name = 'streamflow';
  # snowmelt_method = "D"; approach = "index"
  
  # if(approach == "index"){
  #   if(Dmethod == "s" | Dmethod == "v" | Fmethod == "s" | Fmethod == "v"){
  #   stop("variable and seasonal threhsolds not allowed for indices approaches. Please use a fixed threshold.")
  # }}
  
  #######################
  #section 1: data prep #
  #######################
  #extract ID so you can save files later
  #if(var_name == "Q"){df <- rename(df, Q = streamflow )}
  
  ID <- unique(df$ID)
  
  
  #Calculate the smoothed time series. 
  df[, smoothed := ma(x = get(var_name), n = smoothing_param)]
  df <- df[!is.na(df$smoothed),]
  df$smoothed <- as.numeric(df$smoothed)
  
  #create separate dataframe so you can edit it w/o problems  
  positive_df <- copy(df)
  
  
  
  ##################################################################
  # Drought Threshold Estimation                                   #
  ##################################################################
  if (approach == "index"){
    positive_df$threshold <- drought_percentile
    columns <- c("date", "ID", "threshold", "condition", "rank")  
  }else if(Dmethod == "v"){
    # calculate the variable drought threshold using the TLM methods 
    
    #get the exact day of the year
    positive_df$day <- lubridate::yday(positive_df$date)
    
    #get TLM of all values and select drought periods
    TLM <-  data.table(day = sort(unique(positive_df$day)), threshold = sapply(sort(unique(positive_df$day)),
                                                                               all_roll, n = n, percentile = drought_percentile, 
                                                                               data = positive_df, attribute = est_var))
    
    
    positive_df <- lapply(sort(unique(positive_df$day)), FUN = rolling_value_percentile, n = n,  
                   data = positive_df, attribute = est_var , rank_name = "rank")
    
    positive_df <- rbindlist(positive_df)
    positive_df <- positive_df[order(date)]
    
    positive_df <- merge(TLM, positive_df, by = "day")
    
    columns <- c("date", "ID", "threshold", "condition", "rank", "day")      
  }else if(Dmethod == "f"){
    # calculate the fixed drought threshold                         
    
    positive_df$rank <- positive_df[,sapply(1:length(get(est_var)), rank_ba, get(est_var), dir = 1)]
    
    positive_df$threshold <- positive_df[which.min(abs(drought_percentile-rank))]$smoothed 
    
    columns <- c("date", "ID", "threshold", "condition", "rank", "day")
    
  }else if(Dmethod == "s"){
    # calculate the seasonal drought threshold                     
    
    #type of seasonality is determined by "snowmelt_method" as defined in the forward. 
    if(snowmelt_method == "Bb"){
      #positive_df[, season := NULL]
      positive_df[month %in%  c("09", "10", "11", "12", "01", "02"), season:= "1"]
      positive_df[month %in% c("03", "04", "05", "06", "07", "08"), season:= "2"]
      
    }else if(snowmelt_method == "B"){
      # positive_df <- positive_df[, season := NULL]
      positive_df[month %in%  c("09", "10", "11"), season:= "1"]
      positive_df[month %in% c("03", "04", "05"), season:= "3"]
      positive_df[month %in%  c("12", "01", "02"), season:= "2"]
      positive_df[month %in% c("06", "07", "08"), season:= "4"]
      
    }else if(snowmelt_method == "D"){
      positive_df$day <- lubridate::yday(positive_df$date) # get exact day of year
      
      # calculate variable threshold of temperature using TLM - this is used to get the 28 day rolling median for the doy, 
      #which is then used to determine the average day of melting and freezing. 
      test <- data.table(day = sort(unique(positive_df$day)), threshold = sapply(sort(unique(positive_df$day)), 
                                                                                 all_roll, n = 14, percentile = 0.5,
                                                                                 data = positive_df, 
                                                                                 attribute = "temperature_2m_mean"))
      
      #first doy when the average is over 0 is the melt_day and vice versa for freeze day
      melt <- test %>% 
        summarise(melt_day = first(day[threshold > 0]), 
                  freeze_day = last(day[threshold > 0]))
      
      if(melt$freeze_day - melt$melt_day < 336){
        positive_df[,season := fifelse(day >= melt$freeze_day | day <= melt$melt_day, "1", "2")]
      }else(positive_df[,season := "1"])
      
      rm(test)
      
    }
    # define a fixed threshold per season
    # positive_df <- positive_df %>% 
    #   group_by(season) %>% 
    #   mutate(threshold = quantile(get(var_name), drought_percentile))
    # positive_df <- data.table(positive_df)
    
    # define a fixed threshold per season using smoothed streamflow
    positive_df <- positive_df %>%
      group_by(season) %>%
      mutate(threshold = quantile(get(est_var), drought_percentile))
    positive_df <- data.table(positive_df)
    
    positive_df$rank[positive_df$season == "1"] <- positive_df[positive_df$season == "1",sapply(1:length(get(est_var)), rank_ba, get(est_var), dir = 1)]
    positive_df$rank[positive_df$season == "2"] <- positive_df[positive_df$season == "2",sapply(1:length(get(est_var)), rank_ba, get(est_var), dir = 1)]
    positive_df$rank[positive_df$season == "3"] <- positive_df[positive_df$season == "3",sapply(1:length(get(est_var)), rank_ba, get(est_var), dir = 1)]
    positive_df$rank[positive_df$season == "4"] <- positive_df[positive_df$season == "4",sapply(1:length(get(est_var)), rank_ba, get(est_var), dir = 1)]
    
    columns <- c("date", "ID", "threshold", "condition", "rank", "season", "day")
    
  }
  

  #define as a drought if below threshold, otherwise normal. Reset 0 values to normal??  
  #positive_df[,  condition := ifelse(as.numeric(get(est_var)) <= threshold, "drought", "normal")]
  positive_df[,  condition := ifelse(rank <= drought_percentile, "drought", "normal")] # set condition using percentile converted flow values, instead of variable threshold to facilitate the CDPM method
  
  positive_df[get(var_name) == 0, condition := "normal"] # make 0 values "not drought!"
  
  ############################################
  # Deal with 0 values using the CDPM method #
  ############################################
  #this is equally applied to all drought methods. 
  
  #get associated rank for every value of the smoothed timeseries -- SKIP. WE DID THIS EARLIER IN UPDATE
  #positive_df$rank <- positive_df[,sapply(1:length(get(est_var)), rank_ba, get(est_var), dir = 1)]
  
  #calculate consecutive dry number for all 0 flow periods
  test <- copy(df)
  test <- test[order(date), group := rleid(get(est_var) == 0)] # intermediary step, group 0 values 
  test[order(date), Ndry := cumsum(get(est_var) == 0), by = group] # calculate dry number 
  test$dry_rank <- sapply(1:length(test$Ndry), rank_ba, test$Ndry, dir = 0)  # get associated rank for each dry day 

  
  #combine drought events with 0 flow events 
  positive_df$day <- lubridate::yday(positive_df$date)
  drought_series <- positive_df[,..columns]; rm(positive_df)
  setorder(drought_series, date)
  drought_series <- merge(drought_series, test, by = c("ID", "date"), all = T)
  drought_series[order(date), group := rleid(get(est_var) == 0 | condition == "drought")] # intermediary step, group 0 and drought values 
  drought_series[, Ndd := cumsum(get(est_var) == 0 | condition == "drought"), by = group]

  # drought_series$rank_Ndd <- sapply(1:length(drought_series$Ndd), rank_ba, drought_series$Ndd, dir = 0)  # get associated percentiles for each dry day 

  
  # Interpolation to assign percentiles to Ndd
  if (any(drought_series$smoothed == 0)){
    drought_series[, rank_ndd := ifelse(
    Ndd <= max(Ndry),                   # If Ndd falls within the range of Ndry
    approx(Ndry, dry_rank, Ndd, rule = 2)$y,  # Interpolate dry_rank for Ndd
    1                                 # Assign max percentile if Ndd exceeds Ndry range
  )]
  
    #calculate the scaled percentiles
    # Scale rank_ndd and cap at 100
    drought_series$day <- lubridate::yday(drought_series$date)
    if(Dmethod == "v"){
      drought_series[, scaled_rank_ndd := (min(rank, na.rm = T)/rank_ndd)*rank_ndd, by = day]
    #change rank to dry rank when est_var == 0
    drought_series[get(est_var) == 0, rank := scaled_rank_ndd ]
    
    drought_series[,':='(scaled_rank_ndd = NULL, rank_ndd = NULL)]
    }else {#change rank to dry rank when est_var == 0
      drought_series[get(est_var) == 0, rank := rank_ndd ]
      drought_series[,':='(rank_ndd = NULL)]
      }
    
    
    
  }

  # check if the new series crosses the threshold in order to
  # classify new series as drought or normal
  drought_series[,  condition := ifelse(rank <= drought_percentile, "drought", "normal")]

  
  #Overwite DF and the calculate drought event duration 
  df <- drought_series; rm(drought_series); rm(test)
  
  
  ##########################################################
  #     flood threshold setting. :)                        #
  ##########################################################
  #change var_name for monthly approach
  if(approach == "monthly"){
    var_name <- "monthly_max_flow"
  }  
  
  if(approach == "index"){
    df$flood_threshold <- flood_percentile
  }else if (Fmethod == "v"){
  ########################################################
  #    estimate flood periods from the P time series.    #
  ########################################################
    
    
    
    #get TLM of all values and select flood periods
    df$day <- lubridate::yday(df$date)
    TLM <-  data.table(day = sort(unique(df$day)), 
                       flood_threshold = sapply(sort(unique(df$day)),
                                                all_roll, n = n, percentile = flood_percentile,
                                                data = df,  attribute = var_name))
    
    df <- merge(TLM, df, by = "day")
    
  } else if(Fmethod == "f"){
    # max_series <- df %>% mutate(year = lubridate::year(date)) %>% group_by(year) %>% filter(n_distinct(date) > 336) %>%
    #   group_by(year) %>% summarise(max_var = max(get(var_name)))
    df$flood_threshold <- quantile(df[,get(var_name)], flood_percentile)
  } else if(Fmethod == "s"){
    # calculate the seasonal flood threshold                     
    
    #type of seasonality is determined by "snowmelt_method" as defined in the forward. 
    if(snowmelt_method == "Bb"){
      #df[, season := NULL]
      df[month %in%  c("09", "10", "11", "12", "01", "02"), season:= "1"]
      df[month %in% c("03", "04", "05", "06", "07", "08"), season:= "2"]
      
    }else if(snowmelt_method == "B"){
      #df[, season := NULL]
      df[month %in%  c("09", "10", "11"), season:= "1"]
      df[month %in% c("03", "04", "05"), season:= "3"]
      df[month %in%  c("12", "01", "02"), season:= "2"]
      df[month %in% c("06", "07", "08"), season:= "4"]
      
    }else if(snowmelt_method == "D"){
      df$day <- lubridate::yday(df$date) # get exact day of year
      
      # calculate variable threshold of temperature using TLM - this is used to get the 28 day rolling median for the doy, 
      #which is then used to determine the average day of melting and freezing. 
      test <- data.table(day = sort(unique(df$day)), threshold = sapply(sort(unique(df$day)), 
                                                                        all_roll, n = 14, percentile = 0.5,
                                                                        data = df, 
                                                                        attribute = "temperature_2m_mean"))
      
      #first doy when the average is over 0 is the melt_day and vice versa for freeze day
      melt <- test %>% 
        summarise(melt_day = first(day[threshold > 0]), 
                  freeze_day = last(day[threshold > 0]))
      
      df[,season := fifelse(day >= melt$freeze_day | day <= melt$melt_day, "1", "2")]
      rm(test)
      
    }
    # define a fixed threshold per season
    # max_series <- df %>% mutate(year = lubridate::year(date)) %>% group_by(year) %>% filter(n_distinct(date) > 336) %>%
    #   group_by(year, season) %>% summarise(max_var = max(get(var_name))) %>% 
    #   ungroup() %>%
    #   group_by(season) %>% 
    #   summarise(flood_threshold = quantile(max_var, flood_percentile))
    
    df <- df %>% 
      group_by(season) %>% 
      mutate(flood_threshold = quantile(get(var_name), flood_percentile, na.rm = T))
    #df <- merge(df, max_series, by = "season")
  }
  
  #determine if flood or not
  df <- as.data.table(df)
  df[get(var_name) >= flood_threshold, condition := "flood"]
  
  ########################################################
  # APPLIED TO ALL METHODS                               #
  ########################################################
  
  ####################################
  # Calculate event characteristics  #
  ####################################
  
  #calculate length of event
  df[order(date), group := rleid(condition)]
  df[order(date), first_date := date[1], by = group] # Get the date of the first instance in each group
  df[order(date),  last_date := date[.N], by = group]# Get the date of the last instance in each group
  df[, event_duration := .N, by = group] #get length of each event
  
  
  ###### ****CHECK THIS SECTION!!
  #calculate delay time between drought events, to determine if they should be merged
  test <- df %>% 
    filter(condition == "drought") %>% 
    mutate(group2 = rleid(group))
  
  #add secondary grouping back to original data
  df <- merge(test[,c("group2", "date")], df, by = c("date"), all = T ) #2 is date and 16 is group2
  
  # #merge events if they are close together, determined by threshold defined in the function. 
  df[, merged_condition := condition]
  
  #if the recovery time is too short, change merged condition, dates, and duration
  if(length(unique(df$group2))==1){
    stop(paste0("no drought periods", ID, Dmethod, approach))
  }
  
  for (i in 1:(max(df$group2,  na.rm = T)-1)){
    con <- unique(df[group2 == i]$merged_condition)
    time <- as.numeric(unique(df[group2 == (i+1), first_date])) - as.numeric(unique(df[group2 == i, last_date]))
    if(length(time > 0) ){
      if (time < min_recovery & con == "drought"){
        df[date >= unique(df[group2 == i, first_date]) & date <= unique(df[group2 == (i+1), last_date]), merged_condition := "drought"] # change condition
        df[date >= unique(df[group2 == i, first_date]) & date <= unique(df[group2 == (i+1), last_date]), #change last date
           last_date := unique(df[group2 == (i+1), last_date])]
        df[date >= unique(df[group2 == i, first_date]) & date <= unique(df[group2 == (i+1), last_date]), #change first date
           first_date := unique(df[group2 == i, first_date])]
      }
    }
  }
  
  
  #recalculate length of events now that they have been merged, changed the merged condition.
  #calculate length of event
  df[order(date), group := rleid(merged_condition)]
  df[order(date), first_date := date[1], by = group] # Get the date of the first instance in each group
  df[order(date),  last_date := date[.N], by = group]# Get the date of the last instance in each group
  df[, event_duration := .N, by = group] #get length of each event
  
  #if event is too short, based on predefined threshold, change to normal condition
  df[event_duration < min_length & merged_condition == "drought", merged_condition := "normal"]
  
  ##########################################################
  # get max and cumulative surplus for the flood periods   #
  # get max and cumulative deficit for the drought periods #
  # Negative values are normal or flood conditions where   #
  # we're over the drought threshold and vice versa        #
  ##########################################################
  df <- df %>% 
    group_by(group) %>% 
    mutate(#deficit
      cum_deficit = sum(get(calc_thresholdD) - get(var_name)) , 
      max_deficit =  max(get(calc_thresholdD) - get(var_name)), 
      Qmin = min(get(var_name)), 
      
      #surplus
      cum_surplus = sum(get(var_name) - get(calc_thresholdF)) , 
      max_surplus =  max(get(var_name) - get(calc_thresholdF)), 
      Qpeak = max(get(var_name))) %>% 
    data.table()
  
  ###################################################
  # calculate time between drought and flood events #
  ###################################################
  
  #calculate length of event
  drought_groups <- unique(df$group[df$merged_condition == "drought"])
  flood_groups <- unique(df$group[df$merged_condition == "flood"])
  
  indices <- df[match(drought_groups, df$group),] %>% 
    dplyr::select(group, last_date, first_date, event_duration, cum_deficit, max_deficit, Qmin) %>% 
    mutate(class = "drought")
  
  flood_indices <- df[match(flood_groups, df$group),] %>% 
    dplyr::select(group, last_date, first_date, event_duration, cum_surplus, max_surplus, Qpeak) %>% 
    mutate(class = "flood")
  names(flood_indices) <- names(indices)
  
  indices <- rbind(indices, flood_indices)
  
  indices <- indices[order(first_date),]
  
  indices$group <- as.numeric(row.names(indices))
  
  # Add transition column and time between events
  indices <- indices %>%
    mutate(time_since_previous_event = as.numeric(indices$first_date) - dplyr::lag(as.numeric(indices$last_date)),
           time_since_previous_event = as.numeric(time_since_previous_event), 
           transition = case_when(
             class != dplyr::lag(class) & time_since_previous_event <= seasonal_threshold ~ "transition", #| class != lead(class)
             .default = "independent"
           )) %>% 
    mutate(transition_type = case_when(transition == "transition" & time_since_previous_event <= rapid_threshold ~ "rapid", 
                                       transition == "transition" & time_since_previous_event > rapid_threshold &  time_since_previous_event <= seasonal_threshold  ~ "seasonal", 
                                       #transition == "transition" & time_since_previous_event > seasonal_threshold ~ "long", 
                                       #transition == "independent" & lead(transition) == "transition" & lead(time_since_previous_event) <= seasonal_threshold ~ lead(transition_type),
                                       .default = "independent"), 
           transition_type = case_when(transition == "independent" & dplyr::lead(transition) == "transition" & dplyr::lead(time_since_previous_event) <= seasonal_threshold ~ dplyr::lead(transition_type),
                                       TRUE ~ transition_type))
  
  #make groups
  # Define groups based on conditions
  indices$event_group <- NA
  indices[, event_group := ifelse(transition_type %in% c("seasonal", "rapid"), .I, 0)]
  
  test <- indices[event_group != 0,]
  indices$event_group[indices$event_group != 0] <- cumsum(c(1, diff(test$group) != 1 )) 
  indices <- indices %>% 
    mutate(group = NULL)
  
  # Return event group characteristics
  indices <- indices %>% distinct()
  indices$ID <- ID
  indices$approach <- approach
  indices$Dmethod <- Dmethod
  indices$Fmethod <- Fmethod
  indices$nyears <- length(unique(df$year))
  
  #save the indices 
  fwrite(indices, paste0(processed_events_path, ID, "_", approach, "new_calc_event_indices_", Fmethod, "_", Dmethod, "_", snowmelt_method, ".csv"))
  
  #save the complete timeseries
  df <- merge(df, indices[, c("max_deficit", "Qmin", "cum_deficit") := NULL], by = c("first_date", "last_date", "event_duration", "ID"), all = T)
  
  
  #THIS IS THE WARNING PROBLEM CHILD!! FIX THIS!!
  for (i in 1:max(unique(df$event_group), na.rm = TRUE)){
    df$event_group[df$date >= min(unique(df$first_date[df$event_group == i]), na.rm = TRUE) & df$date <= max(unique(df$last_date[df$event_group == i]), na.rm = TRUE)] <- i
  }
  
  
  df <- df %>% 
    mutate(group = NULL, 
           Ndry = NULL, 
           CDPM_rank = NULL, 
           #rank = NULL, 
           Ndd = NULL, 
           class = NULL, 
           #transition = NULL, 
           group2 = NULL, 
           day = NULL, 
           dry_rank = NULL, 
           rank_Ndd = NULL, 
           Qpeak = NULL, 
           max_surplus = NULL,
           cum_surplus = NULL)
  
  
  #
  #save the complete data 
  #fwrite(df, paste0(processed_events_path, ID, Dmethod, Fmethod, snowmelt_method,"_all.csv"))
  return(df)
}


