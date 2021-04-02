########################################################################################################
#
# Code accompaniment to article in BMC Public Health:
# Comparing Malaria Early Detection Methods in a Declining Transmission Setting in Northwestern Ethiopia
#
# Code Author: Dawn Nekorchuk, dawn.nekorchuk@ou.edu
# Date: 2021-04-02
#
# This file: Functions for comparing various Early Detection algorithms
#   Available:
#     Random alert (based on a weekly probability)
#     CDC EARS (surveillance::earsC()))
#     Weekly statistics (called WHO, but encompasses more weekly statistics)
#     Farrington (surveillance::farringtonFlexible())
#
########################################################################################################

# See each section for more information on the algorithms, or the manuscript

# Input data: 
# In script to run event detection, the input data must be the TWST results WITH
#  observed data pre-pended, otherwise would lose info based on spin up times.
# See ed_comparison_run_public.R

# Assumed fields:
# NewPCODE: ID for geo groups
# obs_date: date of observation

# Note: does NOT confirm correct amount of data exists for b and n_years (Farrington, WHO methods)


#main wrapper function
run_ed_comparison <- function(data_events, 
                              outpath_base = NULL,
                              #number of weeks at end of data to evaluate
                              eval_last_weeks = 52,
                              #allows alarms to be counted if triggered prior to event (up to 3 weeks)
                              pre_week_search = 3, 
                              run_random = FALSE,
                              settings_random = NULL,
                              run_ears = FALSE,
                              settings_ears = NULL,
                              run_who = FALSE, 
                              settings_who = NULL,
                              run_farrington = FALSE,
                              settings_farrington = NULL){
  
  #make output folder
  res_dt <- format(Sys.time(), "%Y%m%d%H%M")
  outpath <- file.path(outpath_base, res_dt) 
  dir.create(outpath, recursive = TRUE, showWarnings = TRUE)
  
  #data prep
  data_dates <- data_events$obs_date %>% unique() %>% sort()
  index_start <- length(data_dates) - eval_last_weeks + 1
  date_start <- data_dates[index_start]
  
  #order is important for dealing with stss that do NOT have ids
  #Hard coded NewPCODE until field passing added
  data_events <- data_events %>% 
    arrange(NewPCODE, obs_date)
  
  geogroups <- unique(data_events$NewPCODE) %>% sort()
  
  
  #metadata
  metadata <- list(res_dt,
                   outpath,
                   eval_last_weeks,
                   pre_week_search,
                   run_random,
                   settings_random,
                   run_ears,
                   settings_ears,
                   run_who, 
                   settings_who,
                   run_farrington,
                   settings_farrington,
                   data_dates,
                   index_start,
                   date_start,
                   geogroups)
  
  #run each algo with each setting
  
  collector_summary <- tibble()
  
  if (run_random){
    ## Random
    rands_summary <- run_random_ed(data_events = data_events, 
                                   outpath = outpath, 
                                   date_start = date_start,
                                   pre_week_search = pre_week_search,
                                   settings_random = settings_random,
                                   res_dt = res_dt,
                                   geogroups = geogroups)
    
    collector_summary <- bind_rows(collector_summary, rands_summary)
  }
  
  if (run_ears){
    ## EARS
    ears_summary <- run_ears_ed(data_events = data_events, 
                                outpath = outpath, 
                                eval_last_weeks = eval_last_weeks,
                                pre_week_search = pre_week_search,
                                settings_ears = settings_ears,
                                res_dt = res_dt,
                                geogroups = geogroups)
    
    collector_summary <- bind_rows(collector_summary, ears_summary)
  }
  
  if(run_who){
    ## WHO
    who_summary <- run_who_ed(data_events = data_events, 
                              outpath = outpath, 
                              eval_last_weeks = eval_last_weeks,
                              pre_week_search = pre_week_search,
                              settings_who = settings_who,
                              res_dt = res_dt,
                              geogroups = geogroups)
    
    collector_summary <- bind_rows(collector_summary, who_summary)
  }
  
  
  if(run_farrington){
    ##Farrington
    farrington_summary <- run_farrington_ed(data_events = data_events, 
                                            outpath = outpath, 
                                            eval_last_weeks = eval_last_weeks,
                                            pre_week_search = pre_week_search,
                                            settings_farrington = settings_farrington,
                                            res_dt = res_dt,
                                            geogroups = geogroups)
    
    collector_summary <- bind_rows(collector_summary, farrington_summary)
  }
  
  
  #add basic meta data to summary
  collector_summary <- collector_summary %>% 
    dplyr::mutate(timestamp = res_dt,
                  eval_last_weeks = eval_last_weeks,
                  pre_week_search = pre_week_search)
  
  write_csv(collector_summary, path = file.path(outpath, "detection_comparisons.csv"))
  saveRDS(metadata, file.path(outpath, "event_detection_run_metadata.RDS"))
  
  #return
  collector_summary
}




#make farrington combinations
#maybe generalize later, farrington important here
make_farrington_combinations <- function(){
  
}


run_random_ed <- function(data_events, 
                          outpath, 
                          date_start,
                          pre_week_search,
                          settings_random,
                          res_dt,
                          geogroups){
  
  
  #set up data
  data_rand <- data_events %>% 
    dplyr::filter(obs_date >= date_start)
  
  #convert into sts for similar formatting to others
  rand_stss <- make_stss_list(geogroups, data_rand, real_dates = TRUE)
  #immediately convert back to df
  rand_dfs <- lapply(rand_stss, as.data.frame)
  
  ## Create random alarms
  
  #collector for set of settings
  rand_res <- tibble()
  
  #run each set of settings
  for (r in seq_along(settings_random)){
    
    this_rand_stss <- rand_stss
    this_rand_dfs <- rand_dfs
    
    #<<>> Add in defaults or stops?
    #no alarm can be within n1 wks of another 
    #alarm runs are capped at n2 wks (no alarms lasting longer than n2 weeks)
    this_p <- settings_random[[r]]$alarm_p 
    this_buffer <- settings_random[[r]]$buffer_wks
    this_max <- settings_random[[r]]$max_alarm_wks 
    
    #create random alarms
    for (w in 1:length(this_rand_dfs)){ #geogroups
      i <- 1
      while (i <= nrow(this_rand_dfs[[w]])){ # weeks, using while because of variable index increments
        #create random number
        rn <- runif(1)
        #compare to prob of alarm
        if (rn < this_p){
          #randomly pick number of weeks (0 to +1 for max b/c of floor and indexing)
          nwks <- floor(runif(1, min = 0, max = this_max + 1))
          #check number of weeks left for geogroup, adjust alarm weeks as needed
          if (i + nwks > nrow(this_rand_dfs[[w]])) {
            a_wks <- nrow(this_rand_dfs[[w]]) - i
          } else a_wks <- nwks 
          for (n in i:(i + a_wks)){
            this_rand_dfs[[w]][n, "alarm"] <- TRUE
          } #end n (alarm weeks) loop
          #check if at end of geogroup, if so, jump out of this level loop
          if (i + a_wks == nrow(this_rand_dfs[[w]])) {break}
          #add buffer weeks
          #check number of weeks left for geogroups, adjust buffer weeks as needed
          if (i + a_wks + this_buffer > nrow(this_rand_dfs[[w]])) {
            buffer <- nrow(this_rand_dfs[[w]]) - (i + a_wks)
          } else {buffer <- this_buffer}
          for (b in (i + a_wks + 1):(i + a_wks + buffer)){
            this_rand_dfs[[w]][b, "alarm"] <- FALSE
          } #end else of if buffer weeks
          #jump week index to after alarm + buffer (and if at end of woreda, break out of this loop level)
          if ((i + a_wks + buffer + 1) <= nrow(this_rand_dfs[[w]])) {
            i <- i + a_wks + buffer + 1
          } else {break}
          #end if alarm
          #else for no alarm  
        } else {
          this_rand_dfs[[w]][i, "alarm"] <- FALSE
          i <- i + 1
        }
      } #end per week loop
    } #end geogroup loop
    
    
    this_rand_overlaps <- calculate_overlaps(this_rand_dfs, pre_week_search)
    
    
    #organize results
    
    this_rand_res <- dplyr::bind_cols(events_caught(this_rand_overlaps[[1]]),
                                      events_timely(this_rand_overlaps[[1]]),
                                      alarms_events(this_rand_overlaps[[2]])) %>% 
      dplyr::mutate(algorithm = "random",
                    timestamp = res_dt,
                    start_data = min(data_events$obs_date),
                    start_eval = min(data_rand$obs_date),
                    end_eval = max(data_rand$obs_date),
                    #input_settings = list(p_perweek = this_p, max_alarm_wks = this_max, buffer_wks = this_buffer),
                    p_perweek = this_p,
                    max_alarm_wks = this_max,
                    buffer_wks = this_buffer) %>% 
      #order
      dplyr::select(algorithm, everything())
    
    #collect summary results
    rand_res <- dplyr::bind_rows(rand_res, this_rand_res)
    
    
    #save out graphs
    #for graphing, have to add alarms into sts
    #note: requires stss & dfs to remain in the same order for data matching and subtitling to work
    
    for (wor in 1:length(this_rand_dfs)){
      this_rand_stss[[wor]]@alarm <- as.matrix(this_rand_dfs[[wor]][,"alarm"])
    }
    graph_file <- file.path(outpath, paste0("random_alarms_graphs_p", this_p, 
                                            "_buff", this_buffer, 
                                            "_max", this_max, ".pdf"))
    pdf(graph_file, width = 11, height = 8.5)
    for (i in 1:length(this_rand_stss)){
      plot(this_rand_stss[[i]], main = "Random Alarms for Skill Test", sub = geogroups[i],
           legend.opts = list(horiz = TRUE, x = "topright"))
    }
    dev.off()
    
    #save out full data results
    #add names of geogroups
    this_rand_out <- setNames(this_rand_dfs, geogroups)
    saveRDS(this_rand_out, file = file.path(outpath, paste0("random_raw_results_p", this_p, 
                                                            "_buff", this_buffer, 
                                                            "_max", this_max, ".RDS")))
    
  } #end each settings
  
  #final collected results
  rand_res
  
} #end random


#EARS function
run_ears_ed <- function(data_events, 
                        outpath, 
                        eval_last_weeks,
                        pre_week_search,
                        settings_ears,
                        res_dt,
                        geogroups){
  
  # data setup 
  max_date <- max(data_events$obs_date)
  
  #collector for set of settings
  ears_res <- tibble()
  
  #run each set of settings
  for (e in seq_along(settings_ears)){
    
    #restrict time ranges: each have different spin-up/burn-in/training times
    #training time (defaults: C1 = 7 weeks, C2 = 9 weeks, C3 = 11 weeks)
    #settings_ears:
    # method = c(C1, C2, C3)
    # baseline
    # alpha
    #By default if alpha is NULL the value 0.001 is assumed for C1 and C2 whereas 0.025 is assumed for C3. These different choices are the one made at the CDC.
    
    #settings
    if (!is.null(settings_ears[[e]][["method"]])){
      this_method <- settings_ears[[e]][["method"]]
    } else {
      this_method <- "C1"
    }
    if (!is.null(settings_ears[[e]][["alpha"]])){
      this_alpha <- settings_ears[[e]][["alpha"]]
    } else {
      this_alpha <- dplyr::case_when(
        this_method == "C1" ~ 0.001, 
        this_method == "C2" ~ 0.001,
        this_method == "C3" ~ 0.025,
        TRUE ~ NA_real_
      ) 
    }
    if (!is.null(settings_ears[[e]][["baseline"]])){
      this_baseline <- settings_ears[[e]][["baseline"]]
    } else {
      this_baseline <- 7
    }
    
    #data set up for consistent evaluation period with different training times
    this_training <- dplyr::case_when(
      this_method == "C1" ~ this_baseline + 1, 
      this_method == "C2" ~ this_baseline + 3,
      this_method == "C3" ~ this_baseline + 5,
      TRUE ~ NA_real_
    ) 
    #last date minus the evaluation period
    # minus the training time
    # plus 1 for offset/include week
    this_dt_start <- max_date - as.difftime(eval_last_weeks, unit = "weeks") - as.difftime(this_training, unit = "weeks") + as.difftime(1, unit = "weeks")
    
    #filter data
    this_ears_data <- data_events %>% 
      dplyr::filter(obs_date > this_dt_start)
    
    #create stss (must be inside settings loop due to date filtering on data for the different methods)
    this_ears_stss <- make_stss_list(geogroups, this_ears_data)
    
    #collector for geogroups loop
    this_ears_raw_results <- vector('list', length(this_ears_stss))
    
    for (i in 1:length(this_ears_stss)){
      this_ears_raw_results[[i]] <- earsC(this_ears_stss[[i]], 
                                          control = list(method = this_method, 
                                                         alpha = this_alpha, 
                                                         baseline = this_baseline))
      
    } #end geogroups loop
    
    #results for this set of settings
    #convert to dfs
    this_ears_overlaps <- calculate_overlaps(lapply(this_ears_raw_results, as.data.frame), 
                                             pre_week_search)
    
    this_ears_res <- bind_cols(events_caught(this_ears_overlaps[[1]]), 
                               events_timely(this_ears_overlaps[[1]]),
                               alarms_events(this_ears_overlaps[[2]])) %>% 
      mutate(algorithm = this_method,
             timestamp = res_dt,
             start_data = min(data_events$obs_date),
             start_eval = max_date - as.difftime(eval_last_weeks, unit = "weeks") + as.difftime(1, unit = "weeks"),
             end_eval = max_date,
             #input_settings = list(alpha = this_alpha, baseline = this_baseline),
             alpha = this_alpha,
             baseline = this_baseline,
             date_start = min(this_ears_data$obs_date),
             date_end = max(this_ears_data$obs_date)) %>% 
      #reorder
      select(algorithm, everything())
    
    #collect results
    ears_res <- dplyr::bind_rows(ears_res, this_ears_res)
    
    #save out graphs
    graph_file <- file.path(outpath, paste0("ears_graphs_", this_method, 
                                            "_base", this_baseline, 
                                            "_alpha", this_alpha, ".pdf"))
    pdf(graph_file, width = 11, height = 8.5)
    for (i in 1:length(this_ears_raw_results)){
      plot(this_ears_raw_results[[i]], 
           main = paste0("EARS Alarms: ", this_method),  
           sub = geogroups[i], #if later add names first, names(this_ears_raw)[[i]]
           legend.opts = list(horiz = TRUE, x = "topright"))
    }
    dev.off()
    
    
    #save out full data results
    #to dataframe
    this_ears_raw_out <- lapply(this_ears_raw_results, as.data.frame)
    #add names of geogroups
    this_ears_raw_out <- setNames(this_ears_raw_out, geogroups)
    
    saveRDS(this_ears_raw_out, 
            file = file.path(outpath, 
                             paste0("ears_raw_results_", this_method, 
                                    "_base", this_baseline, 
                                    "_alpha", this_alpha, ".RDS")))
    
    
  } #end ears settings loop
  
  #return
  ears_res
  
} #end run_ears_ed

# Farrington

run_farrington_ed <- function(data_events, 
                              outpath, 
                              eval_last_weeks,
                              pre_week_search,
                              settings_farrington,
                              res_dt,
                              geogroups){
  
  # data setup 
  
  #grab defaults from surveillance package 
  #  so can report out user given and default parameter values
  far_defaults <- formals(surveillance::farringtonFlexible)$control
  
  #create stss
  far_stss_list <- make_stss_list(geogroups, data_events)
  
  ##set date range for evaluation (row numbers)
  far_range <- seq(nrow(far_stss_list[[1]]) - eval_last_weeks + 1, nrow(far_stss_list[[1]]))
  
  #collector for results of each set of settings
  far_res <- tibble()
  
  #run each set of settings
  for (e in seq_along(settings_farrington)){
    
    #settings
    ## Set up this control list for Farrington (to give to surveillance function)
    # whatever is not user-given will use surveillance package defaults
    # except for b, which if not given, will do the adaptive technique as in epidemiar
    this_far_control <- list()
    
    #same range for all
    this_far_control[["range"]] <- far_range
    
    if (!is.null(settings_farrington[[e]][["w"]])){
      this_far_control[["w"]] <- settings_farrington[[e]][["w"]]
    } else {
      this_far_control[["w"]] <- far_defaults[["w"]]
    }
    if (!is.null(settings_farrington[[e]][["reweight"]])){
      this_far_control[["reweight"]] <- settings_farrington[[e]][["reweight"]]
    } else {
      this_far_control[["reweight"]] <- far_defaults[["reweight"]]
    }
    if (!is.null(settings_farrington[[e]][["weightsThreshold"]])){
      this_far_control[["weightsThreshold"]] <- settings_farrington[[e]][["weightsThreshold"]]
    } else {
      this_far_control[["weightsThreshold"]] <- far_defaults[["weightsThreshold"]]
    }
    if (!is.null(settings_farrington[[e]][["alpha"]])){
      this_far_control[["alpha"]] <- settings_farrington[[e]][["alpha"]]
    } else {
      this_far_control[["alpha"]] <- far_defaults[["alpha"]]
    }
    if (!is.null(settings_farrington[[e]][["trend"]])){
      this_far_control[["trend"]] <- settings_farrington[[e]][["trend"]]
    } else {
      this_far_control[["trend"]] <- far_defaults[["trend"]]
    }
    if (!is.null(settings_farrington[[e]][["pThresholdTrend"]])){
      this_far_control[["pThresholdTrend"]] <- settings_farrington[[e]][["pThresholdTrend"]]
    } else {
      #overriding surveillance default
      #instead of original 0.05, using 1 as in Improved
      #this_far_control[["pThresholdTrend"]] <- far_defaults[["pThresholdTrend"]]
      this_far_control[["pThresholdTrend"]] <- 1
    }
    if (!is.null(settings_farrington[[e]][["limit54"]])){
      this_far_control[["limit54"]] <- settings_farrington[[e]][["limit54"]]
    } #else {
      #this_far_control[["limit54"]] <- far_defaults[["limit54"]]
    #}
    # strange error - Error in k - control$limit54[2] : non-numeric argument to binary operator
    if (!is.null(settings_farrington[[e]][["powertrans"]])){
      this_far_control[["powertrans"]] <- settings_farrington[[e]][["powertrans"]]
    } else {
      this_far_control[["powertrans"]] <- far_defaults[["powertrans"]]
    }
    if (!is.null(settings_farrington[[e]][["fitFun"]])){
      this_far_control[["fitFun"]] <- settings_farrington[[e]][["fitFun"]]
    } else {
      this_far_control[["fitFun"]] <- far_defaults[["fitFun"]]
    }
    if (!is.null(settings_farrington[[e]][["populationOffset"]])){
      this_far_control[["populationOffset"]] <- settings_farrington[[e]][["populationOffset"]]
    } else {
      this_far_control[["populationOffset"]] <- far_defaults[["populationOffset"]]
    }
    if (!is.null(settings_farrington[[e]][["noPeriods"]])){
      this_far_control[["noPeriods"]] <- settings_farrington[[e]][["noPeriods"]]
    } else {
      this_far_control[["noPeriods"]] <- far_defaults[["noPeriods"]]
    }
    if (!is.null(settings_farrington[[e]][["pastWeeksNotIncluded"]])){
      this_far_control[["pastWeeksNotIncluded"]] <- settings_farrington[[e]][["pastWeeksNotIncluded"]]
    } else {
      this_far_control[["pastWeeksNotIncluded"]] <- far_defaults[["pastWeeksNotIncluded"]]
    }
    if (!is.null(settings_farrington[[e]][["thresholdMethod"]])){
      this_far_control[["thresholdMethod"]] <- settings_farrington[[e]][["thresholdMethod"]]
    } else {
      #overriding default ("delta") with thresholdMethod = "nbPlugin" as default instead
      this_far_control[["thresholdMethod"]] <- "nbPlugin"
    }
    
    #test if b settings given, set up for adaptive b per week
    if (!is.null(settings_farrington[[e]][["b"]])){
      adaptive_b <- FALSE
      #take user settings
      this_far_control[["b"]] <- settings_farrington[[e]][["b"]]
    } else {
      #no user value, will adjust to take maximum number of years as possible
      adaptive_b <- TRUE
    }
    
    
    #collector for geogroups loop
    this_far_results <- vector('list', length(far_stss_list))
    
    #each geogroup
    for (i in 1:length(far_stss_list)){
      
      if (adaptive_b) {
        #run per each week in evaluation period
        
        #week setup
        this_week_far_control <- this_far_control
        this_week_collector <- list()
        
        for (wk in seq_along(this_far_control[["range"]])){
          
          #eval 'range' will be just this week
          this_week_far_control[["range"]] <- this_far_control[["range"]][wk]
          
          #adaptive b
          # calculate maximum number of years previous data available
          #get date of current week
          this_week_dt <- as.data.frame(far_stss_list[[i]])[this_week_far_control[["range"]], "epoch"]
          # includes allowance for window value, w # of week
          if (!is.null(this_week_far_control[["w"]])){
            this_adj_dt <- this_week_dt - lubridate::weeks(this_week_far_control[["w"]])
          } else {
            #allow for default w = 3
            this_adj_dt <- this_week_dt - lubridate::weeks(3)
          }
          #calculate number of years difference between earliest available data date 
          # and adjusted date "start"
          #    using interval(), because this allows time_length() to deal with leap years, etc.
          this_geogroup_min_date <- far_stss_list[[i]] %>% as.data.frame() %>% pull(epoch) %>% min()
          #set adaptive b value
          this_week_far_control[["b"]] <- lubridate::interval(this_geogroup_min_date, this_adj_dt) %>%
            #year difference
            lubridate::time_length(unit = "years") %>% 
            #and floor (minimum integer year value)
            floor()
          
          #run this week
          this_week_collector[[wk]] <- farringtonFlexible(far_stss_list[[i]], 
                                                          control = this_week_far_control)
          
          
        } #end week loop
        
        #have sts per week of i geogroup, need to make 1 sts
        this_geogroup_df <- this_week_collector %>% 
          #turn into dataframes to be able to combine
          lapply(as.data.frame) %>% 
          #bind all together into one df
          dplyr::bind_rows()
        
        #turn combined df into one sts
        this_far_results[[i]] <- surveillance::sts(observed = dplyr::select(this_geogroup_df, observed) %>%
                                                     as.matrix(),
                                                   start = far_stss_list[[i]]@start,
                                                   frequency = far_stss_list[[i]]@freq,
                                                   epochAsDate = TRUE,
                                                   epoch = as.numeric(this_geogroup_df$epoch),
                                                   population = dplyr::select(this_geogroup_df, population) %>%
                                                     as.matrix(),
                                                   state = this_geogroup_df[["state"]],
                                                   alarm = this_geogroup_df[["alarm"]],
                                                   upperbound = this_geogroup_df[["upperbound"]])
        
        
      } else {
        
        #just run each geogroup with constant b
        this_far_results[[i]] <- farringtonFlexible(far_stss_list[[i]], 
                                                    control = this_far_control)
        
      } #end else of adaptive_b
    } #end geogroup loop
    
    
    #results for this set of settings
    #convert to dfs
    this_far_overlaps <- calculate_overlaps(lapply(this_far_results, as.data.frame), 
                                            pre_week_search)
    
    this_far_res <- bind_cols(events_caught(this_far_overlaps[[1]]),
                              events_timely(this_far_overlaps[[1]]),
                              alarms_events(this_far_overlaps[[2]])) %>% 
      mutate(algorithm = "Farrington",
             timestamp = res_dt,
             start_data = min(data_events$obs_date),
             start_eval = max(data_events$obs_date) - as.difftime(eval_last_weeks, 
                                                                    unit = "weeks") + as.difftime(1, unit = "weeks"),
             end_eval = max(data_events$obs_date),
             #input settings
             w = this_far_control[["w"]],
             b = this_far_control[["b"]],
             noPeriods = this_far_control[["noPeriods"]],
             trend = this_far_control[["trend"]],
             pThresholdTrend = this_far_control[["pThresholdTrend"]],
             reweight = this_far_control[["reweight"]],
             weightsThreshold = this_far_control[["weightsThreshold"]],
             pastWeeksNotIncluded = this_far_control[["pastWeeksNotIncluded"]],
             alpha_farrington = this_far_control[["alpha"]],
             limit54 = this_far_control[["limit54"]],
             powertrans = this_far_control[["powertrans"]],
             fitFun = this_far_control[["fitFun"]],
             populationOffset = this_far_control[["populationOffset"]],
             thresholdMethod = this_far_control[["thresholdMethod"]],
             b_label = if(is.null(this_far_control[["b"]])) {"ADP"
             } else {as.character(this_far_control[["b"]])}, 
             pastWeeksNotIncluded_label = if(is.null(this_far_control[["pastWeeksNotIncluded"]])) {"w"
             } else {as.character(this_far_control[["pastWeeksNotIncluded"]])},
             farID = paste0(w, b_label, noPeriods, trend, pThresholdTrend, reweight, weightsThreshold,
                            pastWeeksNotIncluded_label, populationOffset, alpha_farrington, thresholdMethod)) %>% 
      #reorder
      select(algorithm, everything())
    
    ### Figure out if can't run: all FALSE? all failures?  
    
    #collect results
    far_res <- bind_rows(far_res, this_far_res)

    #save out graphs
    graph_file <- file.path(outpath, 
                            paste0("farrington_graphs",  
                                   "_set", e, 
                                   "_w", this_far_control[["w"]], 
                                   "_b", (if(is.null(this_far_control[["b"]])) {"ADP"
                                   } else {this_far_control[["b"]]}), 
                                   "_noPeriods", this_far_control[["noPeriods"]],
                                   "_trend", this_far_control[["trend"]],
                                   "_pastexcl", (if(is.null(this_far_control[["pastWeeksNotIncluded"]])) {"w"
                                   } else {this_far_control[["pastWeeksNotIncluded"]]}),
                                   ".pdf"))
    pdf(graph_file, width = 11, height = 8.5)
    for (i in 1:length(this_far_results)){
      plot(this_far_results[[i]], 
           main = "Farrington Alarms",  
           sub = geogroups[i], 
           legend.opts = list(horiz = TRUE, x = "topright"))
    }
    dev.off()
    
    
    #save out full data results
    #to dataframe
    this_far_out <- lapply(this_far_results, as.data.frame)
    #add names of geogroups
    this_far_out <- setNames(this_far_out, geogroups)
    
    saveRDS(this_far_out, 
            file = file.path(outpath, 
                             paste0("farrington_results",  
                                    "_set", e, 
                                    "_w", this_far_control[["w"]], 
                                    "_b", (if(is.null(this_far_control[["b"]])) {"ADP"
                                    } else {this_far_control[["b"]]}), 
                                    "_noPeriods", this_far_control[["noPeriods"]],
                                    "_trend", this_far_control[["trend"]],
                                    "_pastexcl", (if(is.null(this_far_control[["pastWeeksNotIncluded"]])) {"w"
                                    } else {this_far_control[["pastWeeksNotIncluded"]]}),
                                    ".RDS")))
    
    
    
  } #end farrington settings loops
  
  #return
  far_res
  
} #end farrington function


# WHO
#https://apps.who.int/iris/bitstream/handle/10665/272284/9789241565578-eng.pdf?ua=1

run_who_ed <- function(data_events, 
                       outpath, 
                       eval_last_weeks,
                       pre_week_search,
                       settings_who,
                       res_dt,
                       geogroups){
  
  #convert into sts for similar formatting to others
  who_stss_list <- make_stss_list(geogroups, data_events, real_dates = TRUE)
  #immediately convert back to df
  who_df_list <- lapply(who_stss_list, as.data.frame)
  
  #collector for results of each set of settings
  who_res <- tibble()
  
  #run each set of settings
  for (e in seq_along(settings_who)){
    
    if (!is.null(settings_who[[e]][["method"]])){
      this_method <- settings_who[[e]][["method"]]
    } else {
      this_method <- "percentile75"
    }
    
    if (!is.null(settings_who[[e]][["n_years"]])){
      #if years given by user
      this_n_years <- settings_who[[e]][["n_years"]]
      adaptive_n_years <- FALSE
    } else {
      #per week use max years
      this_n_years <- NA
      adaptive_n_years <- TRUE
    }
    
    #collector for geogroups loop
    this_who_results <- vector('list', length(who_df_list))
    this_who_results_full <- vector('list', length(who_df_list))
    
    #each geogroup
    for (i in 1:length(who_df_list)){
      
      #this geogroup data
      this_geo_data <- who_df_list[[i]] %>% 
        dplyr::mutate(year_iso = lubridate::isoyear(epoch),
                      week_iso = lubridate::isoweek(epoch))
      
      #pull last weeks - rows that will be evaluated based on historical data
      this_geo_eval_df <- this_geo_data[(nrow(this_geo_data) - eval_last_weeks + 1):nrow(this_geo_data),]
      
      #week
      this_week_collector <- list()
      
      for (wk in 1:nrow(this_geo_eval_df)){
        
        #eval week data
        this_week_data <- this_geo_eval_df[wk,]
        this_week_date <- this_week_data$epoch
        
        #historical prior to this week
        this_historical_data <- this_geo_data %>% 
          filter(epoch < this_week_date)
        
        #add year filter (last n_years of data)
        if (!adaptive_n_years) {
          this_historical_data <- this_historical_data %>% 
            dplyr::filter(year_iso >= max(year_iso) - this_n_years)
        } #otherwise, no filter so takes all available data
        
        #do all the calculations
        this_hx_stats <- this_historical_data %>% 
          #by week of year  #<<>>handle week 53 better
          dplyr::group_by(week_iso) %>%
          dplyr::summarise(year_count = n(),
                           mean_obs = mean(observed, na.rm=TRUE),
                           sd_obs = sd(observed, na.rm=TRUE),
                           mean2sd = mean_obs + (2 * sd_obs),
                           median_obs = median(observed, na.rm=TRUE),
                           percentile75 = quantile(observed, probs = c(0.75), na.rm=TRUE),
                           percentile85 = quantile(observed, probs = c(0.85), na.rm=TRUE))
        
        #join and calculate alarm
        this_week_results <- this_week_data %>% 
          dplyr::left_join(this_hx_stats, by = c("week_iso")) %>% 
          dplyr::mutate(
            upperbound = case_when(
              this_method == "median" ~ median_obs,
              this_method == "mean2sd" ~ mean2sd,
              this_method == "percentile75" ~ percentile75,
              this_method == "percentile85" ~ percentile85),
            alarm = observed > upperbound)
        
        #append
        this_week_collector <- bind_rows(this_week_collector, this_week_results)
        
      } #end week loop
      
      #collect and turn into sts
      this_who_results[[i]] <- surveillance::sts(observed = dplyr::select(this_week_collector, observed) %>%
                                                   as.matrix(),
                                                 start = who_stss_list[[i]]@start,
                                                 frequency = who_stss_list[[i]]@freq,
                                                 epochAsDate = TRUE,
                                                 epoch = as.numeric(this_week_collector$epoch),
                                                 population = dplyr::select(this_week_collector, population) %>%
                                                   as.matrix(),
                                                 state = this_week_collector[["state"]],
                                                 alarm = this_week_collector[["alarm"]],
                                                 upperbound = this_week_collector[["upperbound"]])
      
      #also keep full data
      this_who_results_full[[i]] <- this_week_collector
      
      
    } #end geogroup loop
    
    #results for this set of settings
    #convert to dfs
    this_who_overlaps <- calculate_overlaps(lapply(this_who_results, as.data.frame), 
                                            pre_week_search)
    
    this_who_res <- bind_cols(events_caught(this_who_overlaps[[1]]),
                              events_timely(this_who_overlaps[[1]]),
                              alarms_events(this_who_overlaps[[2]])) %>% 
      mutate(algorithm = "WHO",
             timestamp = res_dt,
             start_data = min(data_events$obs_date),
             start_eval =   max(data_events$obs_date) - as.difftime(eval_last_weeks, unit = "weeks") + as.difftime(1, unit = "weeks"),
             end_eval = max(data_events$obs_date),
             #input settings
             method = this_method,
             n_years = this_n_years) %>% 
      #reorder
      select(algorithm, everything())
    
    
    #collect settings results
    who_res <- bind_rows(who_res, this_who_res)
    
    
    #save out graphs
    graph_file <- file.path(outpath, 
                            paste0("who_graphs",  
                                   "_method", this_method, 
                                   "_n_years",(if(is.na(this_n_years)) {"ADP"
                                   } else {this_n_years}), 
                                   ".pdf"))
    pdf(graph_file, width = 11, height = 8.5)
    for (i in 1:length(this_who_results)){
      plot(this_who_results[[i]], 
           main = "WHO Alarms",  
           sub = geogroups[i], 
           legend.opts = list(horiz = TRUE, x = "topright"))
    }
    dev.off()
    
    #save raw results
    #add names of geogroups
    this_who_out <- setNames(this_who_results_full, geogroups)
    
    saveRDS(this_who_out, 
            file = file.path(outpath, 
                             paste0("who_results",  
                                    "_method", this_method, 
                                    "_n_years",(if(is.na(this_n_years)) {"ADP"
                                    } else {this_n_years}), 
                                    ".RDS")))
    
  } #end who settings loop
  
  #return
  who_res
  
} #end who function


# functions for any algorithm
#STS surveillance helper functions
#NewPCODE hardcoded for now
make_stss_list <- function(geogroup_list, input_data, real_dates = TRUE){
  stss <- list()
  for (i in 1:length(geogroup_list)){
    this_geogroup <- geogroup_list[i]
    this_data <- dplyr::filter(input_data, NewPCODE == this_geogroup) %>% 
      dplyr::arrange(obs_date)
    this_cases <- as.matrix(dplyr::select(this_data, case_interp))
    this_pop <- as.matrix(dplyr::select(this_data, pop))
    
    if (real_dates){
      stss[[i]] <- surveillance::sts(observed = this_cases,
                                     start = c(isoyear(min(this_data$obs_date)), 
                                               isoweek(min(this_data$obs_date))),
                                     frequency = 52,  #weekly
                                     epochAsDate = TRUE, #have actual dates
                                     epoch = as.numeric(this_data$obs_date),
                                     population = this_pop,
                                     #Event flag in input data is in 'event' column
                                     state = this_data$event)
      
    } else {
      stss[[i]] <- surveillance::sts(observed = this_cases,
                                     start = c(isoyear(min(this_data$obs_date)), 
                                               isoweek(min(this_data$obs_date))),
                                     frequency = 52,  #weekly
                                     epochAsDate = FALSE, 
                                     #epoch = as.numeric(w_data$Date),
                                     population = this_pop,
                                     state = this_data$event)
      #thought this was going to be necessary in Farrington methods, but apparently not. Currently unused option.
      
    } #end else in the if real_dates
  } # end for loop
  stss
}

id_alarms <- function(algo_df, algo_rle){
  #identify alarms: runs of alarm weeks
  algo_df %>%
    #create ID for all runs (True or False)
    mutate(runID = with(algo_rle, { rep(seq_along(lengths), lengths)})) %>%
    #pull out TRUE runs
    filter(alarm == TRUE) %>%
    #renumber for positive runs, not just any run (T or F)
    mutate(alarmID = match(runID, unique(runID))) %>%
    #now group and summarize start and end dates
    group_by(alarmID) %>%
    summarize(start_dt = min(epoch),
              end_dt = max(epoch)) %>%
    mutate(intv = interval(start_dt, end_dt))
}

id_events <- function(algo_df, algo_rle, pre_week_search){
  #identify events: runs of event weeks
  algo_df %>%
    #create ID for all runs (True or False)
    mutate(runID = with(algo_rle, { rep(seq_along(lengths), lengths)})) %>%
    #pull out TRUE runs
    filter(state == TRUE) %>%
    #renumber for positive runs, not just any run (T or F)
    mutate(eventID = match(runID, unique(runID))) %>%
    #now group and summarize start and end dates
    group_by(eventID) %>%
    summarize(start_dt = min(epoch),
              end_dt = max(epoch)) %>%
    mutate(intv = interval(start_dt, end_dt),
           search_start_dt = start_dt - as.difftime(pre_week_search, unit = "weeks"),
           search_intv = interval(search_start_dt, end_dt))
}

calculate_event_overlaps <- function(algo_events, algo_alarms){
  #Calculate for each event, if an alarm overlaps, and get the first date the alarm run was triggered
  #using loops, because easier to walk through logic (3 loops deep: woreda, event, alarm) 
  for(w in 1:length(algo_events)){
    #check that events exist for woreda, if not, skip to next
    if (nrow(algo_events[[w]]) == 0) next
    #add info columns
    algo_events[[w]][, "num_catches"] <- 0
    algo_events[[w]][, "first_alarm"] <- as.Date(NA)
    algo_events[[w]][, "caught"] <- FALSE
    algo_events[[w]][, "timely_wks"] <- difftime(NA, NA, units = "weeks") %>% as.integer() #as.period(units = "week")
    for(e in 1:nrow(algo_events[[w]])){
      #reset 
      alarm_catch <- 0
      first_alarm <- as.Date(NA)
      
      #check that alerts exist for woreda, if so, do alert loop
      if (nrow(algo_alarms[[w]]) > 0) {
        for(a in 1:nrow(algo_alarms[[w]])){
          overlap_this_row <- int_overlaps(algo_events[[w]][[e, "search_intv"]],  algo_alarms[[w]][[a, "intv"]])
          alarm_catch <- alarm_catch + overlap_this_row
          first_this_row <- if (overlap_this_row) algo_alarms[[w]][[a, "start_dt"]] else as.Date(NA)
          first_alarm <- min(first_alarm, first_this_row, na.rm = TRUE)
        } #end alarm loop
      } #end if alerts exist
      
      #update for that event
      algo_events[[w]][e, "num_catches"] <- alarm_catch
      algo_events[[w]][e, "first_alarm"] <- first_alarm
      #calculations from info
      algo_events[[w]][e, "caught"] <- if (algo_events[[w]][[e, "num_catches"]] > 0) TRUE else FALSE
      algo_events[[w]][e, "timely_wks"] <- if (algo_events[[w]][[e, "caught"]]) {
        difftime(algo_events[[w]][[e, "first_alarm"]], algo_events[[w]][[e, "start_dt"]], units = "weeks") %>% as.integer()
      } else {
        NA_integer_
      } #end if
      
    } #end event loop
  } #end woreda loop
  algo_events
}

events_caught <- function(df_event_overlaps){
  df_event_overlaps %>% 
    #combine list of woredas into 1 df
    do.call(rbind, .) %>% 
    #summarize all events
    summarize(num_events = n(),
              num_caught = sum(caught),
              perc_caught = num_caught / num_events * 100,
              perc_missed = 100 - perc_caught) 
}  

events_timely <- function(df_event_overlaps){
  freqs <- df_event_overlaps %>%
    do.call(rbind, .) %>% 
    filter(caught == TRUE) %>% 
    select(timely_wks) %>% 
    #frequencies of each value
    table() %>%
    as.data.frame()
  
  #if no events caught, freqs is 0rows and throws errors below
  if (nrow(freqs) < 1) {
    freqs <- data.frame("weeks" = 0, "freq" = 0)
  } else {
    #table() results end up without a column name ("."). Had to go to data frame to rename via colnames, b/c rename() was not working.
    colnames(freqs) <- c("weeks", "freq")
    #turn table() generated factors into numeric values (of weeks)
    freqs$weeks <- as.numeric(as.character(freqs$weeks))
  }
  
  #manipulate categories of delays
  freqs %>% 
    summarize(preevent = sum(freq[weeks<0]),
              preevent_prc = preevent / sum(freq) * 100,
              no_delay = sum(freq[weeks==0]),
              no_delay_prc = no_delay / sum(freq) * 100,
              minor_delay = sum(freq[weeks==1]),
              minor_daly_prc = minor_delay / sum(freq) * 100,
              delayed = sum(freq[weeks>=2]),
              delayed_prc = delayed / sum(freq) * 100,
              median_delay = median(freq, na.rm = TRUE),
              mean_delay = mean(freq, na.rm = TRUE),
              not_delayed = preevent + no_delay + minor_delay,
              not_delayed_prc = not_delayed / sum(freq) * 100)
}

calculate_alarm_overlaps <- function(algo_alarms, algo_events){
  #Calculate for each alarm, if an event overlapped
  #using loops, because easier to walk through logic (3 loops deep: woreda, alarm, event) 
  for(w in 1:length(algo_alarms)){
    #check that events exist for woreda, if not, skip to next
    if (nrow(algo_alarms[[w]]) == 0) next
    
    #add info columns
    algo_alarms[[w]][, "num_events"] <- 0
    algo_alarms[[w]][, "with_event"] <- FALSE
    
    for(a in 1:nrow(algo_alarms[[w]])){
      #reset 
      event_catch <- 0
      
      #check that events exist for woreda, if so, do alert loop
      if (nrow(algo_events[[w]]) > 0) {
        #for this alarm, check each event
        for(e in 1:nrow(algo_events[[w]])){
          overlap_this_row <- int_overlaps(algo_alarms[[w]][[a, "intv"]],  algo_events[[w]][[e, "search_intv"]])
          event_catch <- event_catch + overlap_this_row
        } #end event loop
      } #end if events exist
      
      #update for that alarm
      algo_alarms[[w]][a, "num_events"] <- event_catch
      #calculations from info
      algo_alarms[[w]][a, "with_event"] <- if (algo_alarms[[w]][[a, "num_events"]] > 0) TRUE else FALSE
      
    } #end alarm loop
  } #end woreda loop
  algo_alarms
}

alarms_events <- function(df_alarm_overlaps){
  
  alarm_over_bind <- df_alarm_overlaps %>% 
    #combine list of woredas into 1 df
    do.call(rbind, .)
  
  #allowance for no alarms & events overlaps
  if (nrow(alarm_over_bind) > 0) {
    alarm_event_results <- df_alarm_overlaps %>% 
      #combine list of woredas into 1 df
      do.call(rbind, .) %>% 
      #summarize all events
      summarize(num_alarms = n(),
                num_truepos = sum(with_event),
                perc_truepos = num_truepos / num_alarms * 100,
                perc_falsepos = 100 - perc_truepos) 
  } else {
    #empty data
    alarm_event_results <- df_alarm_overlaps %>% 
      #combine list of woredas into 1 df
      do.call(rbind, .) %>% 
      #summarize all events
      summarize(num_alarms = NA_real_,
                num_truepos = NA_real_,
                perc_truepos = NA_real_,
                perc_falsepos = NA_real_)
  }
  
  alarm_event_results
  
}  

calculate_overlaps <- function(algo_dfs, pre_week_search){
  
  #RLE for alarm runs
  algo_alarm_rles <- lapply(algo_dfs, function(x) rle(x$alarm))
  algo_event_rles <- lapply(algo_dfs, function(x) rle(x$state))
  
  #use RLEs to create event and alarm dfs
  algo_alarms <- mapply(FUN = id_alarms, 
                        algo_dfs, 
                        algo_alarm_rles, 
                        SIMPLIFY = FALSE)
  algo_events <- mapply(FUN = id_events, 
                        algo_dfs, 
                        algo_event_rles, 
                        MoreArgs = list(pre_week_search), 
                        SIMPLIFY = FALSE)
  
  #calculate overlaps: Events with (or without) alarms
  algo_event_overlaps <- calculate_event_overlaps(algo_events, algo_alarms)
  #calculate overlaps: Alarms with (or without) events
  algo_alarm_overlaps <- calculate_alarm_overlaps(algo_alarms, algo_events)
  
  #return list with event and alarm overlaps
  return(list(algo_event_overlaps, algo_alarm_overlaps))
}


