########################################################################################################
#
# Code accompaniment to article in BMC Public Health:
# Comparing Malaria Early Detection Methods in a Declining Transmission Setting in Northwestern Ethiopia
#
# Code Author: Dawn Nekorchuk, dawn.nekorchuk@ou.edu
# Date: 2021-04-01
#
# This file: Functions for the Trend Weighted Seasonal Thresholds (TWST) 
#  identification of malaria events of interest.
#
########################################################################################################

# TWST algorithm description -------------------------------------------------------------
# The approach was designed to identify events retrospectively in the context of
# seasonal patterns and decreasing long-term trends in disease transmission,
# while allowing for variation in patterns across woredas as well as slight
# time-shifts in seasonal peaks between years.

# The TWST approach identified two thresholds, weekly and yearly, for each
# woreda. This combination of weekly and yearly thresholds has been used in
# other work for defining malaria epidemics [McKelvie 2012].
#McKelvie WR, Haghdoost AA, Raeisi A. Defining and detecting malaria epidemics
#in south-east Iran. Malar J. 2012;11:81.

## Preparation: 

#Raw weekly cases time series smoothed with a 5-week, centered, triangular
#smooth filter, then incidence was calculated.

## Yearly Threshold: 

#Harmonic mean of the entire year plus a multiplier based on the standard
#deviation (`yr_sd_multiplier`)

## Weekly Threshold (3 steps): 

#1. Raw threshold value was the harmonic mean of that week in the year, plus
#multiplier based on the standard deviation (`wk_sd_multiplier`).

#2. Raw thresholds were conditionally weighted based on the trend of the year
#harmonic mean. If there was a declining trend (from the year previous), then
#the weekly threshold values were weighted proportional to the difference
#between the current year harmonic mean and the highest (max) harmonic mean
#using a weighting factor (`thr_weighting`): (max – thr_weighting * (max –
#current)) / max. Note: weight simply based on current / max is too drastic. If
#there was no declining trend, the weekly thresholds were weighted based on the
#previous year mean instead of the current year. I.e., if current year mean is
#greater, use previous year mean for weighting. Note: this does not propagate
#more than 1 year forward. Note: 1st year will always use 1st yr mean / highest
#yr mean as weight. This allows for changing expectations of general decreasing
#(overall) trend in malaria, but also sensitive to localized increases compared
#to previous years

#3. Allowances are made to prevent minor time shifts in increasing and
#decreasing case counts between years from triggering alerts by inflating weekly
#thresholds that were not near in time to peaks. Peak areas were identified
#using a percentile cut-off per year (`yr_peak_percent`), plus short stretches
#(up to 8 weeks) between these high rates. The inflation was based on the
#average of the year and week harmonic means multiplied by an expansion factor
#(`nonpeak_expansion`), which was then added to the trend weighted threshold of
#the previous step to arrive at the final TWST week threshold.

## Event identification

# Anomalies were identified if cases exceeded both the yearly and weekly
# thresholds, and events were identified if anomalies lasted for four or more
# weeks consecutively. Events that were separated by only one or two weeks were
# merged into one event.


# FUNCTIONS ----------------------------------------------

#assumed fields:
#NewPCODE: ID for geo groups
#pop: population count for that observation week - NewPCODE
#R_NAME, Z_NAME, W_NAME: descriptive names (NewPCODE)
#obs_date: date of observation
#week: ISO week of obs_date
#year: ISO year of obs_Date


#master TWST run function
run_twst <- function(df, 
                     field_cases, 
                     out_dir,
                     thr_weighting, 
                     wk_sd_multiplier, 
                     yr_sd_multiplier, 
                     yr_peak_percent, 
                     nonpeak_expansion){
  
  #quosures
  quo_field_cases <- enquo(field_cases)

  #result date/time
  res_dt <- format(Sys.time(), "%Y%m%d%H%M")
  
  #prep data for twst 
  to_twst <- df %>% 
    #interpolate NA values
    mutate(case_interp = !!quo_field_cases) %>% 
    #5 week rolling mean to smooth cases
    rollmean() %>% 
    #calculate original/observed incidence
    mutate(obs_inc = !!quo_field_cases / pop * 1000)
  
  #create week summaries
  twst_wk <- create_twst_week(to_twst, 
                              wk_sd_multiplier)
  #create year summaries
  twst_yr <- create_twst_year(to_twst, 
                              thr_weighting,
                              yr_sd_multiplier) 

  #join all together
  twst_final <- join_all_twst(to_twst, 
                              twst_wk, 
                              twst_yr) %>% 
    #expand non-peaks
    expand_nonpeak(yr_peak_percent, 
                   nonpeak_expansion) %>% 
    #calc events
    calc_alerts_events()
  
  #save out
  save_twst(twst_final, 
            res_dt, 
            quo_field_cases,
            out_dir, 
            thr_weighting, 
            wk_sd_multiplier, 
            yr_sd_multiplier, 
            yr_peak_percent, 
            nonpeak_expansion)
  
  #create and save graphs
  plot_twst(twst_final, 
            res_dt, 
            quo_field_cases,
            out_dir, 
            thr_weighting, 
            wk_sd_multiplier, 
            yr_sd_multiplier, 
            yr_peak_percent, 
            nonpeak_expansion)
  
  #rles
  run_twst_rles(twst_final, 
                res_dt, 
                quo_field_cases, 
                out_dir)
  
  #dot plots
  #Note: dot plots use R_NAME, W_NAME, Z_NAME instead of NewPCODE/field_id
  make_twst_stackedevents(twst_final, 
                          res_dt, 
                          quo_field_cases, 
                          out_dir)
  
  return(twst_final)
}


# Rolling mean, weighted, centered, 5 time points
roll_mean5_weight <- function(x) {
  weights <- c(1/9, 2/9, 3/9, 2/9, 1/9)
  as.numeric(stats::filter(x, weights, sides=2))
}

# Geometric mean, vectorized, zero and NA tolerant. With optional zero propagation.
gm_mean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  # https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
}

# Function to look for runs of a minimum length 
is_run <- function(x, look_for = 1, min_length = 4) {
  rle.x <- rle(x)
  meet_criteria <-
    rle.x$lengths >= min_length &
    rle.x$values == look_for
  rle.x$values <- rep(FALSE, length(rle.x$values))
  rle.x$values[meet_criteria] <- TRUE
  inverse.rle(rle.x)
}

#capture just 1 legend as grob
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
  # https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
}

#TWST internal functions
# #Fill NAs values by interpolation
# fill_NA <- function(df, 
#                     quo_field_cases){
#   df %>% 
#     group_by(NewPCODE) %>% 
#     arrange(NewPCODE, obs_date) %>% 
#     dplyr::mutate(case_interp = zoo::na.approx(!!quo_field_cases, rule=1:1, na.rm = FALSE))
# }

#Smoothing: moving average / rolling mean
rollmean <- function(df){
  df %>% 
    #group
    group_by(NewPCODE) %>% 
    #make sure ordering is correct before moving average
    arrange(NewPCODE, obs_date) %>% 
    #moving average/rolling mean
    mutate(smooth = roll_mean5_weight(case_interp),
           smooth_inc = smooth / pop * 1000) %>% 
    #trim off ends (where rolling average cannot be calculated properly)
    filter(obs_date > min(obs_date)+2*7 & obs_date < max(obs_date)-2*7)
}

#calculate week summary values for TWST
create_twst_week <- function(twst_df, 
                             wk_sd_multiplier){
  #calc smoothed weekly log harmonic means and (arithmetic) standard deviations
  twst_df %>% 
    group_by(NewPCODE, week) %>%
    # harmonic mean zero = FALSE to calc mean on non-zero items (and not just return zero if any zeros present)
    summarize(wk_inc_hmean = psych::harmonic.mean(smooth_inc, na.rm = TRUE, zero = FALSE),
              wk_inc_sd = sd(smooth_inc, na.rm = TRUE)) %>% 
    #b/c iso wk 53 only happens rarely, fill in with wk 52 (via making NA, and fill down)
    mutate(wk_inc_hmean = ifelse(week == 53, NA, wk_inc_hmean),
           wk_inc_sd = ifelse(week == 53, NA, wk_inc_sd)) %>%
    tidyr::fill(wk_inc_hmean, wk_inc_sd, .direction = "down") %>% 
    # set raw thresholds (hmean + k SD)
    mutate(wk_inc_raw_thr = wk_inc_hmean + wk_sd_multiplier * wk_inc_sd)
}

#calculate year summary values for TWST
create_twst_year <- function(twst_df, 
                             thr_weighting,
                             yr_sd_multiplier){
  #calc harmonic annual mean, yearly thresholds
  #and appropriate year operational hmeans for weighting
  
  twst_df %>%
    #confirming ordering for later lag call
    arrange(NewPCODE, year) %>% 
    #group to yr to calc
    group_by(NewPCODE, year) %>% 
    # harmonic mean zero = FALSE to calc mean on non-zero items (and not just return zero if any zeros present)
    summarize(yr_inc_hmean = psych::harmonic.mean(smooth_inc, na.rm = TRUE, zero = FALSE),
              yr_inc_sd = sd(smooth_inc, na.rm = TRUE),
              #handle if data starts (or ends) mid-year
              yr_completeness = n() / 52) %>%
    #adjust for partial year data
    mutate(yr_inc_hmean = yr_inc_hmean * yr_completeness, 
           yr_inc_sd = yr_inc_sd * yr_completeness) %>% 
    # #<<>> 2012 OVERRIDE, since half year.
    # # TODO: build in more appropriate ways of handling partial year data
    # mutate(yr_inc_hmean = if(year == 2012) {yr_inc_hmean / 2} else {yr_inc_hmean},
    #        yr_inc_sd = if(year == 2012) {yr_inc_sd / 2} else {yr_inc_sd}) %>% 
    #ungroup and regroup to just geo id
    ungroup() %>% 
    group_by(NewPCODE) %>% 
    mutate(prev_inc_hmean = lag(yr_inc_hmean),   # lag 1 of the hmeans
           # calculate appropriate operational hmeans (current year if < previous, else previous)
           opt_inc_hmean = case_when(
             yr_inc_hmean <= prev_inc_hmean ~ yr_inc_hmean,
             yr_inc_hmean > prev_inc_hmean ~ prev_inc_hmean,
             TRUE ~ yr_inc_hmean),
           #calc max of yr_inc_means (used later for weight ratio)
           max_inc_hmean = max(yr_inc_hmean),
           #calc weight (proportional to ratio of opt to max)
           opt_inc_weight = (max_inc_hmean - (max_inc_hmean - opt_inc_hmean) * thr_weighting) / max_inc_hmean,
           #yearly thresholds
           yr_inc_thr = yr_inc_hmean + yr_sd_multiplier * yr_inc_sd)
}

#join original (weekly), week summary, and year summary tbls
join_all_twst <- function(twst_df, 
                          twst_wk_df, 
                          twst_yr_df){
  twst_df %>%
    #week info
    left_join(twst_wk_df, 
              #NSE
              by = c("NewPCODE", "week")) %>% 
    #year info
    left_join(twst_yr_df,
              by = c("NewPCODE", "year")) %>% 
    #group
    group_by(NewPCODE) %>% 
    # weight weekly thresholds by appropriate yearly value (trend weighting seasonal)
    #incidence
    mutate(wk_inc_thr_opt = wk_inc_raw_thr * opt_inc_weight)
}

#non-peak expansion
expand_nonpeak <- function(twst_all_df, 
                           yr_peak_percent, 
                           nonpeak_expansion){
  twst_all_df %>% 
    #Peak-related: if series is above certain percentile (calc per year)
    #group for calc
    group_by(NewPCODE, year) %>% 
    #percentile for part of peak identification
    mutate(yr_inc_pk = quantile(wk_inc_thr_opt, probs = yr_peak_percent, names = FALSE, na.rm = TRUE)) %>% 
    #ungroup and regroup back to just geo id
    ungroup() %>% 
    group_by(NewPCODE) %>% 
    #peak
    mutate(peak = ifelse(wk_inc_thr_opt >= yr_inc_pk, TRUE, FALSE),
           ##also peak if surrounded by high values on BOTH sides (w/in 8 weeks)
           #leading edges
           near_peak_lead = ifelse(lead(peak, default = FALSE) | lead(peak, n = 2, default = FALSE) |
                                     lead(peak, n = 3, default = FALSE) | lead(peak, n = 4, default = FALSE) |
                                     lead(peak, n = 5, default = FALSE) | lead(peak, n = 6, default = FALSE) |
                                     lead(peak, n = 7, default = FALSE) | lead(peak, n = 8, default = FALSE),
                                   TRUE, FALSE),
           #lagging edges
           near_peak_lag = ifelse(lag(peak, default = FALSE) | lag(peak, n = 2, default = FALSE) |
                                    lag(peak, n = 3, default = FALSE) | lag(peak, n = 4, default = FALSE) |
                                    lag(peak, n = 5, default = FALSE) | lag(peak, n = 6, default = FALSE) |
                                    lag(peak, n = 7, default = FALSE) | lag(peak, n = 8, default = FALSE),
                                  TRUE, FALSE),
           #if not already a peak, and BOTH lead and lag, then make part of 'peak'; otherwise use original peak value
           peak = ifelse(near_peak_lead & near_peak_lag & !(peak), TRUE, peak)) %>% 
    
    ##shoulder expansion (if not peak) 
    mutate(wk_inc_thr = ifelse(peak, 
                               wk_inc_thr_opt, 
                               ((yr_inc_hmean + wk_inc_hmean) / 2) * nonpeak_expansion + wk_inc_thr_opt))
}

#calculate alerts and events
calc_alerts_events <- function(twst_all_df_exp){
  twst_all_df_exp %>% 
    #confirm grouping
    group_by(NewPCODE) %>% 
    #add new fields for if smoothed week value exceeded the week, yearly, both thresholds
    mutate(wk_alert = ifelse(smooth_inc > wk_inc_thr, TRUE, FALSE),
           yr_alert = ifelse(smooth_inc > yr_inc_thr, TRUE, FALSE),
           alert = ifelse((wk_alert & yr_alert), TRUE, FALSE)) %>% 
    #now flag those with runs of 4 or more
    arrange(NewPCODE, obs_date) %>% 
    mutate(event_run = is_run(alert, TRUE, 4),
           #merge events if separated by 2 weeks or less
           split_event_lead = ifelse(lead(event_run, n = 1, default = FALSE) & 
                                       (lag(event_run, n = 1, default = FALSE) | lag(event_run, n = 2, default = FALSE)),
                                     TRUE, FALSE),
           split_event_lag = ifelse(lag(event_run, n = 1, default = FALSE) & 
                                      (lead(event_run, n = 1, default = FALSE) | lead(event_run, n = 2, default = FALSE)),
                                    TRUE, FALSE),
           split_event = ifelse(split_event_lead | split_event_lag, TRUE, FALSE),
           #event if alert run, or split event
           event = ifelse(event_run | split_event, TRUE, FALSE)) %>% 
    ungroup()
}

#saving function
save_twst <- function(twst_final, 
                      res_dt, 
                      quo_field_cases,
                      out_dir, 
                      thr_weighting, 
                      wk_sd_multiplier, 
                      yr_sd_multiplier, 
                      yr_peak_percent, 
                      nonpeak_expansion){
  #labels 
  fcase <- as_label(quo_field_cases)

  write_csv(twst_final, 
            file.path(out_dir, 
                      paste0(fcase, "_", res_dt,".csv")))
  
  saveRDS(twst_final, 
          file.path(out_dir, 
                    paste0(fcase, "_", res_dt,".RDS")))
  
  param_df <- data.frame(matrix(nrow = 7, ncol = 2, 
                                dimnames = list(c(), c("Parameter", "Value"))))
  param_df[1,] <- cbind("res_dt", paste0("YmdHM_", res_dt))
  param_df[2,] <- cbind("field_cases", fcase)
  param_df[3,] <- cbind("thr_weighting", thr_weighting)
  param_df[4,] <- cbind("wk_sd_multiplier", wk_sd_multiplier)
  param_df[5,] <- cbind("yr_sd_multiplier", yr_sd_multiplier)
  param_df[6,] <- cbind("yr_peak_percent", yr_peak_percent)
  param_df[7,] <- cbind("nonpeak_expansion", nonpeak_expansion)
  write_csv(param_df, file.path(out_dir, paste0(fcase, "_", res_dt, "_log.csv")))
}

#save plots of event graphs per woreda
plot_twst <- function(twst_final, 
                      res_dt, 
                      quo_field_cases,
                      out_dir, 
                      thr_weighting, 
                      wk_sd_multiplier, 
                      yr_sd_multiplier, 
                      yr_peak_percent, 
                      nonpeak_expansion){
  
  cf <- as_label(quo_field_cases)

  #all codes
  all_ids <- unique(twst_final$NewPCODE) %>% sort()

  twst_plots <- list()
  for (i in seq_along(all_ids)) {
    this_id <- all_ids[i]
    this_data <- filter(twst_final, NewPCODE == this_id)
    #max per geo id of either threshold or observed
    this_data$max_yval <- max(max(select(this_data, obs_inc), na.rm = TRUE),
                           max(select(this_data, wk_inc_thr), na.rm = TRUE))
    this_event <- filter(this_data, event == TRUE)
    this_alert <- filter(this_data, alert == TRUE)
    
    #labels (all same in this_data, just taking first for 1 value)
    this_region <- first(this_data$R_NAME)
    this_zone <- first(this_data$Z_NAME)
    this_woreda <- first(this_data$W_NAME)

    #plot
    twst_plots[[i]] <- ggplot(data=this_data) +
      theme_bw() +
      labs(title = this_woreda, subtitle = paste0("Region: ", this_region, ", Zone: ", this_zone)) +
      ylab("Incidence per 1000") +
      theme(axis.title.x = element_blank()) +
      #ensure that each year gets a tick and a label (and doesn't push the plot too far away from the edges of the axis)
      scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.02, 0)) +
      #not all woredas have outbreaks/alerts. with empty datasets, errors with Aesthetics lengths not matching
      #need ugly if statements (completely enclosed in {} to deal with +), and fake points for legends to work right
      {if (nrow(this_alert) > 0) geom_point(data = this_alert, aes(obs_date, 0.025*max_yval+max_yval, color = "Alert Marker"), size = 1)} +
      {if (nrow(this_alert) == 0) geom_point(aes(min(obs_date), 0, color = "Alert Marker"), alpha = 0)} +
      {if (nrow(this_event) > 0) geom_point(data = this_event, aes(obs_date, 0.06*max_yval+max_yval, color = "Event Marker"), size = 1)} +
      {if (nrow(this_event) == 0) geom_point(aes(min(obs_date), 0, color = "Event Marker"), alpha = 0)} +
      #note use of aes_ to work with quosures, add tilde to others
      #geom_line(aes_(~obs_date, quo_incidence, color = "Observed Incidence"), linetype = 1) +
      geom_line(aes(obs_date, obs_inc, color = "Observed Incidence"), linetype = 1) +
      geom_line(aes(obs_date, smooth_inc, color = "Smoothed Incidence"), linetype = 1) +
      geom_line(aes(obs_date, wk_inc_thr, color = "Expanded Week Threshold"), linetype = 2, size = 0.4) +
      geom_line(aes(obs_date, yr_inc_thr, color = "Year Threshold"), linetype = 2, size = 0.4) +
      scale_color_manual("", values = c("Event Marker" = "red",
                                        "Alert Marker" = "goldenrod2",
                                        "Observed Incidence" = "gray60",
                                        "Smoothed Incidence" = "black",
                                        "Expanded Week Threshold" = "seagreen",
                                        "Year Threshold" = "royalblue"),
                         breaks = c("Event Marker",
                                    "Alert Marker",
                                    "Observed Incidence",
                                    "Smoothed Incidence",
                                    "Expanded Week Threshold",
                                    "Year Threshold"),
                         guide = guide_legend(override.aes = list(
                           linetype = c(0, 0, 1, 1, 2, 2),
                           shape = c(16, 16, rep(NA, 4))))) 
  }
  #capture just 1 legend as grob
  common_legend <- g_legend(twst_plots[[1]])
  #remove legends from all
  for (i in 1:length(twst_plots)) {
    twst_plots[[i]] <- twst_plots[[i]] +
      theme(legend.position = 'none')
  }
  #put common legend at beginning
  twst_plots[2:(length(twst_plots)+1)] <- twst_plots[1:(length(twst_plots))]
  twst_plots[[1]] <- common_legend
  #add details in textgrob
  details <- textGrob(label = paste0("Week Threshold: ", wk_sd_multiplier, " SD. ", 
                                     "Year Threshold: ", yr_sd_multiplier, " SD. ", 
                                     "Threshold Year Weighting at ", thr_weighting, 
                                     " of difference between max and raw weight. ",
                                     "Expanded nonpeak factor ", nonpeak_expansion, 
                                     " below peak percentile ", 
                                     yr_peak_percent, "."),
                      gp = gpar(fontsize = 8))
  #arrange all the plots
  multi_11 <- marrangeGrob(twst_plots, nrow = 1, ncol = 1, top = NULL, bottom = details)
  filenm <- file.path(out_dir, paste0(cf, "_", res_dt, "_plots.pdf"))
  ggsave(filenm, multi_11, width=11, height=8.5)
}

#rles for event summaries
run_twst_rles <- function(twst_final, 
                          res_dt, 
                          quo_field_cases,
                          out_dir){
  
  cf <- as_label(quo_field_cases)

  #function so that rle to honor group_by                  
  get_rle_info <- function(x){
    x_rle <- rle(x)
    
    runID <- rep(seq_along(x_rle$lengths), times = x_rle$lengths)
    run_length <- rep(x_rle$lengths, times = x_rle$lengths)
    cond <- rep(x_rle$values, times = x_rle$lengths)
    
    bind_cols("runID" = runID, "run_length" = run_length, "cond" = cond)
  } #end get_rle_info()
  
  # per event
  twst_events <- twst_final %>% 
    #get rle results, creating ID for all runs (TRUE or FALSE runs)
    dplyr::do(cbind(., get_rle_info(.$event))) %>% 
    #pull out TRUE runs (events)
    filter(cond == TRUE) %>% 
    #renumber events (to events, not just runs)
    mutate(eventID = match(runID, unique(runID))) %>% 
    #select relevant fields
    select(NewPCODE, R_NAME, Z_NAME, W_NAME,
           obs_date, year, eventID, run_length, 
           !!quo_field_cases, obs_inc, smooth_inc, alert, event) 
  
  #create summary of event
  twst_events_stats <- twst_events %>% 
    group_by(NewPCODE, R_NAME, Z_NAME, W_NAME, eventID) %>% 
    summarize(start_date = min(obs_date),
              end_date = max(obs_date),
              num_weeks = unique(run_length),
              total_cases = sum(!!quo_field_cases, na.rm = TRUE),
              peak_cases = max(!!quo_field_cases, na.rm = TRUE),
              peak_inc = max(obs_inc, na.rm = TRUE),
              event_yr = min(year)) %>% 
    arrange(NewPCODE)
  write_csv(twst_events_stats, 
            file.path(out_dir, paste0(cf, "_", res_dt, "_per_event_stats.csv")))
  
  #events per year
  twst_events_yr_sum <- twst_events_stats %>% 
    ungroup() %>% 
    group_by(event_yr) %>% 
    summarize(nevents = n(),
              ave_wks = mean(num_weeks, na.rm = TRUE),
              median_wks = median(num_weeks, na.rm = TRUE),
              ave_cases = mean(total_cases, na.rm = TRUE),
              median_cases = median(total_cases, na.rm = TRUE))
  write_csv(twst_events_yr_sum, 
            file.path(out_dir, paste0(cf, "_", res_dt, "_per_year_event_stats.csv")))
}

#create stacked dot plot showing all events in all woredas over time
make_twst_stackedevents <- function(twst_final, 
                                    res_dt, 
                                    quo_field_cases,
                                    out_dir){
  cf <- as_label(quo_field_cases)
  
  #sorting/organizing first
  dot_data <- twst_final %>%
    mutate(rzw_fx = paste0(R_NAME, "-", Z_NAME, "-", W_NAME),
           zw_fx = paste0(W_NAME, " (", Z_NAME, ")")) %>%
    arrange(rzw_fx)
  #set the order we just did (by creating factor)
  dot_data$rzw_fx <- factor(dot_data$rzw_fx, levels=unique(as.character(dot_data$rzw_fx)))
  dot_data$zw_fx <- factor(dot_data$zw_fx, levels=unique(as.character(dot_data$zw_fx)))
  

  dot_data_ev <- filter(dot_data, event == TRUE)
  
  #need date lines on January and July during the time series
  #beginning of needed sequence
  min_month <- min(dot_data$obs_date, na.rm = TRUE) %>% month()
  min_year <- min(dot_data$obs_date, na.rm = TRUE) %>% year()
  #get January or July months, whichever comes next, with appropriate year
  st_dt_line <- if_else(min_month < 7, 
                        lubridate::make_date(year = min_year,
                                             month = 7,
                                             day = 1L),
                        lubridate::make_date(year = min_year + 1,
                                             month = 1,
                                             day = 1L))
  #end of needed sequence
  max_month <- max(dot_data$obs_date, na.rm = TRUE) %>% month()
  max_year <- max(dot_data$obs_date, na.rm = TRUE) %>% year()
  #get previous January or July, whichever came just before
  end_dt_line <- if_else(max_month < 7,
                         lubridate::make_date(year = max_year,
                                              month = 1,
                                              day = 1L),
                         lubridate::make_date(year = max_year,
                                              month = 7,
                                              day = 1L))
  #sequence of Jan/Jul dates in data range
  date_lines <- seq(st_dt_line, end_dt_line, by = "6 months")
  
  d_plot <- ggplot() +
    #background dashed line all the way across
    geom_segment(data = dot_data, 
                 aes(x = zw_fx, xend = zw_fx, y = min(obs_date), yend = max(obs_date)),
                 linetype = "dashed", size = 0.1, color = "gray60") +
    #event markers
    geom_point(data = dot_data_ev,
               aes(x = zw_fx, y = obs_date), color = "red3", size = 1.75) +
    #Date line
    geom_hline(yintercept = date_lines, color = "grey50", size = 0.2) +
    #labs(title = (expression(paste(italic("P. falciparum"), " Events")))) +
    xlab("Woreda (Zone)") +
    #ylab("Date") +
    scale_x_discrete(limits = rev(levels(dot_data$zw_fx))) +
    #ensure that each year gets a tick mark, and removes excess space
    scale_y_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.02, 0)) +
    theme_classic() +
    #axis text
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          axis.text.y = element_text(size = 6)) +
    coord_flip()
  #d_plot ##reaches elapsed time limit. -- won't display in RStudio, look at saved plot (below)
  
  ggsave(d_plot,
         filename = file.path(out_dir, paste0(cf, "_", res_dt, "_events_stacked.jpg")),
         width=8.5, height=11)
}




