
############################################################
#
# Demo script for running Event Detection Comparison
# 
# On SYNTHETIC data for DEMONSTRATION ONLY
#
############################################################


outpath <- file.path("output", "ed_demo")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)


# Libraries ---------------------------------------------------------------

library(tidyverse) #readr, dplyr, tidyr, ggplot2
library(rlang)
library(lubridate) #date manip
library(zoo)
library(surveillance) #time series and outbreak functions


#due to experimental dplyr::summarise() parameter
options(dplyr.summarise.inform=F)

#Event Detection custom functions
source("R/functions_event_detection_compare_demo.R")


# Data & Preprocessing ----------------------------------------------------

#read in TWST results
# usually from output folder, e.g.
# twst_run <- readRDS("output/twst_demo/test_pf_tot_202104021300.RDS")
# however, for data included in demo
twst_run <- readRDS("data_demo/twst_demo/test_pf_tot_202104021300.RDS")

#pull list of NewPCODEs  
pcodes <- unique(twst_run$NewPCODE) %>% sort()

#Need to pre-pend with case counts (otherwise lose time because of spin up)
#surveillance data
demo_data <- readRDS("data_demo/synthetic_data.RDS") %>% 
  #rename
  rename(pop = pop_at_risk)

#prepend ALL (needed for long term calculations)
data_pre <- demo_data %>% 
  dplyr::filter(obs_date < min(twst_run$obs_date)) %>% 
  #need NA-interpolated case field, event (fill FALSE to not throw off metrics)
  dplyr::mutate(case_interp = zoo::na.approx(test_pf_tot, na.rm = FALSE),
                event = FALSE)

#full data set to give to early detection algorithm comparison functions
twst_extd <- dplyr::bind_rows(data_pre, twst_run) %>% 
  dplyr::arrange(NewPCODE, obs_date)


# Settings sets -------------------


settings_rand_demo <- list(list(alarm_p = 0.2, buffer_wks = 4, max_alarm_wks = 5),
                        list(alarm_p = 0.1, buffer_wks = 4, max_alarm_wks = 5),
                        list(alarm_p = 0.05, buffer_wks = 4, max_alarm_wks = 5),
                        list(alarm_p = 0.025, buffer_wks = 4, max_alarm_wks = 5),
                        list(alarm_p = 0.012, buffer_wks = 4, max_alarm_wks = 5),
                        list(alarm_p = 0.006, buffer_wks = 4, max_alarm_wks = 5))

settings_ears_demo <- list(list(method = "C1", alpha = 0.001, baseline = 7),
                        list(method = "C2", alpha = 0.001, baseline = 7),
                        list(method = "C3", alpha = 0.025, baseline = 7),
                        list(method = "C1", alpha = 0.001, baseline = 14),
                        list(method = "C2", alpha = 0.001, baseline = 14),
                        list(method = "C3", alpha = 0.025, baseline = 14),
                        list(method = "C1", alpha = 0.001, baseline = 28),
                        list(method = "C2", alpha = 0.001, baseline = 28),
                        list(method = "C3", alpha = 0.025, baseline = 28),
                        list(method = "C1", alpha = 0.001, baseline = 56),
                        list(method = "C2", alpha = 0.001, baseline = 56),
                        list(method = "C3", alpha = 0.025, baseline = 56))


settings_far_demo <- list(
  #original defaults, no seasonality, no pop
  list(w = 3, b = 5, noPeriods = 1, trend = TRUE, weightsThreshold = 1, pThresholdTrend = 0.05, thresholdMethod = "delta", populationOffset = FALSE),
  #original defaults, w/ noPeriods = 4, no pop
  list(w = 3, b = 5, noPeriods = 4, trend = TRUE, weightsThreshold = 1, pThresholdTrend = 0.05, thresholdMethod = "delta", populationOffset = FALSE),
  
  #improved defaults, no seasonality, no pop
  list(w = 3, b = 5, noPeriods = 1, trend = TRUE, weightsThreshold = 2.58, pThresholdTrend = 1, pastWeeksNotIncluded = 26, populationOffset = FALSE, thresholdMethod = "nbPlugin"),
  #improved defaults, w/ noPeriods = 4, no pop
  list(w = 3, b = 5, noPeriods = 4, trend = TRUE, weightsThreshold = 2.58, pThresholdTrend = 1, pastWeeksNotIncluded = 26, populationOffset = FALSE, thresholdMethod = "nbPlugin"),
  
  #all others (npPlugin default in my function if not given)
  list(w = 3, b = 5, noPeriods = 8, trend = TRUE, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 8, trend = TRUE, pastWeeksNotIncluded = 26, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 8, trend = FALSE, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 8, trend = FALSE, pastWeeksNotIncluded = 26, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 4, trend = TRUE, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 4, trend = TRUE, pastWeeksNotIncluded = 26, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 4, trend = FALSE, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 4, trend = FALSE, pastWeeksNotIncluded = 26, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 1, trend = TRUE, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 1, trend = TRUE, pastWeeksNotIncluded = 26, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 1, trend = FALSE, populationOffset = TRUE),
  list(w = 3, b = 5, noPeriods = 1, trend = FALSE, pastWeeksNotIncluded = 26, populationOffset = TRUE),
  
  list(w = 3, b = 3, noPeriods = 8, trend = TRUE, populationOffset = TRUE),
  #etc.

  list(w = 3,        noPeriods = 8, trend = TRUE, populationOffset = TRUE),
  #etc.
  
  list(w = 5, b = 5, noPeriods = 8, trend = TRUE, populationOffset = TRUE),
  #etc.
  
  list(w = 5, b = 3, noPeriods = 8, trend = TRUE, populationOffset = TRUE),
  #etc.
  
  list(w = 5,        noPeriods = 8, trend = TRUE, populationOffset = TRUE)
  #etc.
)


settings_who_demo <- list(list(method = "median", n_years = 5),
                       list(method = "median"),
                       list(method = "mean2sd", n_years = 5),
                       list(method = "mean2sd"),
                       list(method = "percentile75", n_years = 5),
                       list(method = "percentile75"),
                       list(method = "percentile85", n_years = 5),
                       list(method = "percentile85"))



# Run ------------------------------------------------------------------------

# *** NOTE: Will take several minutes (~20-30 min) to run ALL sets *** 

# Will save RDS, PDF of time series, for each result
# Will save csv with statistics for all (detection_comparisons.csv)
# Will save an RDS with metadata (event_detection_run_metadata.RDS)

set.seed(42) #for reproducibility of random alarms

#2018 - last 52 weeks of dataset
#REM: demo dataset runs mid-2012 through 2018
demo_ed <- run_ed_comparison(data_events = twst_extd,
                             outpath_base = outpath,
                             #number of weeks at end of data to evaluate
                             eval_last_weeks = 52,
                             #allows alarms to be counted if triggered prior to event (up to 3 weeks)
                             pre_week_search = 3,
                             run_random = TRUE,
                             settings_random = settings_rand_demo,
                             run_ears = TRUE,
                             settings_ears = settings_ears_demo,
                             run_farrington = TRUE,
                             settings_farrington = settings_far_demo, 
                             run_who = TRUE,
                             settings_who = settings_who_demo)
demo_ed


