############################################################
#
# Demo script for running TWST
# 
# On SYNTHETIC data for DEMONSTRATION ONLY
#
############################################################


outpath <- file.path("output", "twst_demo")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# LIBRARIES ----------------------------------------------

library(tidyverse)
library(rlang) #quosure tidy verse
library(lubridate) #date manip
library(zoo)
library(psych) #other means functions (harmonic)
library(grid) #plotting multiple
library(gridExtra) #plotting multiple
library(gtable) #manipulating grobs

#TWST functions
source("R/functions_twst_demo.R")

# DATA ---------------------------------------------------

#demo data
demo_data <- readRDS("data_demo/synthetic_data.RDS") %>% 
  #rename
  rename(pop = pop_at_risk) %>% 
  #add ISO info
  mutate(week = lubridate::isoweek(obs_date),
         year = lubridate::isoyear(obs_date))


# Trend Weighted Seasonal Thresholds (TWST) -------------------------------

# Will save out:
#  TWST results as csv and RDS
#  csvs of various statistics generated
#  jpg of stacked dot plot of events
#  pdf of each group (by id) times series with thresholds and events


#DEMO P. falciparum (not true case data)
twst_pfm <- run_twst(df = demo_data, 
                     field_cases = test_pf_tot, 
                     out_dir = outpath,
                     thr_weighting = 0.5, 
                     wk_sd_multiplier = 1, 
                     yr_sd_multiplier = 1.5, 
                     yr_peak_percent = 0.85, 
                     nonpeak_expansion = 1.2)


#DEMO P.vivax (not true case data)
twst_pv <- run_twst(df = demo_data, 
                    field_cases = test_pv_only, 
                    out_dir = outpath,
                    thr_weighting = 0.25, 
                    wk_sd_multiplier = 0.75, 
                    yr_sd_multiplier = 1, 
                    yr_peak_percent = 0.9, 
                    nonpeak_expansion = 1.1)

