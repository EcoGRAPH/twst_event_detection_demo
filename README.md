
Project contains a code accompaniment to the article in BMC Public Health: "Comparing Malaria Early Detection Methods in a Declining Transmission Setting in Northwestern Ethiopia" by Nekorchuk et al. (2021). 

The project contains functions to perform the following two analyses which were presented in the paper:

1) Our novel Trend Weighted Seasonal Thresholds (TWST) approach which was designed to identify malaria events retrospectively in the context of seasonal patterns and decreasing long-term trends in disease transmission, while allowing for variation in patterns across geographical districts as well as slight time-shifts in seasonal peaks between years.

2) Event Detection Comparison: Comparing various Early Detection algorithms used in the paper: Random alarms (naive model), weekly statistics-based thresholds (e.g. WHO, Cullen), CDC EARS, and Farrington algorithms. 

This project also contains demo run scripts for TWST, 'run_twst_demo.R', and event detection comparison, 'run_ed_compare_demo.R', with synthetic data for demonstration purposes ONLY. 

---
This work is part of a larger project, Epidemic Prognosis Incorporating Disease and Environmental Monitoring for Integrated Assessment (EPIDEMIA). The EPIDEMIA Forecasting System integrate surveillance and environmental data to model and create short-term forecasts for environmentally-mediated diseases. 

For more information, please see the demo project based on malaria in Ethiopia (with demo data): https://github.com/EcoGRAPH/epidemiar-demo

EPIDEMIA project: http://ecograph.net/epidemia/

