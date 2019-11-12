library(bfast)
library(zoo)
library(lubridate)
library(stringr)

sim_dir = "dir/containing/simulations" # Directory containing CSV files for analysis
output_dir = "./bfast_monitor" # Directory for output

# Minimum number of observations needed in time series as a training period
min_obs <- 46

run_bfast_monitor <- function(sim_file) {
  
  # Read in CSV file
  sim_data = read.csv(sim_file)
  
  # Remove NA values
  sim_data <- na.omit(sim_data)
  
  sim_data$datetime <- as.Date(sim_data$datetime)
  
  # Strip out unique name
  uq_name = gsub('.csv', '', tail(strsplit(sim_file, '/')[[1]], n=1))
  
  # Set up name for output file
  results_file <- sprintf('%s/%s%s', output_dir, uq_name, '_bfast_monitor.csv') 
  
  # Set up data frame for output
  results <- data.frame(start_date=as.Date(character()), end_date=as.Date(character()), rmse=numeric(), change_date=as.Date(character()), 
                        magnitude=numeric(), intercept=numeric(), trend=numeric(), cos1=numeric(), sin1=numeric(), cos2=numeric(), sin2=numeric())
  
  # Repeat monitoring until there is not enough data remaining
  # Need at least 1 observation to detect change, so anything > min_obs is fine
  while(nrow(sim_data) > min_obs) {
    
    rownames(sim_data) <- NULL
    
    # Create irregular BFAST time series
    ts_data <- bfastts(sim_data$band_1, sim_data$datetime, type=c("irregular"))
    
    # Get date of first non-missing value after the minimum number of obs required for training, which is when monitoring will start
    last_obs <- tail(head(sim_data, min_obs+1), 1)$datetime

    # run BFAST Monitor - can include NA values
    # The Reverse Order CUSUM method will find any breaks in the history period and only use data up to that point
    monitor_res = bfastmonitor(ts_data, decimal_date(last_obs), history="ROC", order=2)
    
    # Get residuals for whole data set
    all_resids <- sim_data$band_1 - monitor_res$tspp$prediction
    
    # Get date of start of history period
    # The ROC test can discount unstable values at the start of the time series
    start_date <- as.Date(date_decimal(monitor_res$history[1]))
    
    # Some dates are incorrect by 1 day when returned from a ts object - find closest real date
    start_date <- sim_data$datetime[which.min(abs(sim_data$datetime-start_date))]

    if(!is.na(monitor_res$breakpoint)) {
      
      # Get change date
      change_date <- as.Date(date_decimal(monitor_res$breakpoint))
      change_date <- sim_data$datetime[which.min(abs(sim_data$datetime-change_date))]
      
      end_date <- sim_data$datetime[which(sim_data$datetime == change_date) - 1]
      
      # Get residuals for stable period
      stable_ix <- which(sim_data$datetime >= start_date & sim_data$datetime <= end_date)
      resids <- all_resids[stable_ix]
      
      # Calculate RMSE
      rmse <- sqrt(mean(resids**2))
      
      change_mag <- monitor_res$magnitude # Get magnitude of change
      
      # Get coefficients
      intercept <- monitor_res$model$coefficients[["(Intercept)"]]
      trend <-  monitor_res$model$coefficients[["trend"]]
      cos1 <- monitor_res$model$coefficients[["harmoncos1"]]
      sin1 <- monitor_res$model$coefficients[["harmonsin1"]]
      cos2 <- monitor_res$model$coefficients[["harmoncos2"]]
      sin2 <- monitor_res$model$coefficients[["harmonsin2"]]
      
      results <- rbind(results, data.frame(start_date=start_date, end_date=end_date,rmse=rmse,change_date=change_date,
                              magnitude=change_mag, intercept=intercept, trend=trend, cos1=cos1, sin1=sin1, cos2=cos2, sin2=sin2))
      
      # Subset data after break for next iteration
      sim_data <- sim_data[as.Date(sim_data$datetime) >= change_date,]
      
    } else {
      
      # No breakpoint found, end of data reached
      
      # Get RMSE over whole data set
      rmse <- sqrt(mean(all_resids**2))
      end_date <- max(sim_data$datetime)
      
      # Get coefficients
      intercept <- monitor_res$model$coefficients[["(Intercept)"]]
      trend <-  monitor_res$model$coefficients[["trend"]]
      cos1 <- monitor_res$model$coefficients[["harmoncos1"]]
      sin1 <- monitor_res$model$coefficients[["harmonsin1"]]
      cos2 <- monitor_res$model$coefficients[["harmoncos2"]]
      sin2 <- monitor_res$model$coefficients[["harmonsin2"]]
      
      results <- rbind(results, data.frame(start_date=start_date, end_date=end_date,rmse=rmse,change_date=NA, magnitude=NA,
                                           intercept=intercept, trend=trend, cos1=cos1, sin1=sin1, cos2=cos2, sin2=sin2))
      
      break
    }
    
  }
  
  # Write out file once done
  write.csv(results, results_file)
  
}

# Get list of files in directory
files <- list.files(path=sim_dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)

# Run BFAST on all files
lapply(files, run_bfast_monitor)




