library(bfast)
library(zoo)
library(lubridate)
library(stringr)

sim_dir = "dir/containing/simulations" # Location of directory containing CSV files for analysis
output_dir = "./bfast" # Output directory

run_bfast <- function(sim_file) {
  
  # Read in CSV file
  sim_data = read.csv(sim_file)
  
  # Strip out unique name
  uq_name = gsub('.csv', '', tail(strsplit(sim_file, '/')[[1]], n=1))
  
  output_file <- sprintf('%s/%s%s', output_dir, uq_name, '_bfast.csv') # File for results
  
  # Remove any NA values at the start and end of the time series since they can't be reliably interpolated
  # Start
  real_start <- as.numeric(rownames(head(na.omit(sim_data), 1)))
  sim_data <- sim_data[real_start:nrow(sim_data),]
  rownames(sim_data) <- NULL
  
  # End
  real_end <- as.numeric(rownames(tail(na.omit(sim_data), 1)))
  sim_data <- sim_data[1:real_end,]
  rownames(sim_data) <- NULL
  
  # Get first and last years in data
  start_year <- as.numeric(strftime(min(as.Date(sim_data$datetime)), "%Y"))
  end_year <- as.numeric(strftime(max(as.Date(sim_data$datetime)), "%Y"))
  
  # Calculate start of time series in ts units
  start_doy <- as.numeric(strftime(sim_data$datetime[1], format = "%j")) # Convert first date to DOY
  start_ob_count <- ceiling(start_doy / 16) # Divide by 16 and round up to get number of 16-day intervals into the year
  
  # Create time series object 
  ts_data <- ts(sim_data$band_1, start=c(start_year, start_ob_count), frequency=23)
  
  # Interpolate any missing values (linear interpolation) - needed to run BFAST
  interpd_ts <- na.approx(ts_data)
  
  h_dist <- 46 / length(interpd_ts) # Allow at least 2 years between breaks

  # Run BFAST decomposition
  bfast_res <- bfast(interpd_ts, season="harmonic", type="OLS-MOSUM", h=h_dist, max.iter=5)
  
  # Get number of iterations needed to find breaks
  num_iters <- length(bfast_res$output)
  
  if(bfast_res$nobp$Vt && bfast_res$nobp$Wt) { # If no breaks were found 
    
    start_date <- as.character(sim_data$datetime[1])
    end_date <- as.character(sim_data$datetime[nrow(sim_data)])
    
    # Get trend
    trend <- diff(bfast_res$output[[num_iters]]$Tt[1:2]) / diff(time(bfast_res$output[[num_iters]]$Tt)[1:2])
    
    output_row <- data.frame(start_date=start_date, end_date=end_date, change_date=NA, trend=trend, magnitude=NA, break_type=NA)
    
    write.csv(output_row, output_file, row.names=FALSE)
    
  } else { # Some breaks were found
    
    # Set up data frame for results
    all_results <- data.frame(start_date=character(), end_date=character(), change_date=character(), trend=numeric(), magnitude=numeric(), break_type=character(), stringsAsFactors=FALSE)
    
    # Get trend output
    trend_output <- bfast_res$output[[num_iters]]$bp.Vt
    
    # If breaks in trend were found, get change dates
    # bfast_res$nobp$Vt tells us if there were no breakpoints detected in the trend component
    if(!bfast_res$nobp$Vt) {
      
      # BFAST records breaks as being the last date before the current model no longer fits
      # Need to add 1 to make this comparable to the other methods
      breakpoints <- trend_output$breakpoints + 1
      
      change_dates <- as.character(sim_data[breakpoints,]$datetime)
      
      # Get magnitudes of trend breaks
      break_mags <- bfast_res$Mags[,3]
      
      curr_start <- 1 
      start_dates <- c()
      end_dates <- c()
      trends <- c()
      
      for(t_break in breakpoints) {
        
        start_dates <- append(start_dates,  as.character(sim_data$datetime[curr_start]))
        end_dates <- append(end_dates, as.character(sim_data$datetime[t_break-1]))
        
        # Get trend
        trend_end <- curr_start+1
        trend <- diff(bfast_res$output[[num_iters]]$Tt[curr_start:trend_end]) / diff(time(bfast_res$output[[num_iters]]$Tt)[curr_start:trend_end])
        
        trends <- append(trends, trend)
        
        curr_start <- t_break
      }
      
      # Last segment
      start_dates <- append(start_dates, as.character(sim_data$datetime[curr_start]))
      end_dates <- append(end_dates, as.character(sim_data$datetime[nrow(sim_data)]))
      
      change_dates <- append(change_dates, NA)
      break_mags <- append(break_mags, NA)
      
      trend_end <- curr_start+1
      trend <- diff(bfast_res$output[[num_iters]]$Tt[curr_start:trend_end]) / diff(time(bfast_res$output[[num_iters]]$Tt)[curr_start:trend_end])
      
      trends <- append(trends, trend)
      
      # Vector for break type
      break_types <- rep('trend', length(start_dates))
      
      # Set up data frame for trend results
      break_details <- data.frame(start_date=start_dates, end_date=end_dates, change_date=change_dates, trend=trends, magnitude=break_mags, break_type=break_types)
      
      # Append trend results to overall results
      all_results <- rbind(all_results, break_details)
    }
    
    # Get season output
    
    # Needs to also output start/end dates
    
    season_output <- bfast_res$output[[num_iters]]$bp.Wt
    
    # If breaks in season were found, get change dates
    # bfast_res$nobp$Wt tells us if there were no breakpoints detected in the trend component
    if(!bfast_res$nobp$Wt) {
      
      breakpoints <- season_output$breakpoints
      change_dates <- as.character(sim_data[breakpoints,]$datetime)
      
      curr_start <- 1
      start_dates <- c()
      end_dates <- c()
      
      for(t_break in breakpoints) {
        
        start_dates <- append(start_dates,  as.character(sim_data$datetime[curr_start]))
        end_dates <- c(end_dates, as.character(sim_data$datetime[t_break-1]))
        
        curr_start <- t_break
      }
      
      # Last segment
      start_dates <- append(start_dates, as.character(sim_data$datetime[curr_start]))
      end_dates <- append(end_dates, as.character(sim_data$datetime[nrow(sim_data)]))
      
      change_dates <- append(change_dates, NA)
      
      # Vector for break type
      break_types <- rep('season', length(start_dates))
      
      # Vector for change magnitude - seasonal breaks have no magnitude
      break_mags <- rep(NA, length(start_dates))
      
      # Vector for trend - seasonal breaks have no trend
      trends <- rep(NA, length(start_dates))
      
      # Set up data frame for seasonal results
      break_details <- data.frame(start_date=start_dates, end_date=end_dates, change_date=change_dates, trend=trends, magnitude=break_mags, break_type=break_types)
    
      # Append seasonal results to overall results
      all_results <- rbind(all_results, break_details)
    }
    
    # Output results
    write.csv(all_results, output_file, row.names=FALSE)
  }
}

# Get list of files in directory
files <- list.files(path=sim_dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)

# Run BFAST on all files
lapply(files, run_bfast)
