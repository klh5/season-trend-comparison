############################################################################################
### Exponentially Weighted Moving Average Change Detection (EWMACD) for R
### v 1.3.0
### Author: Evan B. Brooks
### Author email: evbrooks@vt.edu

### Citation: Brooks, E.B., R.H. Wynne, V.A. Thomas, C.E. Blinn, and J.W. Coulston. 2014. 
#### On-the-fly massively multitemporal change detection using statistical quality control 
#### charts and Landsat data. IEEE Transactions on Geosciences and Remote Sensing 
#### 52(6):3316-3332. doi: dx.doi.org/10.1109/TGRS.2013.2272545

### November 2018
### Edited by Katie Awty-Carroll @ Aberystwyth University (klh5@aber.ac.uk) for use with CSV files

############################################################################################

library(zoo)
library(stringr)
library(ggplot2)
library(scales)

sim_dir = "/path/to/inputs" # Directory containing CSV files for input
output_dir = "./ewmacd" # Directory for outputs

###################################################################
###### EWMACD functions for use with the rasterEngine paradigm

## Base function (pixel-level), identical to base function for rasterEngine-based HREG except that this one does not incorporate additional arguments explicitly
EWMACD.pixel.for.rasterEngine=function(myPixel,datetime,nc,ns,DOYs,Years,historybound,xBarLimit1,xBarLimit2,lowthresh,lambdasigs,lambda,rounding,persistence,trainingStart,trainingEnd,testingEnd,...){
  
  myPixel=c(myPixel) ### Ensuring vector format
  Dates=length(myPixel) ### Convenience object
  tmp=rep(-2222, Dates) ### Coded 'No data' output, fills the output as an initial value
  Beta=cbind(rep(NA,(ns+nc+1))) ### Coded other 'No data' output for the coefficients
  tmp2=-4 # Dummy value
  myPixel00=myPixel ### Backup value for myPixel
  
  ind00=c(1:Dates) ### Index list for original data
  myPixel01=myPixel[1:historybound] ### Training data
  myPixel02=myPixel[(historybound+1):Dates] ### Testing data
  
  bkgd.ind00=which(is.na(myPixel00)==F) ### Index for all non-missing data
  myPixel0=myPixel[bkgd.ind00] ### 'Present' data
  Dates00=length(myPixel0) ### Convenience object for number of dates
  bkgd.ind01=which(is.na(myPixel00)==F & ind00<=historybound) ### Index for non-missing training data
  historybound01=length(bkgd.ind01) ### Adjustment of training cutoff to reflect present data only
  
  myPixel1=myPixel00[bkgd.ind01] ### Present training data
  timedat01=DOYs[bkgd.ind01] ### Present training dates, note the implicit dependence on DOYS
  timedat1=timedat01*2*pi/365 ### Conversion of training dates only to [0,2pi]
  timedatALL=DOYs[bkgd.ind00]*2*pi/365 ### Conversion of all present dates to [0,2pi]

  if(length(myPixel1)>0){ ### Checking if there is data to work with...
    
    ### Harmonic regression component
    X=cbind(rep(1,length(timedat1)),sin(t(matrix(rep(c(1:ns),length(timedat1)),ncol=length(timedat1)))*timedat1),cos(t(matrix(rep(c(1:nc),length(timedat1)),ncol=length(timedat1)))*timedat1))
    
    XAll=cbind(rep(1,Dates00),sin(t(matrix(rep(c(1:ns),Dates00),ncol=Dates00))*timedatALL),cos(t(matrix(rep(c(1:nc),Dates00),ncol=Dates00))*timedatALL))
    
    if(length(myPixel1)>(ns+nc+1) & abs(det(t(X)%*%X))>=0.001){ # Ensuring design matrix is sufficient rank and nonsingular
      
      Preds1 = (X%*%solve(t(X)%*%X)%*%t(X)%*%cbind(c(myPixel1)))
      
      ## Block for X-bar chart anomaly filtering
      Resids1 <<- myPixel1-Preds1 # Calculate residuals
      std=sd(Resids1) # Calculate SD
      screen1=(abs(Resids1)>(xBarLimit1*std))+0
      keeps=which(screen1==0)
      
      if(length(keeps)>(nc+ns+1)) {
        Beta = solve(t(X[keeps,])%*%X[keeps,])%*%t(X[keeps,])%*%myPixel0[keeps]
      }
    }
    
    ### EWMA component
    if(is.na(Beta[1])==F) { ### Checking for present Beta
      
      y0 <<- as.numeric(myPixel0-t(XAll%*%Beta)) ### Residuals for all present data, based on training coefficients
      
      # Plot fitted model
      #output_df <- data.frame(fitted=XAll%*%Beta, actual=na.omit(myPixel), real_dates=as.Date(datetime[!is.na(myPixel)]))
      #output_plt <<- ggplot(output_df, aes(real_dates, fitted, group=1)) + geom_path(color='blue') + 
                    #geom_point(aes(real_dates, actual), stroke=0, alpha=0.5, size=2)

      y01=y0[1:historybound01] ### Training residuals only

      y02=c() 
      
      if(length(y0)>length(y01)){
        y02=y0[(historybound01+1):length(y0)]
      } ### Testing residuals
      
      mu=mean(y01) ### First estimate of historical mean (should be near 0)
      histsd=sd(y01) ### First estimate of historical SD.  
      ind0=c(1:length(y0)) ### Index for residuals
      ind01=ind0[1:historybound01] ### Index for training residuals
      
      ind02=c()
      
      if(length(y0)>length(y01)){
        ind02=ind0[(historybound01+1):length(y0)]
      } ### Index for testing residuals
      
      ### Creating date information in linear form (days from a starting point instead of Julian days of the year)
      eaYear=c(0,(rep(365,length(c(trainingStart:(testingEnd-1))))+1*(c(trainingStart:(testingEnd-1))%%4==0))) 
      cuYear=cumsum(eaYear)
      x0=(cuYear[Years-trainingStart+1]+DOYs)[bkgd.ind00] ### Compare this to the DateMaker function we wrote

      ### Modifying SD estimates based on anomalous readings in the training data
      UCL0=c(rep(xBarLimit1,length(ind01)),rep(xBarLimit2,length(ind02)))*histsd ### Note that we don't want to filter out the changes in the testing data, so xBarLimit2 is much larger!
      x=x0[myPixel0>lowthresh & abs(y0)<UCL0] ### Keeping only dates for which we have some vegetation and aren't anomalously far from 0 in the residuals
      y=y0[myPixel0>lowthresh & abs(y0)<UCL0] ### Keeping only dates for which we have some vegetation and aren't anomalously far from 0 in the residuals
      ind=ind0[myPixel0>lowthresh & abs(y0)<UCL0] ### Keeping only dates for which we have some vegetation and aren't anomalously far from 0 in the residuals

      histsd=sd(y01[which(myPixel1>lowthresh & abs(y01)<UCL0[1:historybound01])]) ### Updating the training SD estimate.  This is the all-important driver of the EWMA control limits.

      if(is.na(histsd)==1){
        return(tmp)
        break
      }
      
      Totals=y0*0 ### Future EWMA output
      tmp2=rep(-2222,length(y)) ### Coded values for the 'present' subset of the data
      
      ewma=y[1] ### Initialize the EWMA outputs with the first present residual

      for(i in (2):length(y)){ 
        ewma=c(ewma,(ewma[(i-1)]*(1-lambda)+lambda*y[i])) ### Appending new EWMA values for all present data.
      }
      
      UCL=histsd*lambdasigs*sqrt(lambda/(2-lambda)*(1-(1-lambda)^(2*c(1:length(y))))) ### EWMA upper control limit.  This is the threshold which dictates when the chart signals a disturbance.
      if(rounding==TRUE){tmp2=sign(ewma)*floor(abs(ewma/UCL))} ### Integer value for EWMA output relative to control limit (rounded towards 0).  A value of +/-1 represents the weakest disturbance signal
      if(rounding==FALSE){tmp2=round(ewma,0)} ### EWMA outputs in terms of resdiual scales.  

      ###  Keeping only values for which a disturbance is sustained, using persistence as the threshold
      if(persistence>1 & length(tmp2)>3){ ### Ensuring sufficent data for tmp2
        tmpsign=sign(tmp2) # Disturbance direction
        shiftpoints=c(1,which(tmpsign[-1]!=tmpsign[-length(tmpsign)]),length(tmpsign)) # Dates for which direction changes
        
        tmp3=rep(0,length(tmpsign))
        for(i in 1:length(tmpsign)){  # Counting the consecutive dates in which directions are sustained
          tmp3lo=0
          tmp3hi=0
          
          while(((i-tmp3lo)>0)){if((tmpsign[i]-tmpsign[i-tmp3lo])==0){tmp3lo=tmp3lo+1} else{break}}
          while(((tmp3hi+i)<=length(tmpsign))){if((tmpsign[i+tmp3hi]-tmpsign[i])==0){tmp3hi=tmp3hi+1} else{break}}
          tmp3[i]=tmp3lo+tmp3hi-1
          
        }

        tmp4=rep(0, length(tmp3)) 
        for(i in 1:length(tmp3)){ # If sustained dates are long enough, keep; otherwise set to previous sustained state
          if(tmp3[i]>=persistence){
            tmp4[i]=tmp2[i]
          } else {
              tmp4[i] = max(tmp2[which.max(tmp3[which(tmp3[1:i] >= persistence)])], 0)
              }
        }

        tmp2=tmp4
        
      }

      tmp[bkgd.ind00[ind]]=tmp2 ### Assigning EWMA outputs for present data to the original template.  This still leaves -2222's everywhere the data was missing or filtered.
      
      if(tmp[1]==-2222){ ### If the first date of myPixel was missing/filtered, then assign the EWMA output as 0 (no disturbance).
        tmp[1]=0
      }

      if(tmp[1]!=-2222){ ### If we have EWMA information for the first date, then for each missing/filtered date in the record, fill with the last known EWMA value
        for(stepper in 2:Dates){
          if(tmp[stepper]==-2222){tmp[stepper]=tmp[stepper-1]}
        }
      }
      
    }
  }
  return(as.integer(tmp)) ### Final output.  All -2222's if data were insufficient to run the algorithm, otherwise an EWMA record of relative (rounded=T) or raw-residual (rounded=F) format.
}

get_vertices <- function(subset, dist_thresh) {
  
  # Recursive function to find list of potential vertices for re-initialisation
  
  first_last <- subset[c(1, nrow(subset)),] # Just get first and last points
  
  seg_line <- lm(mags ~ ix, data=first_last) # Fit a line between the first and last points
  
  subset['resids'] <- (subset$mags - as.numeric(predict(seg_line, subset)))**2 # Calculate squared residuals for all results in subset
  subset$resids <- round(subset$resids, 4)
  
  if(any(subset$resids > 0)) { # If residuals are all 0, there is no vertex since all points fit on the line
    
    # Get point of maximum deviation
    max_resid <- subset[subset$resids == max(subset$resids),]$ix[1]
    
    if(all(abs((max_resid-vertices)) > dist_thresh)) { # Check that vertex isn't within a certain distance of any current candidates
      
      vertices <<- append(vertices, max_resid) # Append to list of current candidates
      
      resid_ix <- which(subset$ix == max_resid)
      
      # Divide series up into two new segments 
      get_vertices(subset[1:resid_ix,], dist_thresh)
      get_vertices(subset[resid_ix:nrow(subset),], dist_thresh)
      
    }
  }
}

# Run EWMACD on simulated time series
# Get 46 observations and round up to neareast year

run_ewmacd <- function(sim_file) {
  
  min_obs <- 46
  
  persistence <- 6 # Number of times control limit has to be exceeded for a change to count as persistent - 6 matches CCDC
  
  ts_data <- read.csv(sim_file) # Read in data file
  doy_data <- read.csv("DOY.csv") # Read in DOY file - all simulations use the same file
  
  # Strip out unique name
  uq_name = gsub('.csv', '', tail(strsplit(sim_file, '/')[[1]], n=1))
  
  # Set up name for output file
  results_file <- sprintf('%s/%s%s', output_dir, uq_name, '_ewmacd.csv') 
  
  # Set up data frame for output
  results <- data.frame(start_date=as.Date(character()), end_date=as.Date(character()), RMSE=numeric(), change_date=as.Date(character()), magnitude=numeric(), ewmacd_magnitude=numeric())
  
  while(nrow(na.omit(ts_data)) > min_obs) { # Continue until not enough data is left in the time series

    start_train <- min(doy_data$Year) # Get first year in time series
    
    # Reset row indices to 0
    rownames(ts_data) <- NULL
    rownames(doy_data) <- NULL

    # Get index of min_obs non-missing observation
    end_ix <- rownames(tail(na.omit(ts_data)[1:min_obs,], 1))
    
    # Get end of training period (non-inclusive)
    end_train <- doy_data[end_ix,]$Year + 1

    # Get last year in dataset
    end_test <- max(doy_data$Year)

    # historybound represents the number of observations in the history period (including missing values)
    hb <- as.numeric(end_ix)

    # Run change detection with 2 sine and 2 cosine terms
    ewmacd_res <- EWMACD.pixel.for.rasterEngine(ts_data$band_1,ts_data$datetime,nc=2,ns=2,DOYs=doy_data$DOY,Years=doy_data$Year,historybound=hb,xBarLimit1=1.5,xBarLimit2=20,lowthresh=-1,lambdasigs=3,lambda=0.3,rounding=T,persistence=persistence,trainingStart=start_train,trainingEnd=end_train,testingEnd=end_test)

    res_df <- data.frame(datetime=as.Date(ts_data$datetime), magnitude=ewmacd_res)
    
    # Add magnitude of disturbance to plot
    #output_plt <- output_plt + geom_point(data=res_df, aes(datetime, magnitude/10, colour='red'), show.legend=FALSE) +
      #scale_y_continuous(sec.axis = sec_axis(~ . *10, name = "EWMACD change magnitude")) + theme_bw() + 
      #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            #axis.line = element_line(colour = "black"), axis.text=element_text(size=12)) + xlab('Year') + ylab('NDVI') +
      #scale_x_date(date_breaks="1 year", labels=date_format("%Y"))
    
    #print(output_plt)
      
    if(!all(ewmacd_res == -2222)) { # All -2222 values mean the algorithm failed
      
      if(any(ewmacd_res[1:hb] != 0)) { # Check for stability in training period

        # Remove one observation from the start, then try again
        ts_data <- ts_data[-1,]
        doy_data <- doy_data[-1,]

      } else {
      
        start_date <- min(as.Date(na.omit(ts_data)$datetime)) # Start date of current model
        
        if(!all(ewmacd_res == 0)) { # If a change was found
          
          change_ix <- head(which(ewmacd_res != 0), 1)       # Get index of point where EWMACD flags a change
          ewmacd_magnitude <- ewmacd_res[change_ix]         # Get value of change magnitude in control limits
          change_date <- as.character(ts_data[change_ix,]$datetime)  # Get date of change
          end_date <- as.character(ts_data[change_ix-1,]$datetime) # End date of current model if there is a change
          
          # y0 is residuals for whole model, minus NA values
          # To get the RMSE for the stable model, we need the residuals for all non-NA values up to the change point
          # To do this we need to find the index of the change date when NA's are removed
          no_na <- na.omit(ts_data)
          rownames(no_na) <- NULL # Reset index
          last_date_ix <- which(no_na$datetime == change_date)-1 # Last date before change
          rmse <- sqrt(mean(y0[1:last_date_ix]**2)) # RMSE for stable period
          
          # Magnitude is the residual value for the date of change
          magnitude <- y0[last_date_ix+1]
  
          dist_thresh <- persistence / 2 # Distance new vertices must be from previous vertices
          
          per_change <- ewmacd_res[change_ix:length(ewmacd_res)] # Subset all results from first vertex onwards
          
          # Check if all change magnitudes are identical
          # If they are, the next vertex will be the first date after the change
          if(length(unique(per_change)) == 1) {
            ts_data <- ts_data[change_ix:nrow(ts_data),]
            doy_data <- doy_data[change_ix:nrow(doy_data),]
          } else {
            
            subset <- data.frame(mags=per_change, ix=c(change_ix:length(ewmacd_res))) # Convert to data frame
            
            vertices <<- c() # List of potential second vertices
            
            # Get candidate vertices
            get_vertices(subset, dist_thresh)
            
            second_vertex <- min(vertices)
            
            ts_data <- ts_data[second_vertex:nrow(ts_data),]
            doy_data <- doy_data[second_vertex:nrow(doy_data),]
          }
  
          # Append results to data frame
          results <- rbind(results, data.frame(start_date=as.Date(start_date), end_date=as.Date(end_date), RMSE=rmse, change_date=as.Date(change_date), magnitude=magnitude, ewmacd_magnitude=ewmacd_magnitude))
        
          } else { # No change found
            
            rmse <- sqrt(mean(y0**2)) # RMSE over whole data set
            end_date <- max(as.Date(ts_data$datetime)) # End date of current model if no change
            
            # Append results to data frame with no change data
            results <- rbind(results, data.frame(start_date=as.Date(start_date), end_date=as.Date(end_date), RMSE=rmse, change_date=NA, magnitude=NA, ewmacd_magnitude=NA))
  
            break
          }
      }

    } else {
      # Algorithm failed to run
      
      # Check for any results from previous runs
      if(dim(results)[1] > 0) {
        write.csv(results, results_file)
      }
      
      break
    }
  } 
  
  # Not enough data left
  # Write out file once done
  write.csv(results, results_file)
}

# Get list of files in directory
files <- list.files(path=sim_dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)

# Run BFAST on all files
lapply(files, run_ewmacd)


