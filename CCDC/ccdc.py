import pandas as pd
import numpy as np
import sys
import xarray as xr
import argparse
import csv
import os
import fnmatch
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
from makeModel import MakeCCDCModel
from datetime import datetime
from scipy.interpolate import interp1d
from multiprocessing import Pool
from multiprocessing import Manager

np.set_printoptions(precision=4)
matplotlib.rcParams.update({'font.size': 16})

# Set up list that can be accessed by multiple processes at once
manager = Manager()
rows = manager.list()

def addChangeMarker(num_bands, change_date, obs_data, axs, model_list):

    """ Adds vertical lines to each plot every time change is detected """
       
    for i in range(num_bands):
        y_min = np.amin(obs_data[:,i+1])
        y_max = np.amax(obs_data[:,i+1])

        axs[i, 0].plot([change_date, change_date], [y_min, y_max], 'r', linewidth=2)

        interp = interp1d(model_list[i].datetimes, model_list[i].predicted, kind='cubic')
        xnew = np.linspace(model_list[i].getMinDate(), model_list[i].getMaxDate(), 500)
        
        axs[i, 0].plot(xnew, interp(xnew), 'm-', linewidth=2)
        
def setupModels(all_band_data, num_bands, init_obs, cv, alpha, bands, model_list):
    
    """Creates a model for each band and stores it in model_list"""
    
    # Get column of datetime values
    datetimes = all_band_data[:,0]
    
    # Create a model for each band and store it in model_list
    for i in range(1, all_band_data.shape[1]):
        
        band_data = all_band_data[:,i]
   
        ccdc_model = MakeCCDCModel(datetimes, init_obs, bands[i-1])
            
        ccdc_model.fitModel(band_data, cv, alpha)
        
        model_list[i-1] = ccdc_model

def getNumYears(date_list):

    """Get number of years spanned by the dataset (from Python/Rata Die date)"""
    
    try:
        last_date = datetime.fromordinal(np.amax(date_list).astype(int)).strftime('%Y')
        first_date = datetime.fromordinal(np.amin(date_list).astype(int)).strftime('%Y')
            
        num_years = int(last_date) - int(first_date)
        
    except ValueError as err:
        print("ValueError when trying to find number of years covered by dataset: {}".format(err))
        print(date_list)
        
        return 0

    return num_years

def datesToNumbers(dates):
    
    dates_as_ordinal = np.array([pd.Timestamp(x).toordinal() for x in dates])
    
    return dates_as_ordinal

def transformToArray(dataset_to_transform):

    """Transforms xarray Dataset object into a Numpy array"""
    
    ds_to_array = datesToNumbers(dataset_to_transform.time.data).reshape(-1, 1)
    
    for var in dataset_to_transform.data_vars:
        ds_to_array = np.hstack((ds_to_array, dataset_to_transform[var].values.reshape(-1, 1)))
        
    # Remove NaNs and sort by datetime    
    ds_to_array = tidyData(ds_to_array)

    return ds_to_array

def setupPredictionFile(output_file, num_bands, band_names):
    
    """Creates an output CSV file for the pixel being predicted, with column names"""
           
    with open(output_file, 'w') as output:
        writer = csv.writer(output)
        
        row = []
        
        for band in band_names:
            row.append(band)
 
        writer.writerow(row) 
        
def writeOutPrediction(output_file, end_date, model_list):
    
    with open(output_file, 'a') as output:
        writer = csv.writer(output)
        row = []
        
        for model in model_list:
            prediction = model.getPrediction(end_date)[0]
            row.append(prediction)
            
        writer.writerow(row)
      
def tidyData(pixel_ts):
    
    """Takes a single pixel time series, removes NaN values, and sorts by date."""
     
    # Remove NaNs
    pixel_nan_mask = np.any(np.isnan(pixel_ts), axis=1)
    pixel_ts = pixel_ts[~pixel_nan_mask]
    
    # Sort by date
    pixel_ts = pixel_ts[np.argsort(pixel_ts[:,0])]
                                             
    return pixel_ts        

def initModel(pixel_data, num_bands, init_obs, cv, alpha, bands, model_list):

    """Finds a sequence of 6/12/18/24 consecutive clear observations without any change, to initialize the model"""

    # Subset first n clear observations for model initialisation
    curr_obs_list = pixel_data[:init_obs,]
    
    # Start off with the model uninitialized
    model_init = False
    num_iters = 0
    
    # The next observation to be added to the model to detect change
    init_end = None

    # Model initialization sequence - keeps going until a clear set of observations is found
    while(model_init == False):

        # Check if there are not enough data points left to initialise the models
        num_data_points = len(curr_obs_list)
        
        if(num_data_points < init_obs):
            #print("Could not find a period of no change for model initialization.")
            return None
    
        # Re-initialize the models
        setupModels(curr_obs_list, num_bands, init_obs, cv, alpha, bands, model_list)

        # Get total time used for model initialization
        total_time = np.max(curr_obs_list[:,0]) - np.min(curr_obs_list[:,0])
        
        total_slope_eval = 0
        total_start_eval = 0
        total_end_eval = 0
        
        # Check for change during the initialization period. We need 12 observations with no change
        for band_model in model_list: # For each model
            
            slope_val = np.absolute(band_model.coefficients[0]) / (3 * band_model.RMSE / total_time)
            total_slope_eval += slope_val
        
            start_val = np.absolute(band_model.start_val - band_model.predicted[0]) / (3 * band_model.RMSE)
            total_start_eval += start_val
            
            end_val = np.absolute(band_model.end_val - band_model.predicted[init_obs-1]) / (3 * band_model.RMSE)
            total_end_eval += end_val
 
        if((total_slope_eval / num_bands) > 1 or (total_start_eval / num_bands) > 1 or (total_end_eval / num_bands) > 1):
            num_iters += 1
            curr_obs_list = pixel_data[0+num_iters:init_obs+num_iters,:] # Shift along 1 row
        
        else:
            model_init = True
            init_end = init_obs + num_iters
            #print("Model initialized. Iterations needed: {}".format(num_iters))

    return curr_obs_list, init_end

def findChange(pixel_data, change_file, num_bands, init_obs, x, y, args, axs, model_list):
    
    """Continues to add data points to the model until either a new breakpoint is detected, or there
        are not enough observations remaining."""
        
    # this is where all of the csv output gets written
    # needs to now go to the shared list
    global rows
    
    try:
        model_data, next_obs = initModel(pixel_data, num_bands, init_obs, args.cross_validate, args.alpha, args.bands, model_list)
    except TypeError:
        return []
    
    if(args.output_mode == "normal" and args.outtype == 'csv'):
        model_output = [[x, y, m.band, m.getMinDate(), m.getMaxDate(), m.start_val, m.end_val, ["{:0.5f}".format(c) for c in m.coefficients], m.RMSE, m.lasso_model.intercept_, m.alpha] for m in model_list]

    # Detect change
    change_flag = 0
    change_time = None
    change_mags = None
    
    prev_date = None

    while((next_obs+1) <= len(pixel_data)):

        change_eval = 0
        
        # Get next row from the time series
        new_obs = pixel_data[next_obs,]
        
        # Get date
        new_date = new_obs[0]
        
        if not prev_date:
            prev_date = new_date
        
        # Calculate difference between real and predicted value for the new observation in each band
        residuals = [new_obs[i] - model_list[i-1].getPrediction(new_date)[0] for i in range(1, num_bands+1)]
        
        for i, residual in enumerate(residuals):
            rval = np.absolute(residual) / (2 * model_list[i].getRMSE(new_date)) # RMSE will be seasonally adjusted if there are >=24 observations
            change_eval += rval            
            
        if((change_eval / num_bands) <= 1): # No deviation from model detected
            #print("Adding new data point")
            model_data = np.append(model_data, [new_obs], axis=0)
            
            date_diff = new_date - prev_date # Will be 0 on first pass
            
            if(date_diff >= args.re_init):
                setupModels(model_data, num_bands, init_obs, args.cross_validate, args.alpha, args.bands, model_list)
                
                if(args.output_mode == "normal" and args.outtype == 'csv'):
                        model_output = [[x, y, m.band, m.getMinDate(), m.getMaxDate(), m.start_val, m.end_val, ["{:0.5f}".format(c) for c in m.coefficients], m.RMSE, m.lasso_model.intercept_, m.alpha] for m in model_list]
                
                prev_date = new_date
                
            change_flag = 0 # Reset change flag because we have an inlier

        else: # Deviation from model detected
            change_flag += 1 # Don't add the new pixel to the model

            if(change_flag == 1): # If this is the first observed possible change point
                change_time = new_date # Log this as the new possible date of change
                change_mags = residuals # Log residuals as new possible magnitudes of change
    
        if(change_flag == 6):
            #print("Change detected!")
            
            if(args.output_mode == "normal"):
                
                if(args.outtype == 'plot'):
                    addChangeMarker(num_bands, change_time, pixel_data, axs, model_list)
    
                else:  
                    for model_ix in range(len(model_output)):
                        model_output[model_ix].append(change_time)
                        model_output[model_ix].append(change_mags[model_ix])
                        rows.append(model_output[model_ix])
                                               
            return pixel_data[next_obs-5:,] # Return index of date when change was first flagged
        
        # Need to get the next observation
        next_obs += 1
    
    # No change detected, end of data reached

    # Write model details out to file if needed
    if(args.output_mode == "normal" and args.outtype == 'csv'):
                        
        for model_ix in range(len(model_output)):
            rows.append(model_output[model_ix])
            
    return []
    
def runCCDC(input_data, num_bands, x_val, y_val, args):

    """The main function which runs the CCDC algorithm. Loops until there are not enough observations
        left after a breakpoint to attempt to initialize a new model."""

    # The algorithm needs at least 12 observations to start
    if(len(input_data) >= 12):

        output_file = None
        axs = None
        model_list = [None for i in range(num_bands)] # Set up list of models
        
        if(args.output_mode == "normal"):
                              
            if(args.outtype == 'plot'):

                output_file = os.path.join(args.outdir, "{}_{}.png".format(x_val, y_val))
                  
                fig, axs = plt.subplots(num_bands, 1, sharex=True, figsize=(20, 10), squeeze=False)
                
                # Set up plots with original data and screened data
                for i in range(num_bands):
                    axs[i, 0].plot(input_data[:,0], input_data[:,i+1], 'o', color='k', label='Original data', markersize=4)                    
                    axs[i, 0].set_ylabel(args.bands[i], fontdict={"size": 18})
    
        if(args.output_mode == "predictive"):
            
            # Convert stopping date to ordinal so that it can easily be predicted
            end_date = datetime.strptime(args.date_to_predict, "%Y-%m-%d").toordinal()

            # Set up output file
            output_file = os.path.join(args.outdir, "{}_{}.csv".format(x_val, y_val))
            setupPredictionFile(output_file, num_bands, args.bands)                 
        
        # Get total number of clear observations for this pixel
        num_clear_obs = len(input_data)
        
        # Decide window size - this is based on a minimum of 6 obs, plus 6 to detect change
        # A larger window size means a more complex model will be fitted
        
        if(num_clear_obs >=12 and num_clear_obs < 18):
            # Use simple model with initialization period of 6 obs
            window = 6
        
        elif(num_clear_obs >= 18 and num_clear_obs < 24):
            # Use simple model with initialization period of 12 obs
            window = 12

        elif(num_clear_obs >= 24 and num_clear_obs < 30):
            # Use advanced model with initialization period of 18 obs
           window = 18
        
        elif(num_clear_obs >= 30):
            # Use full model with initialisation period of 24 obs
            window = 24
        
        # Remaining data length needs to be smaller than window size
        while(len(input_data) >= window):

            input_data = findChange(input_data, output_file, num_bands, window, x_val, y_val, args, axs, model_list)
                    
            if(args.output_mode == "predictive"):
                
                # input_data always contains the next set of values, i.e. after a break. If the first observation in this
                # set is greater than the date to predict, the current set of models are the ones to use.
                if(len(input_data) > 0): # If there is still data left to process                        
                    if(input_data[0][0] > end_date):
                        writeOutPrediction(output_file, end_date, model_list)
                        return
                        
                else: # End of data has been reached without finding the date to predict
                    writeOutPrediction(output_file, end_date, model_list)
                
        if(args.output_mode == "normal"):
            if(args.outtype == 'plot'):

                for i in range(num_bands):
                    interp = interp1d(model_list[i].datetimes, model_list[i].predicted, kind='cubic')
                    xnew = np.linspace(model_list[i].getMinDate(), model_list[i].getMaxDate(), 500)
                    axs[i, 0].plot(xnew, interp(xnew), 'm-', linewidth=2) # Plot fitted model              
    
                # Plot empty datasets so start/end of change is included in legend
                axs[0, 0].plot([], [], 'r', label='Change')
                axs[0, 0].plot([], [], 'm', label='Fitted model')
                
                axs[0, 0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=4, fancybox=True, shadow=True)
                plt.xlabel("Date", fontdict={"size": 18})
                myFmt = mdates.DateFormatter('%m/%Y') # Format dates as month/year rather than ordinal dates
                axs[0, 0].xaxis.set_major_formatter(myFmt)
                plt.tight_layout()
                plt.savefig(output_file)
                plt.close(fig)
                                                                                               
def runOnCSV(num_bands, args):
    
    global rows

    # Check for valid directory 
    if(os.path.isdir(args.csv_dir)):
        
        # For every file in the directory
        for file in os.listdir(args.csv_dir):

            if fnmatch.fnmatch(file, '*.csv'): # Check if it's a CSV file

                try:
                    csv_file_path = os.path.join(args.csv_dir, file)
        
                    ts_data = pd.read_csv(csv_file_path)
                    
                    # Remove missing data
                    ts_data = ts_data.dropna(axis=0, how='any')
                    
                    # Convert dates to ordinal
                    ts_data.datetime = datesToNumbers(ts_data.datetime)
                    
                    # Generate unique name for output file
                    uq_name = "{}_change.csv".format(file.split('/')[-1].replace('.csv', '')) 
                
                    output_file = os.path.join(args.outdir, uq_name)                    
                
                    runCCDC(ts_data.values, num_bands, 0, 0, args)       
                                  
                    # Write headers to file
                    headers = ["x", "y", "band", "start_date", "end_date", "start_val", "end_val", "coeffs", "RMSE", "intercept", "alpha", "change_date", "magnitude"]
                      
                    with open(output_file, 'w') as output:
                        writer = csv.writer(output)
                        writer.writerow(headers)
                        writer.writerows(rows)
     
                except AttributeError:
                    print("Could not process CSV file {}. Check column names".format(file))
                    
            rows = []
                    
    else:
        print("CSV directory invalid")
        
def runOnNetCDF(num_bands, args):
    
    global rows

    # Calculate the right number of columns to be returned from the data cube
    input_num_cols = num_bands + 1    

    input_data = xr.open_dataset(args.input_file)
    
    ccdc_args = []
                
    # We want to process each pixel seperately
    for i in range(len(input_data.x)):
        for j in range(len(input_data.y)):

            input_ts = input_data.isel(x=i, y=j) # Get just one pixel
            
            x_val = float(input_ts.x)
            y_val = float(input_ts.y)

            input_ts = transformToArray(input_ts) # Transform the time series into a numpy array                    
                                      
            if(input_ts.shape[0] > 0 and input_ts.shape[1] == input_num_cols):
                
                argslist = (input_ts, num_bands, x_val, y_val, args)
                ccdc_args.append(argslist)
                                  
    # Do some tidying up
    del input_data   
  
    # Run processes for this key                                        
    with Pool(processes=args.num_procs) as pool:
        pool.starmap(runCCDC, ccdc_args)

    # Generate output file name for this key
    output_file = os.path.join(args.outdir, "{}.csv".format(args.output_file)) 

    # Write headers to file
    headers = ["x", "y", "band", "start_date", "end_date", "start_val", "end_val", "coeffs", "RMSE", "intercept", "alpha", "change_date", "magnitude"]
  
    with open(output_file, 'w') as output:
        writer = csv.writer(output)
        writer.writerow(headers)
        writer.writerows(rows)

    # Reset shared list
    rows = []                              
      
def main(args):
    
    """Program runs from here"""

    os.makedirs(os.path.dirname(args.outdir), exist_ok=True)
       
    if(not args.csv_dir and not args.input_file):
        print("Either a directory containing CSV files or a path to a NetCDF file is required for input.")
        sys.exit()
       
    num_bands = len(args.bands)
    
    # Check output mode
    if(args.output_mode == "predictive"): # Predictive mode requires some extra parameters
        if(args.date_to_predict):
        
            new_dir = os.path.join(args.outdir, args.date_to_predict)

            if not(os.path.isdir(new_dir)):
                os.makedirs(new_dir)
                
            args.outdir = new_dir
        
        else:
            print("Date to predict must be specified if output mode is predictive.")
            sys.exit()
               
    if(args.process_mode == "csv" and args.csv_dir is not None):
        runOnCSV(num_bands, args)
        
    elif(args.process_mode == "nc" and args.input_file is not None): # NetCDF file provided
        runOnNetCDF(num_bands, args)
        
    else:
        print("You did not provide the correct process mode, input directory, or input file.")
    

if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description="Run CCDC algorithm using Data Cube.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--outdir', default="./", help="The output directory for any produced csv files/images/plots.")
    parser.add_argument('-pm', '--process_mode', choices=['csv', 'nc'], default='all', help="Whether the algorithm should be run on a specified area, a subsample of a (specified) area, a specific tile, or all available data. You can also provide a CSV file or NetCDF file to analyse.")
    parser.add_argument('-om', '--output_mode', choices=['normal','predictive'], default="normal", help="Whether the algorithm should generate change output (normal) or output a prediction for the area specified.")
    parser.add_argument('-pdate', '--date_to_predict', help="The date to predict for, if output_mode is predictive. Must be in format YYYY-MM-DD")  
    parser.add_argument('-ot', '--outtype', choices=['plot', 'csv'], default='csv', help="Specifies the format of the output data. Either a plot or a CSV file will be produced for each pixel.")
    parser.add_argument('-i', '--re_init', type=int, default=0, help="The number of days that should pass before the model is refitted. By default the model is refitted after every new observation.")
    parser.add_argument('-p', '--num_procs', type=int, default=1, help="The number of processes to use.")
    parser.add_argument('-b', '--bands', nargs='+', required=True, help="List of band names to use in the analysis.")    
    parser.add_argument('-csv', '--csv_dir', help="The directory to search for CSV files, if process mode is CSV.")
    parser.add_argument('-a', '--alpha', type=float, default=1.0, help="Alpha parameter for Lasso regression. Defaults to 1.0.")
    parser.add_argument('-cv', '--cross_validate', default=False, action='store_true', help="Whether to use cross validation to find alpha parameter. If set the --alpha argument will be ignored.")
    parser.add_argument('-if', '--input_file', help="Path to input NetCDF file, if process mode is nc.")
    parser.add_argument('-of', '--output_file', help="The output file name to use if using CSV output, e.g. tile number.")

    args = parser.parse_args()
      
    main(args)








