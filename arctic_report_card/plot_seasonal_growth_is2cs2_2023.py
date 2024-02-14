
import matplotlib, sys
#matplotlib.use('Agg')
#from mpl_toolkits.basemap import Basemap
import numpy as np
from pylab import *
import numpy.ma as ma
import xarray as xr
import pandas as pd
from scipy.interpolate import griddata
from netCDF4 import Dataset
import netCDF4 as nc4
import time
import dask.array as da
import sys
import os
from glob import glob
from datetime import datetime
#import rasterio
#from rasterio.fill import fillnodata
#sys.path.append('../../ICESat-2-sea-ice-thickness-v4/Code/')
#import common_functions as cF
#cF.reset_matplotlib()

def reset_matplotlib():
    """
    Reset matplotlib to a common default.

    """
    
    # Force agg backend.
    plt.switch_backend('agg')
    # These settings must be hardcoded for running the comparision tests and
    # are not necessarily the default values.
    plt.rcParams['ytick.major.size'] = 2
    plt.rcParams['axes.linewidth'] = .6
    plt.rcParams['lines.linewidth'] = .6
    plt.rcParams['patch.linewidth'] = .6
    plt.rcParams['ytick.labelsize']=8
    plt.rcParams['xtick.labelsize']=8
    plt.rcParams['legend.fontsize']=8
    plt.rcParams['font.size']=8
    plt.rcParams['pcolor.shading']='auto'
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})


def add_time_dim_v3(xda):
    """ dummy function to just set current time as a new dimension to concat files over, change later! """
    xda = xda.set_coords(["latitude","longitude", "x", "y"])
    xda = xda.expand_dims(time = [datetime.now()])
    return xda

def read_IS2SITMOGR4(data_type='zarr-s3', version='V3', local_data_path="/data/IS2SITMOGR4/", 
    bucket_name="icesat-2-sea-ice-us-west-2", persist=True, download=False): 
    """ Read in IS2SITMOGR4 monthly gridded thickness dataset from local netcdf files, download the netcdf files from S3 storage, or read in the aggregated zarr dataset from S3. Currently supports either Version 2 (V2) or Version 3 (V3) data. 
    
    Note than in Version 3 there was a change in the xgrid/ygrid coordinates to x/y.
    
    Args: 
        data_type (str, required): (default to "zarr-s3", but also "netcdf-s3" or "netcdf-local" which is a local version of the netcdf files)
        version (str, required): IS2SITMOGR4 version (default "V2")
        local_data_path (str, required): local data directory
        bucket_name (str, required): S3 bucket name
        persist (boleen, required): if zarr option decide if you want to persist (load) data into memory
        download (boleen, required): download from s3 bucket to local storage

    Returns: 
        is2_ds (xr.Dataset): aggregated IS2SITMOGR4 xarray dataset, dask chunked/virtually allocated in the case of the zarr option (or allocated to memory if persisted). 
        
    Version History: 
        (November 2023 for V3 data release):  
            - Moved the download code to it's own section at the start of the function
            - Changed local paths
            - Baked in the date_str label as that is just a function of the dataset version anyway
            - Adapted the netcdf reader to use open_mfdataset, required a preprocessing data dimension step. Much more elegant!
    
    """
    
    if download==True:
        print("download from S3 bucket: ", bucket_name)

        # Download netCDF data files
        s3_path = 's3://'+bucket_name+'/IS2SITMOGR4_'+version+'/netcdf/'
        fs = s3fs.S3FileSystem(anon=True)

        #files references the entire bucket.
        files = fs.ls(s3_path)
        for file in files:
            print('Downloading file from bucket to local storage...', file)
            fs.download(file, local_data_path+version+'/')
        
    if data_type=='zarr-s3':
        if version=='V2':
            date_str='201811-202204'
        else:
            date_str='201811-202304'
        print('load zarr from S3 bucket: ', bucket_name)
        s3_path = 's3://'+bucket_name+'/IS2SITMOGR4_'+version+'/zarr/IS2SITMOGR4_'+version+'_'+date_str+'.zarr/all/'
        print('zarr_path:', s3_path)
        s3 = s3fs.S3FileSystem(anon=True)
        store = s3fs.S3Map(root=s3_path, s3=s3, check=False)
        is2_ds = xr.open_zarr(store=store)
        
        # Had a problem with these being loaded as dask arrays which cartopy doesnt like
        is2_ds = is2_ds.assign_coords(longitude=(["y","x"], is2_ds.longitude.values))
        is2_ds = is2_ds.assign_coords(latitude=(["y","x"], is2_ds.latitude.values))

        if persist==True:
            is2_ds = is2_ds.persist()

    elif data_type=='netcdf':
        print('Looking in...', local_data_path+version+'/netcdf/')
        filenames = glob(local_data_path+version+'/netcdf/*.nc')
        if len(filenames) == 0: 
            raise ValueError("No files, exit")
            return None
        
        dates = [pd.to_datetime(file.split("IS2SITMOGR4_01_")[1].split("_")[0], format = "%Y%m")  for file in filenames]
        # Add a dummy time then add the dates I want, seemed the easiest solution
        if version=='V2':
            is2_ds = xr.open_mfdataset(filenames, preprocess = add_time_dim_v2, engine='netcdf4')
        else:
            is2_ds = xr.open_mfdataset(filenames, preprocess = add_time_dim_v3, engine='netcdf4')
                
        is2_ds["time"] = dates

        # Sort by time as glob file list wasn't!
        is2_ds = is2_ds.sortby("time")
        if version=='V2':
            is2_ds = is2_ds.set_coords(["latitude","longitude","xgrid","ygrid"]) 
        else:
            is2_ds = is2_ds.set_coords(["latitude","longitude","x","y"])
        
        is2_ds = is2_ds.assign_coords(longitude=(["y","x"], is2_ds.longitude.values))
        is2_ds = is2_ds.assign_coords(latitude=(["y","x"], is2_ds.latitude.values))
        
        is2_ds = is2_ds.assign_attrs(description="Aggregated IS2SITMOGR4 "+version+" dataset.")

    
    return is2_ds


def getWinterDateRangev3(start_year, end_year): 
    """ Gets date range for winter season/s
    Args: 
        start_year (int): start year 
        end_year (int): end year 
        
    Returns: 
        winters (list): list of dates for all winter seasons in the input range (i.e: ['1980-11','1980-12','1981-01',
         '1981-02','1981-03','1981-04')
    """
    winters = []
    winters += pd.date_range(start = str(2018) + '-11', end = str(2019) + '-04', freq = 'MS')
    winters += pd.date_range(start = str(2019) + '-10', end = str(2020) + '-04', freq = 'MS')
    winters += pd.date_range(start = str(2020) + '-10', end = str(2021) + '-04', freq = 'MS')
    winters += pd.date_range(start = str(2021) + '-10', end = str(2022) + '-04', freq = 'MS')
    winters += pd.date_range(start = str(2022) + '-10', end = str(2023) + '-04', freq = 'MS')
        
    return winters


def calcMeans(dataset, extraStr='_int'): 
    """Calculates mean monthly uncertainity and mean monthly ice thicknesses for all ice types combined, multi year ice, and first year ice, and adds data as coordinates to the dataset.
    
    Args: 
        dataset (xr Dataset): dataset generated by Load_IS2 notebook
    
    Returns: 
        datasetMeans (xr Dataset): dataset with means added as data variables to the input dataset
    """  

    #calculate mean thickness for all ice types
    meanThickness = dataset['ice_thickness'+extraStr].mean(dim = ['x','y'], skipna = True)
    meanThickness.attrs = {'description': 'mean monthly ice thickness for all ice types', 'units': 'meters'}

    return meanThickness


def get_cs2():
    cs2 = pd.read_csv(cs2DataPath+'cs2smos-timeseries-v205-monthly-is2domain.csv')
    #cs2['time']= pd.to_datetime(cs2['time'], format='%Y-%m-%d')
    cs2_dt = pd.to_datetime(cs2['time'], format='%Y-%m-%d')
    cs2=cs2.to_xarray()
    cs2=cs2.mean_sea_ice_thickness
    cs2=cs2.assign_coords(index = cs2_dt.values)
    cs2=cs2.rename({'index': 'time'})
    return cs2


reset_matplotlib()

#----IS2

version='V3'
#ancDataPath = '/Users/akpetty/Data/raw_data/OTHER/'
is2DataPath='/Users/akpetty/Data/ICESat-2/IS2SITMOGR4/'
cs2DataPath='/Users/akpetty/Data/CS2/AWI_SMOS/'
figPath='../figures/'

startYear = 2018
endYear = 2023
winters = getWinterDateRangev3(startYear, endYear) #get date range for winter 18-19 and winter 19-20

IS2SITMOGR4_v3 = read_IS2SITMOGR4(data_type='netcdf', 
                                   local_data_path=is2DataPath, version=version, 
                                  download=False, persist=True) 

innerArctic = [1,2,3,4,5,6]
regionsText = 'Inner_Arctic'
IS2SITMOGR4_v3_innerArctic = IS2SITMOGR4_v3.where(IS2SITMOGR4_v3.region_mask.isin(innerArctic))


int=True
if int==True:
    var='ice_thickness_int'
    is2means = calcMeans(IS2SITMOGR4_v3_innerArctic, extraStr='_int')
else:
    var='ice_thickness'
    is2means = calcMeans(IS2SITMOGR4_v3_innerArctic, extraStr='')

cs2=get_cs2()


#initialize figure
fig = plt.figure(figsize=(5, 3))
ax=plt.gca()

colors=['y', '#a6cee3', '#b2df8a', '#beaed4', 'k']
alphas=[0.5, 0.5, 0.5, 0.5, 0.8]
num_winters=5
linewidths=[0.4, 0.4, 0.4, 0.4, 1.0]
# IS-2
for x in range(num_winters):
    print(x)
    winterStr =  str(startYear+x)[0:4] + '-' + str(startYear + x + 1)[0:4]
    winterMonths = [2, 3,4,5,6,7,8]
    winterMonthsStr = ['Oct', 'Nov','Dec','Jan','Feb','Mar','Apr']
    datasetWinter_cs2 = cs2.sel(time = slice('Oct ' + str(startYear+x), 'Apr' + str(startYear + x + 1)))
    ax.plot(winterMonths, datasetWinter_cs2, color = colors[x], linestyle = '-', linewidth=linewidths[x], marker = 'v', markersize=4, alpha=alphas[x], label = 'CS-2/SMOS')
    
for x in range(num_winters):
    if x==0:
        # First year starts in November for IS-2
        datasetWinter_is2 = is2means.sel(time = slice('Nov ' + str(startYear+x), 'Apr' + str(startYear + x + 1))).values
        # add Nan to represent missing Oct
        datasetWinter_is2 = np.insert(datasetWinter_is2, 0, np.nan, axis=0)
    else:
        datasetWinter_is2 = is2means.sel(time = slice('Oct ' + str(startYear+x), 'Apr' + str(startYear + x + 1))).values

    ax.plot(winterMonths, datasetWinter_is2, color = colors[x], linestyle = '--', linewidth=linewidths[x], marker = 'o', markersize=4, alpha=alphas[x], label = 'IS-2 ('+winterStr+')')
    
    datasetWinter_cs2 = cs2.sel(time = slice('Oct ' + str(startYear+x), 'Apr' + str(startYear + x + 1)))
    ax.fill_between(winterMonths, 
    np.amin([datasetWinter_cs2, datasetWinter_is2], axis=0), 
    np.amax([datasetWinter_cs2, datasetWinter_is2], axis=0), 
    facecolor = colors[x], alpha = 0.35, edgecolor = 'none')

ax.set_xticks([2, 3,4,5,6,7,8])
ax.set_xticklabels(['Oct','Nov','Dec','Jan','Feb','Mar','Apr']
)

#add legend & labels
ax.legend(loc =2, ncol =2, fontsize = 7, frameon=False)
ax.set_ylabel('Sea ice thickness (m)')

ax.set_ylim([0.6, 2.4 ])
ax.set_yticks(np.arange(0.5, 2.6, 0.5))


#ax.set_xlabel('Month')
    
plt.subplots_adjust(bottom=0.09, left=0.09, top = 0.98, right=0.98)

plt.savefig(figPath+'seasonal_growth_is2_cs2_'+str(startYear)+'-'+str(endYear)+var+regionsText+'.pdf', dpi = 300)


