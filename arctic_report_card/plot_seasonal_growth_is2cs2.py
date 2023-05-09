
import matplotlib, sys
#matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
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
#import rasterio
#from rasterio.fill import fillnodata
sys.path.append('../')
import common_functions as cF
cF.reset_matplotlib()


def getWinterDateRange(start_year, end_year): 
    """ Gets date range for winter season/s
    Args: 
        start_year (int): start year 
        end_year (int): end year 
        
    Returns: 
        winters (list): list of dates for all winter seasons in the input range (i.e: ['1980-11','1980-12','1981-01',
         '1981-02','1981-03','1981-04')
    """
    winters = []
    for year in range(start_year, end_year, 1):
        winters += pd.date_range(start = str(year) + '-11', end = str(year + 1) + '-04', freq = 'MS')
    return winters

def getWinterDateRangev2(start_year, end_year): 
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
        
    return winters

def getIS2Data(dataPath, dates, atl10rel='004', version='002'): 
    """Gets ICESat-2 data for provided date range
    
    Args: 
        dataPath (str): path to local directory of ICESat-2 data 
        dates (list): pandas Timestamp objects generated by getWinterDateRange
        
    Returns: 
        is2 (xarray dataset): ICESat-2 data or NONE if file does not exist for inputted date range
    """
    is2List = [] #empty list for compiling xarray DataArray objects
    for date in dates: 
        try:
            print(dataPath+ 'IS2SITMOGR4_01_*' + date.strftime('%y') + date.strftime('%m') +'_'+ atl10rel+'_'+version+'.nc')
            filename = glob(dataPath+ 'IS2SITMOGR4_01_*' + date.strftime('%y') + date.strftime('%m') +'_'+ atl10rel+'_'+version+'.nc')[0]
        except: 
            print('Cannot find files; check date range, filepath or if glob is imported')
            return None
        is2 = xr.open_dataset(filename)
        is2 = is2.assign_coords({'time': date})
        is2List.append(is2)
    is2Data = xr.concat(is2List, dim = 'time') #concatenate all DataArray objects into a single DataArray
    
    return is2Data

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

def calcMeans_old(dataset, extraStr='_int'): 
    """Calculates mean monthly uncertainity and mean monthly ice thicknesses for all ice types combined, multi year ice, and first year ice, and adds data as coordinates to the dataset.
    
    Args: 
        dataset (xr Dataset): dataset generated by Load_IS2 notebook
    
    Returns: 
        datasetMeans (xr Dataset): dataset with means added as data variables to the input dataset
    """  
    #calculate mean uncertainites 
    meanUnc = dataset['ice_thickness_unc'+extraStr].mean(dim = ['x','y'], skipna = True)
    meanUnc.attrs = {'description': 'mean monthly ice thickness uncertainty', 'units': 'meters'}
    
    #calculate mean thickness for all ice types
    meanThickness = dataset['ice_thickness'+extraStr].mean(dim = ['x','y'], skipna = True)
    meanThickness.attrs = {'description': 'mean monthly ice thickness for all ice types', 'units': 'meters'}

    #calculate mean freeboard for all ice types
    meanFreeboard = dataset['freeboard'].mean(dim = ['x','y'], skipna = True)
    meanFreeboard.attrs = {'description': 'mean monthly ice freeboard for all ice types', 'units': 'meters'}

    #calculate mean snow depth for all ice types
    meanSnowDepth = dataset['snow_depth'+extraStr].mean(dim = ['x','y'], skipna = True)
    meanSnowDepth.attrs = {'description': 'mean monthly snow depth for all ice types', 'units': 'meters'}
    
    #calculate mean thickness for multi year ice
    MYIThickness = dataset['ice_thickness'+extraStr].where(dataset['ice_type'+extraStr] == 1)
    meanMYIThickness = MYIThickness.mean(dim = ['x','y'], skipna = True)
    meanMYIThickness.attrs = {'description': 'mean monthly multi year ice thickness', 'units': 'meters'}
    
    #calculate mean thickness for first year ice 
    FYIThickness = dataset['ice_thickness'+extraStr].where(dataset['ice_type'+extraStr] == 0)
    meanFYIThickness = FYIThickness.mean(dim = ['x','y'], skipna = True)
    meanFYIThickness.attrs = {'description': 'mean monthly first year ice thickness', 'units': 'meters'}
    
    #add means as coordinates to dataset
    datasetMeans = dataset.assign_coords(coords = {'mean_ice_thickness_unc': meanUnc, 'mean_freeboard': meanFreeboard, 'mean_snow_depth':meanSnowDepth,
        'mean_ice_thickness': meanThickness, 'mean_MYI_thickness': meanMYIThickness, 'mean_FYI_thickness': meanFYIThickness})
    
    return datasetMeans

def getRegionMask(dataPath): 
    """Gets NSIDC region mask for map projection
    
     Args:
         dataPath (str): path to NSIDC region mask
         
     Returns: 
         shapedMask (numpy array): NSIDC arctic region mask gridded to shape [448, 304]
         shapedLons (numpy array): longitudes gridded to shape [448, 304]
         shapedLats (numpy array): latitudes gridded to shape [448, 304] 
    """ 
    gridShape = [448, 304] #shape of grid to reshape data to 
    
    regionMask = open(dataPath + '/sect_fixed_n.msk', 'rb') #open region mask 
    shapedMask = np.reshape(np.fromfile(file = regionMask, dtype='uint8'), gridShape) #reshape mask to grid shape
    
    maskLons = open(dataPath + '/psn25lons_v3.dat', 'rb') #open region mask longitudes
    maskLats = open(dataPath + '/psn25lats_v3.dat', 'rb') #open region mask latitudes
    shapedLons = np.reshape(np.fromfile(file = maskLons, dtype='<i4')/100000., gridShape) #reshape longitudes to grid shape
    shapedLats = np.reshape(np.fromfile(file = maskLats, dtype='<i4')/100000., gridShape) #reshape latitudes to grid shape

    return shapedMask, shapedLons, shapedLats

#function from regional_analysis notebook
def restrictRegionally(dataset, regionKeyList): 
    """Restrict dataset to input regions.
    
    Args: 
        dataset (xr Dataset): dataset generated by Load_IS2 notebook
        regionKeyList (list): list of region keys to restrict data to 
        
    Returns: 
        regionalDataset (xr Dataset): dataset with restricted data to input regions
    """
    
    def checkKeys(regionKeyList, regionTbl): 
        """Check that regionKeyList was defined correctly

        Raises: 
            ValueError if regionKeyList was not defined correctly 
            warning if all data was removed from the dataset
        """
        if type(regionKeyList) != list: #raise a ValueError if regionKeyList is not a list 
            raise ValueError('regionKeyList needs to be a list. \nFor example, if you want to restrict data to the Beaufort Sea, define regionKeyList = [13]')

        for key in regionKeyList: 
            if key not in list(regionTbl['key']): 
                raise ValueError('Region key ' + str(key) + ' does not exist in region mask. \n Redefine regionKeyList with key numbers from table')

        if len(regionKeyList) == 0: 
            warnings.warn('You removed all the data from the dataset. Are you sure you wanted to do this? \n If not, make sure the list regionKeyList is not empty and try again. \n If you intended to keep data from all regions, set regionKeyList = list(tbl[\"key\"])')
 
    #create a table of keys and labels
    regionMask = dataset.region_mask.attrs
    regionTbl = pd.DataFrame({'key': regionMask['keys'], 'label': regionMask['labels']})
    
    #call function to check if regionKeyList was defined correctly
    checkKeys(regionKeyList, regionTbl)
    
    #keys to remove (all keys that are note listed in regionKeyList)
    keysToRemove = [key for key in list(regionTbl['key']) if key not in regionKeyList]
    
    #filter elements from the ice thickness DataArray where the region is the desired region
    regionalDataset = dataset.copy()
    for var in dataset.data_vars: 
        if var != 'seaice_conc_monthly_cdr':
            regionalVar = regionalDataset[var]
            for key in keysToRemove: 
                regionalVar = regionalVar.where(regionalVar['region_mask'] != key)
            regionalDataset[var] = regionalVar
    
    #find name of labels 
    labels = [regionTbl[regionTbl['key'] == key]['label'].item() for key in regionKeyList]
    
    #add new attributes describing changes made to the dataset
    if len(labels) < len(regionTbl['key']): 
        if set(regionKeyList) == set([10,11,12,13,15]): #convert to sets so unordered lists are compared
            regionalDataset.attrs['regions with data'] = 'Inner Arctic'
        else:    
            regionalDataset.attrs['regions with data'] = ('%s' % ', '.join(map(str, labels)))
        print('Regions selected: ' + regionalDataset.attrs['regions with data'])
    else: 
        regionalDataset.attrs['regions with data'] = 'All'
        print('Regions selected: All \nNo regions will be removed')
    
    return regionalDataset

def plotThicknessCompare(dataset, var='mean_ice_thickness', var_str = 'sea ice thickness', yearStart = 2018, regionsText = None, figPath = None, label=''):
    """ Plots mean monthly ice thickness for two winter seasons. 

    Args: 
        dataset (xr Dataset): dataset generated by Load_IS2 notebook
        yearStart (int): year to start plot at (default to 2018)
        regionsText(str, optional): string describing regions containing data, if user wants to define the region themself (default to None)
        figPath (str, optional): path to save fig (default to None)
        
    Returns: 
        Figure displayed in notebook 
     
    Restrictions: 
        dataset input needs to contain the following coordinates: mean_ice_thickness_unc, mean_ice_thickness, mean_MYI_thickness, mean_FYI_thickness 
    """
        
    #initialize figure
    fig = plt.figure(figsize=(4, 2.3))
    ax=plt.gca()
    #ax = plt.axes([0, 0, 1, 1]) 
    #title = plt.title('Arctic sea ice thickness', y = 1.11, x = 0.5, fontsize = 'x-large', horizontalalignment = 'center')
    #gridlines = plt.grid(b = True, linestyle = '--', alpha = 0.4) #add gridlines 

    #add title describing regions with data 
    #if 'regions with data' in list(dataset.attrs.keys()): 
        #regionsText = dataset.attrs['regions with data']
        #regionsTitle = ax.text(x = 0.5, y = 1.05, s = 'Region/s: ' + regionsText, size = 12, transform=ax.transAxes, fontsize = 'large', horizontalalignment = 'center')

    #get list of months for plotting x axis 
    winterMonths1 = [3,4,5,6,7,8]
    winterMonthsStr1 = ['Nov','Dec','Jan','Feb','Mar','Apr']
    winterMonths2 = [3,4,5,6,7,8]
    winterMonthsStr2 = [ 'Nov','Dec','Jan','Feb','Mar','Apr']

    #plot data for winter 1
    datasetWinter1 = dataset.sel(time = slice('Nov ' + str(yearStart), 'Apr' + str(yearStart + 1)))
    winter1Str =  str(yearStart)[0:4] + '/' + str(yearStart + 1)[0:4]
    #ax.plot(winterMonths, datasetWinter1['mean_ice_thickness'].values, color = color1, linestyle = '-', marker = 'o', label = '2018-2019')
    #ax.fill_between(winterMonths, datasetWinter1['mean_ice_thickness'].values - datasetWinter1['mean_ice_thickness_unc'].values, datasetWinter1['mean_ice_thickness'].values + datasetWinter1['mean_ice_thickness_unc'].values, facecolor = color1, alpha = 0.1, edgecolor = 'none')
    #ax.errorbar(winterMonths, datasetWinter1['mean_ice_thickness'].values, yerr = datasetWinter1['mean_ice_thickness_unc'].values, color = color1, alpha = 0.7, capsize = 3, zorder = 2, label = winter1Str)        
    #ax.errorbar(winterMonths, datasetWinter1['mean_ice_thickness'].values, yerr=datasetWinter1['mean_ice_thickness_unc'].values, color = 'c', alpha = 0.7, capsize = 3, marker='o', zorder = 2, label = winter1Str)

    ax.plot(winterMonths1, datasetWinter1[var].values, color = 'm', linestyle = '-', marker = 'o', markersize=4, alpha=0.9, label = winter1Str)
    if var=='mean_ice_thickness':
        ax.fill_between(winterMonths1, datasetWinter1[var].values - datasetWinter1['mean_ice_thickness_unc'].values, datasetWinter1['mean_ice_thickness'].values + datasetWinter1['mean_ice_thickness_unc'].values, facecolor = 'm', alpha = 0.07, edgecolor = 'none')

    #plot data for winter 2
    datasetWinter2 = dataset.sel(time = slice('Sep ' + str(yearStart + 1), 'Apr' + str(yearStart + 2)))
    winter2Str = str(yearStart + 1)[0:4] + '/' + str(yearStart + 2)[0:4]

    ax.plot(winterMonths2, datasetWinter2[var].values, color = 'c', linestyle = '-', marker = 'o', markersize=4, label = winter2Str)
    if var=='mean_ice_thickness':
        ax.fill_between(winterMonths2, datasetWinter2[var].values - datasetWinter2['mean_ice_thickness_unc'].values, datasetWinter2[var].values + datasetWinter2['mean_ice_thickness_unc'].values, facecolor = 'c', alpha = 0.09, edgecolor = 'none')

    #plot data for winter 3
    datasetWinter3 = dataset.sel(time = slice('Sep ' + str(yearStart + 2), 'Apr' + str(yearStart + 3)))
    winter3Str = str(yearStart + 2)[0:4] + '/' + str(yearStart + 3)[0:4]
    color3 = 'k'
    ax.plot(winterMonths2, datasetWinter3[var].values, color = 'k', linestyle = '-', marker = 'o', markersize=4, label = winter3Str)
    if var=='mean_ice_thickness':
        ax.fill_between(winterMonths2, datasetWinter2[var].values - datasetWinter2['mean_ice_thickness_unc'].values, datasetWinter2[var].values + datasetWinter3['mean_ice_thickness_unc'].values, facecolor = 'k', alpha = 0.09, edgecolor = 'none')


    #ax.errorbar(winterMonths, datasetWinter2['mean_ice_thickness'].values, yerr=datasetWinter2['mean_ice_thickness_unc'].values, color = 'b', alpha = 0.5, capsize = 3, marker='o', zorder = 2, label = winter2Str)
    #ax.errorbar(winterMonths, datasetWinter2['mean_ice_thickness'].values, yerr = datasetWinter2['mean_ice_thickness_unc'].values, color = color2, alpha = 0.7, capsize = 3, zorder = 2, label = winter2Str)        
    

    ax.set_xticks(winterMonths2)
    ax.set_xticklabels(winterMonthsStr2)

    #ax.annotate('ICESat-2/NESOSIM Winter Arctic sea ice thickness', xy=(0.5, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='center',color='k', zorder=11)

    if var=='mean_freeboard':
        ax.set_ylim([0.1, 0.5])
        ax.set_yticks(np.arange(0.1, 0.51, 0.1))
        ax.set_ylabel('Freebaord (m)')

    if var=='mean_snow_depth':
        ax.set_ylim([0.1, 0.3])
        ax.set_yticks(np.arange(0.1, 0.31, 0.05))
        ax.set_ylabel('Snow depth (m)')

    if var=='mean_ice_thickness':
        #add legend & labels
        ax.legend(loc =2, fontsize = 9, frameon=False)
        ax.set_ylabel('Sea ice thickness (m)')

        ax.set_ylim([0.5, 2.5])
        ax.set_yticks(np.arange(0.5, 2.6, 0.5))

    
    #ax.set_xlabel('Month')
    	
    plt.subplots_adjust(bottom=0.1, left=0.13, top = 0.98, right=0.98)

    plt.savefig(figPath+'seasonal_growth_is2'+label+var+regionsText+'.pdf', dpi = 300)


#----IS2

mapProj = Basemap(epsg=3411,resolution='l', llcrnrlon=279.26, llcrnrlat=48., urcrnrlon=100, urcrnrlat=55.37)

relStr='rel004'
runStr='run4'

is2DataPath='/sea_ice_pso/aapetty/thickness_data/'
cs2DataPath='/sea_ice_pso/aapetty/raw_data/CS2/CS2SMOS_INNERARCTIC/'
figPath='../../Figures/'
savePathIS2=is2DataPath+'/'+relStr+'/'+runStr+'/final_data_gridded/002/'


startYear = 2018
endYear = 2021
winters = getWinterDateRangev2(startYear, endYear) #get date range for winter 18-19 and winter 19-20

is2 = getIS2Data(savePathIS2, winters)

#drop projection variable 
is2 = is2.drop('projection')

#get lat and lon
is2Lats = is2.latitude.isel(time = 0).values
is2Lons = is2.longitude.isel(time = 0).values
is2LonsAttrs = is2.longitude.attrs 
is2LatsAttrs = is2.latitude.attrs

#assign lat and lon as coordinates to dataset
is2 = is2.assign_coords(coords = {'latitude': (('y','x'), is2Lats), 'longitude': (('y','x'), is2Lons)})

regionMask, maskLons, maskLats = getRegionMask('../../AncData/')

#coords and attributes for Region Mask
regionMaskCoords = {'region_mask': (('y','x'), regionMask)}
regionMaskKeys = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 21])
regionMaskLabels = np.array(['non-region oceans', 'Sea of Okhotsk and Japan','Bering Sea','Hudson Bay','Gulf of St. Lawrence',
                    'Baffin Bay, Davis Strait & Labrador Sea','Greenland Sea', 'Barents Seas','Kara Sea','Laptev Sea','East Siberian Sea',
                    'Chukchi Sea','Beaufort Sea','Canadian Archipelago','Arctic Ocean','Land','Coast'])
regionMaskAttrs = {'description': 'NSIDC region mask for the Arctic', 'keys': regionMaskKeys, 'labels' : regionMaskLabels, 'note': 'keys and labels ordered to match by index'}

#add region mask as coordinate to dataset
is2 = is2.assign_coords(coords = regionMaskCoords)

#add descriptive attributes 
is2.region_mask.attrs = regionMaskAttrs
is2.longitude.attrs = is2LonsAttrs
is2.latitude.attrs = is2LatsAttrs
#define a list of keys corresponding to the region of interest
regionKeyList = [10, 11, 12, 13, 15] #Inner Arctic

#restrict data to that region
is2 = restrictRegionally(is2, regionKeyList)


int=True
if int==True:
    var='ice_thickness_int'
    is2means = calcMeans(is2, extraStr='_int')
else:
    var='ice_thickness'
    is2means = calcMeans(is2, extraStr='')


#plotThicknessCompare(is2means, var='mean_'+variable, var_str = variable, yearStart = 2018, regionsText = 'InnerArctic_incKaBa', figPath=figPath, label=relStr+'_'+runStr+str(endYear))


cs2 = pd.read_csv(cs2DataPath+'cs2smos-timeseries-v203-monthly-is2domain.csv')


regionsText = 'Inner_Arctic'
label=relStr+'_'+runStr+str(endYear)


#initialize figure
fig = plt.figure(figsize=(5, 3))
ax=plt.gca()
#ax = plt.axes([0, 0, 1, 1]) 
#title = plt.title('Arctic sea ice thickness', y = 1.11, x = 0.5, fontsize = 'x-large', horizontalalignment = 'center')
#gridlines = plt.grid(b = True, linestyle = '--', alpha = 0.4) #add gridlines 

#add title describing regions with data 
#if 'regions with data' in list(dataset.attrs.keys()): 
    #regionsText = dataset.attrs['regions with data']
    #regionsTitle = ax.text(x = 0.5, y = 1.05, s = 'Region/s: ' + regionsText, size = 12, transform=ax.transAxes, fontsize = 'large', horizontalalignment = 'center')

#get list of months for plotting x axis 



colors=['#a6cee3', '#b2df8a', 'k']
alphas=[0.75, 0.75, 0.8]
num_winters=3
linewidths=[0.8, 0.8, 1.2]
# IS-2
for x in range(num_winters):
    if x==0:
        datasetWinter = is2means.sel(time = slice('Nov ' + str(startYear+x), 'Apr' + str(startYear + x + 1)))
        winterMonths = [3,4,5,6,7,8]
        winterMonthsStr = ['Nov','Dec','Jan','Feb','Mar','Apr']
    else:
        datasetWinter = is2means.sel(time = slice('Oct ' + str(startYear+x), 'Apr' + str(startYear + x + 1)))
        winterMonths = [2, 3,4,5,6,7,8]
        winterMonthsStr = ['Oct','Nov','Dec','Jan','Feb','Mar','Apr']
    
    winterStr =  str(startYear+x)[0:4] + '/' + str(startYear + x + 1)[0:4]
    ax.plot(winterMonths, datasetWinter.values, color = colors[x], linestyle = '-', linewidth=linewidths[x], marker = 'o', markersize=6, alpha=alphas[x], label = 'ICESat-2 ('+winterStr+')')
    
    winterMonths = [2, 3,4,5,6,7,8]
    winterStr =  str(startYear+x)[0:4] + '/' + str(startYear + x + 1)[0:4]
    ax.plot(winterMonths, cs2['mean_sea_ice_thickness'].values[x*7:(x+1)*7], color = colors[x], linestyle = '-', linewidth=linewidths[x], marker = 'v', markersize=6, alpha=alphas[x], label = 'CryoSat-2/SMOS')


#if x==2:
    #    ax.fill_between(winterMonths, datasetWinter[var].values - datasetWinter['mean_ice_thickness_unc'].values, datasetWinter['mean_ice_thickness'].values + datasetWinter['mean_ice_thickness_unc'].values, facecolor = colors[x], alpha = alpha[x], edgecolor = 'none')

# CS-2
#for x in range(num_winters):
    

#ax.errorbar(winterMonths, datasetWinter2['mean_ice_thickness'].values, yerr=datasetWinter2['mean_ice_thickness_unc'].values, color = 'b', alpha = 0.5, capsize = 3, marker='o', zorder = 2, label = winter2Str)
#ax.errorbar(winterMonths, datasetWinter2['mean_ice_thickness'].values, yerr = datasetWinter2['mean_ice_thickness_unc'].values, color = color2, alpha = 0.7, capsize = 3, zorder = 2, label = winter2Str)        


ax.set_xticks(winterMonths)
ax.set_xticklabels(winterMonthsStr)

#add legend & labels
ax.legend(loc =2, fontsize = 8, frameon=False)
ax.set_ylabel('Sea ice thickness (m)')

ax.set_ylim([0.6, 2.4 ])
ax.set_yticks(np.arange(0.5, 2.6, 0.5))


#ax.set_xlabel('Month')
    
plt.subplots_adjust(bottom=0.1, left=0.1, top = 0.98, right=0.98)

plt.savefig(figPath+'seasonal_growth_is2'+label+var+regionsText+'cs2.pdf', dpi = 300)


