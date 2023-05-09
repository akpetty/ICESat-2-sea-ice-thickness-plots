import matplotlib, sys
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import numpy as np
from pylab import *
import numpy.ma as ma
import xarray as xr
import pandas as pd
import os
from glob import glob
import netCDF4 as nc4
from scipy.interpolate import griddata
from netCDF4 import Dataset
import scipy.stats as st
import matplotlib.ticker as ticker
import time
import sys
sys.path.append('../../ICESat-2-sea-ice-thickness-v2/Code/')
import common_functions as cF
import seaborn as sns
import dask.array as da

cF.reset_matplotlib()

def get_hists(dataOutPathT, start_date, end_date, varT, binVals, binWidth, maxValue, beam, min_seg_len=4, max_seg_len=200, factor=1.):
        
        vars=[varT, 'seg_length', 'region_flag', 'ssh_flag']
        #if (monStr=='12'):
        #    dataOutPath=dataPathIS2+releaseStr+'/'+runStr+'/raw/'
        #else:
        #    dataOutPath=dataPathIS2+releaseStr+'/'+runStr+'/raw/'
        print(dataOutPathT)
        dFbeams = cF.getProcessedATL10ShotdataNCDF_daterange(dataOutPathT, start_date, end_date, ssh_mask=0, minseg=min_seg_len, maxseg=max_seg_len, concat=1, vars=vars, beamStr=beam)
        print('Got data')
        dFbeams=dFbeams.where(dFbeams[varT]>0.0, drop=True)
        dFbeams=dFbeams.where(dFbeams[varT]<30, drop=True)
        dFbeams=dFbeams.where(~np.isnan(dFbeams[varT]), drop=True)
        dFbeams=dFbeams.where(dFbeams.seg_length>4, drop=True)

        dFbeams=dFbeams.where(dFbeams.seg_length<200, drop=True)
        
        vals=dFbeams[varT][np.isin(dFbeams.region_flag, regions)]*factor
        segs=dFbeams['seg_length'][np.isin(dFbeams.region_flag, regions)]

        weights=segs/segs.sum().values

        countsT=vals.count().values
        meansT=(vals*segs).sum().values/segs.sum().values

        h, bins = da.histogram(vals.data, bins=size(binVals)-1, range=[0, maxValue], weights=weights.data)

        histValsT=h.compute()

        return histValsT, meansT, countsT


figPath='./figures/'
if not os.path.exists(figPath):
    os.makedirs(figPath)
baseDataPath='/sea_ice_pso/aapetty/thickness_data/'

#labelStrs=['r002', 'r003', 'r004']

beams=['bnum1', 'bnum3', 'bnum5']
beams_str=['beam 1', 'beam 3', 'beam 5']

releaseStr='rel005'
runStr='run1'

snowVar='NPdist'
#var='freeboard'
#variables=['freeboard', 'snow_depth_'+snowVar, 'ice_thickness_'+snowVar]
#varStrs=['freeboard', 'snow depth', 'sea ice thickness']

variable='freeboard'
unit='cm'

if variable =='freeboard':
    maxValue=80
    factor = 100.
    formatter = '%.1f'
    varStr='total freeboard'
    

elif variable =='snow_depth_'+snowVar:
    maxValue=0.6
    varStr='snow depth'

elif variable =='ice_thickness_'+snowVar:
    maxValue=5
    varStr='sea ice thickness'
    #c_colors = ['k', '#1f78b4', '#33a02c', '#b2df8a']


cmap = cm.magma
c_colors = cmap(np.linspace(0, 0.8, size(beams)))

start_dates=['20181101', '20190101', '20190401']
end_dates=['20181130', '20190131', '20190430']

region_label='Arctic Ocean (inc Kara)'
regions=[9, 10, 11, 12, 13, 15]

if releaseStr=='rel005':
    regions=[1,2,3,4,5,6]

numBins=40
histVals=ma.masked_all((np.size(start_dates), size(beams), numBins))
counts=ma.masked_all((np.size(start_dates), size(beams)))
means=ma.masked_all((np.size(start_dates), size(beams)))
binVals=[]
binWidths=[]
#maxValue=[0.8, 0.6, 5]

for x in range(np.size(start_dates)):
    #maxValue=np.percentile((dFbeams[var].values), 95)
    
    binValsT=np.linspace(0, maxValue, numBins+1)

    binVals.append(binValsT)
    binWidthT=binValsT[1]-binValsT[0]
    binWidths.append(binWidthT)

    for m in range(size(beams)):

        dataOutPath=baseDataPath+releaseStr+'/'+runStr+'/raw/'

        histVals[x, m], means[x, m], counts[x, m] = get_hists(dataOutPath, start_dates[x], end_dates[x], variable, binValsT, binWidthT, maxValue, beams[m], factor=factor)

    savetxt(figPath+'/stats/beamcompcomp_means'+start_dates[x]+'-'+end_dates[x]+variable+'.txt', means, fmt='%1.3f')
    savetxt(figPath+'/stats/beamcomp_counts'+start_dates[x]+'-'+end_dates[x]+variable+'.txt', counts, fmt='%.0f')



fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(3.5, 5))

for x in range(np.size(start_dates)):
    ax=axs[x]
    sca(ax)
    #ax.annotate('Means:', xy=(0.98, 0.9), color='k', xycoords='axes fraction', horizontalalignment='right', verticalalignment='bottom')

    for m in range(np.size(beams)):
        im = plt.plot(binVals[x][0:-1]+binWidths[x]/2., histVals[x, m],linestyle='-', linewidth=1.4, label=beams[m], color=c_colors[m], alpha=0.8)
        
        #meanStr=str(np.round(means[x, m], 2))
        ax.axvline(x=means[x, m], linestyle='--', linewidth=1.4, color=c_colors[m])
        mean_str=formatter %(means[x, m])

        ax.annotate(beams_str[m]+': '+mean_str+' '+unit, xy=(0.97, 0.89-(0.1*m)), color=c_colors[m], xycoords='axes fraction', horizontalalignment='right', verticalalignment='bottom')

    ax.annotate('('+chr(97+x+3)+') '+start_dates[x][0:4]+'-'+start_dates[x][4:6], xy=(0.01, 1.01), color='k', xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom')

    ax.set_xlim(0, maxValue)
    ax.set_ylim(0, 0.13)
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.04))
    #ax.grid(which='major', axis='both', linestyle='--')

    if (x==1):
        ax.set_ylabel('Probability density') 

    if (x==2):
        ax.set_xlabel(varStr+' ('+unit+')', labelpad=3)
    else:
        ax.set_xticklabels([]) 
        ax.set_xticks([])

subplots_adjust(bottom=0.075, left=0.15, right=0.97, top=0.97, hspace=0.15)
plt.savefig(figPath+'/distributions'+variable+'_beams'+str(size(beams))+'rels'+releaseStr+runStr+start_dates[0]+end_dates[-1]+region_label+'.pdf', dpi=300)



