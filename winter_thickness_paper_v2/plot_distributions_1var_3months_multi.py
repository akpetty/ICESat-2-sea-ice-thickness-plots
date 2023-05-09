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

def get_hists(dataOutPathT, start_date, end_date, varT, binVals, binWidth, maxValue, min_seg_len=4, max_seg_len=200, factor=1):
        
        vars=[varT, 'seg_length', 'region_flag', 'ssh_flag']
        #if (monStr=='12'):
        #    dataOutPath=dataPathIS2+releaseStr+'/'+runStr+'/raw/'
        #else:
        #    dataOutPath=dataPathIS2+releaseStr+'/'+runStr+'/raw/'
        print(dataOutPathT)
        dFbeams = cF.getProcessedATL10ShotdataNCDF_daterange(dataOutPathT, start_date, end_date, ssh_mask=0, minseg=min_seg_len, maxseg=max_seg_len, concat=1, vars=vars, beamStr=beam)
        print('Got data')
        dFbeams=dFbeams.where(dFbeams[varT]>0.0, drop=True)
        dFbeams=dFbeams.where(dFbeams[varT]<maxValue*5, drop=True)
        dFbeams=dFbeams.where(~np.isnan(dFbeams[varT]), drop=True)
        dFbeams=dFbeams.where(dFbeams.seg_length>min_seg_len, drop=True)

        dFbeams=dFbeams.where(dFbeams.seg_length<max_seg_len, drop=True)
        
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


#labelStrs=['r002', 'r003', 'r004']


beam='bnum1'
snowVar='NPdist'
#var='freeboard'
#variables=['freeboard', 'snow_depth_'+snowVar, 'ice_thickness_'+snowVar]
#varStrs=['freeboard', 'snow depth', 'sea ice thickness']

variable='freeboard'

cmap = cm.viridis
#c_colors = cmap(np.linspace(0.1, 0.9, size(releaseStrs)))
#color=plt.get_cmap("viridis", 10)(range(0,8))

#c_colors = ['#377eb8', '#ff7f00', '#4daf4a',
#                  '#984ea3', '#a65628','#f781bf',
#                  '#999999', '#e41a1c', '#dede00']

c_colors = cm.viridis(np.linspace(0, 0.9, 4))


if variable =='freeboard':
    maxValue=80
    factor = 100.
    formatter = '%.1f'
    varStr='total freeboard'
    #c_colors = ['tab:orange','tab:blue', 'tab:green', 'tab:gray']
    unit='cm'
    

    releaseStrs=['rel002', 'rel003', 'rel004', 'rel005']
    runStrs=['run1', 'run1', 'run1', 'run1']
    labelStrs=['rel002', 'rel003', 'rel004', 'rel005']

    #cmap = cm.viridis
    #cmap=plt.cm.get_cmap('Set2', 8)
    #c_colors = cmap(np.linspace(0, 0.9, 5))
    
elif variable =='snow_depth_NPdist':
    maxValue=60
    factor = 100.
    unit='cm'
    formatter = '%.1f'
    varStr='snow depth'
    #c_colors = ['k', '#1f78b4', '#33a02c', '#b2df8a']

    releaseStrs=['rel002','rel004', 'rel004', 'rel005']
    runStrs=['run1','run1_Nv1', 'run1', 'run1']
    labelStrs=['rel002/Nv1.0','rel004/Nv1.0', 'rel004/Nv1.1', 'rel005/Nv1.1']

    #cmap = cm.viridis
    #c_colors = cmap(np.linspace(0.1, 0.9, size(releaseStrs)))

elif variable =='ice_thickness_NPdist':
    maxValue=5
    factor=1.
    unit='m'
    formatter = '%.2f'
    varStr='sea ice thickness'

    releaseStrs=['rel002','rel004', 'rel004', 'rel005']
    runStrs=['run1','run1_Nv1', 'run1', 'run1']
    labelStrs=['rel002/Nv1.0','rel004/Nv1.0', 'rel004/Nv1.1', 'rel005/Nv1.1']

    #cmap = cm.viridis
    #c_colors = cmap(np.linspace(0.1, 0.9, size(releaseStrs)))

    #c_colors = ['k', '#1f78b4', '#33a02c', '#b2df8a']

    #releaseStrs=['rel002', 'rel004', 'rel004', 'rel004']
    #runStrs=['run12', 'run4_Nv1', 'run4', 'run4_clim']
    #labelStrs=['r002/Nv1.0', 'r004/Nv1.0', 'r004/Nv1.1', 'r004/Nv1.1c']

    #c_colors = ['k', '#a6cee3', '#1f78b4']
    #releaseStrs=['rel002', 'rel003', 'rel004']
    #runStrs=['run1', 'run1', 'run1_Nv1']
    #labelStrs=['r002', 'r003', 'r004']

    #c_colors = ['#1f78b4', '#33a02c', '#b2df8a']
    #releaseStrs=['rel004', 'rel004', 'rel004']
    #runStrs=['run4_Nv1', 'run4', 'run4_clim']
    #labelStrs=['r004/Nv1.0', 'r004/Nv1.1', 'r004/Nv1.1c']


start_dates=['20181101', '20190101', '20190401']
end_dates=['20181130', '20190131', '20190430']

region_label='Arctic Ocean (inc Kara)'


numBins=40
histVals=ma.masked_all((np.size(start_dates), size(releaseStrs), numBins))
counts=ma.masked_all((np.size(start_dates), size(releaseStrs)))
means=ma.masked_all((np.size(start_dates), size(releaseStrs)))
binVals=[]
binWidths=[]


for x in range(np.size(start_dates)):
    #maxValue=np.percentile((dFbeams[var].values), 95)
    
    binValsT=np.linspace(0, maxValue, numBins+1)

    binVals.append(binValsT)
    binWidthT=binValsT[1]-binValsT[0]
    binWidths.append(binWidthT)

    for m in range(size(releaseStrs)):

        #baseDataPath='/cooler/sea_ice_pso/aapetty/thickness_data/'
        #regions=[9, 10, 11, 12, 13, 15]

        if (releaseStrs[m]=='rel005'):
            baseDataPath='/cooler/sea_ice_pso/aapetty/thickness_data/'
            regions=[1,2,3,4,5,6]

        elif (releaseStrs[m]=='rel004'):
            baseDataPath='/cooler/sea_ice_pso/aapetty/thickness_data/'
            regions=[9, 10, 11, 12, 13, 15]

        else:
            baseDataPath='/cooler/freezer/thickness_data/'
            regions=[9, 10, 11, 12, 13, 15]

        dataOutPath=baseDataPath+releaseStrs[m]+'/'+runStrs[m]+'/raw/'

        histVals[x, m], means[x, m], counts[x, m] = get_hists(dataOutPath, start_dates[x], end_dates[x], variable, binValsT, binWidthT, maxValue, factor=factor)

    savetxt(figPath+'/stats/innerarctic_4vars_releasecomp_means'+start_dates[x]+'-'+end_dates[x]+variable+'.txt', means, fmt='%1.3f')
    savetxt(figPath+'/stats/innerarctic_4vars_releasecomp_counts'+start_dates[x]+'-'+end_dates[x]+variable+'.txt', counts, fmt='%.0f')



fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(3.5, 5))

for x in range(np.size(start_dates)):
    ax=axs[x]
    sca(ax)
    #ax.annotate('Means:', xy=(0.98, 0.9), color='k', xycoords='axes fraction', horizontalalignment='right', verticalalignment='bottom')

    for m in range(np.size(releaseStrs)):
        im = plt.plot(binVals[x][0:-1]+binWidths[x]/2., histVals[x, m],linestyle='-', linewidth=1.2, label=labelStrs[m], color=c_colors[m], alpha=0.8)
        
        #meanStr=str(np.round(means[x, m], 2))
        ax.axvline(x=means[x, m], linestyle='--', linewidth=1.2, color=c_colors[m], alpha=0.8)
        mean_str=formatter %(means[x, m])

        ax.annotate(labelStrs[m]+': '+mean_str+' '+unit, xy=(0.97, 0.89-(0.1*m)), color=c_colors[m], xycoords='axes fraction', horizontalalignment='right', verticalalignment='bottom')

    ax.annotate('('+chr(97+x)+') '+start_dates[x][0:4]+'-'+start_dates[x][4:6], xy=(0.01, 1.01), color='k', xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom')

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

subplots_adjust(bottom=0.075, left=0.15, right=0.97, top=0.97, hspace=0.15)
plt.savefig(figPath+'/distributions'+variable+'_'+str(size(releaseStrs))+'rels'+start_dates[0]+end_dates[-1]+region_label+'.pdf', dpi=300)



