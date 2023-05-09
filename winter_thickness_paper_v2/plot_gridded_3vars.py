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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colorbar as mcbar

import time
import sys
sys.path.append('../../ICESat-2-sea-ice-thickness-v2/Code/')
import common_functions as cF
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean
cF.reset_matplotlib()

mapProj = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l' , round=False)


relStr='rel004'
runStr='run4'
version_str='002'
figPath='./figures/'
#figPath='/sea_ice_pso/aapetty/figures/IS2/'+relStr+'/'+runStr+'/Maps/'
if not os.path.exists(figPath):
    os.makedirs(figPath)

IS2_path='/sea_ice_pso/aapetty/thickness_data/'+relStr+'/'+runStr+'/final_data_gridded/'+version_str+'/'
concDataPath='/sea_ice_pso/aapetty/raw_data/ICECONC/CDR/monthly/v4/'
iceTypePath='/sea_ice_pso/aapetty/raw_data/ICETYPE/OSISAF/'
region_mask_data='/sea_ice_pso/aapetty/raw_data/OTHER/'

yearStrs = ['2020', '2021']
monStr='04'
monLabel=cF.monLabels(int(monStr)-1)

strs=[monLabel+' '+yearStrs[0], monLabel+' '+yearStrs[1]]
print(strs)

#var='freeboard'
#var='snow_depth'
var='freeboard'

if var =='freeboard':
    factor=100.
    minValue = 0
    maxValue=80
    minValueDiff=-20
    maxValueDiff=20
    tickint= 5
    varStr='sea ice freeboard'
    units='(cm)'
    cmap = plt.cm.YlOrRd
    labels=['(a) '+ strs[0], '(b) '+strs[1], '(c) '+strs[1][4:]+' - '+strs[0][4:]]


elif var =='snow_depth_int':
    factor=100.
    minValue = 0
    maxValue=60
    minValueDiff=-20
    maxValueDiff=20
    tickint= 5
    varStr='snow depth'
    units='(cm)'
    labels=['(d) ', '(e) ', '(f)']
    cmap = plt.cm.BuPu


elif var =='ice_thickness_int':
    factor=1.
    minValue = 0
    maxValue=6
    minValueDiff=-1.5
    maxValueDiff=1.5
    tickint= 7
    varStr='sea ice thickness'
    units='(m)'
    labels=['(g) ', '(h) ', '(i)']
    #labels=['(a) '+ strs[0], '(b) '+strs[1], '(c) '+strs[1][4:]+' - '+strs[0][4:]]
    #cmap=cmocean.cm.thermal
    cmap=plt.cm.viridis


minval=[minValue, minValue, minValueDiff]
maxval=[maxValue, maxValue, maxValueDiff]

#---declare empty arrays
data=np.empty((3,448, 304))
data[:] = np.nan

iceconcs=np.empty((2,448, 304))
iceconcs[:] = np.nan

icetypes=np.empty((2,1120, 760))
icetypes[:] = np.nan

for i in range(len(yearStrs)):
    
    iceconcs[i]=cF.get_cdr_conc_v2(concDataPath, yearStrs[i], monStr)
    xpts_it, ypts_it, icetypes[i]=cF.getIceTypeRaw(iceTypePath, mapProj, '15', monStr, yearStrs[i], res=1)
    xpts, ypts, _, _, dataT = cF.getIS2gridded(IS2_path, '01_'+yearStrs[i]+monStr+'_'+relStr[3:]+'_'+version_str, mapProj, variable=var, poleHole=100)

    data[i]=dataT*factor

data[2] = data[1] - data[0]

#regions = [9, 10, 11, 12, 13, 15] #Inner Arctic inc Kara
#regionMask_mask = np.isin(regionMask, regions).astype(int)

cbarlabels=[varStr+' '+units, varStr+' '+units, varStr+' difference '+units]
cmaps=[cmap,  cmap, plt.cm.RdBu_r]


fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(8, 4.3))
plt.subplots_adjust(bottom=0.01, wspace=0.03, left=0.02, top=0.94, right=0.98)

for i in range(len(data)):
    ax=axs.flatten()[i]
    plt.sca(ax)
    
    im1 = mapProj.contourf(xpts , ypts, data[i] , levels=np.linspace(minval[i], maxval[i]+0.1, tickint*5), cmap=cmaps[i], vmin=minval[i], vmax=maxval[i], extend='both', shading='gouraud', edgecolors='None', zorder=2, rasterized=True)
    
    if i<2:
        im11=mapProj.pcolormesh(xpts, ypts, iceconcs[i], vmin=0, vmax=2, cmap=plt.cm.gray_r,  zorder=1)
    
    #if i<1:
    #    im12 = mapProj.contour(lon2d_greater , lat2d_greater, regionMask_mask_g, colors='m', linewidths=1, levels=[0.5],transform=ccrs.PlateCarree(), zorder=3)

    if i<2:
        im13 = mapProj.contour(xpts_it , ypts_it, icetypes[i], colors='k', linewidth=3, levels=[0.5], zorder=3)

    #ax.coastlines(resolution='50m', linewidth=0.22, zorder=3)
    #ax.gridlines(draw_labels=False,
    #          linewidth=0.22, color='gray', alpha=0.5, linestyle='--',zorder=4)
    #ax.add_feature(cfeature.LAND, facecolor='0.95', zorder=2)

    mapProj.fillcontinents(color='0.95',lake_color='grey', zorder=2)
    mapProj.drawcoastlines(linewidth=0.25, zorder=5)

    #ax.set_extent([-179, 179, 48, 90], ccrs.PlateCarree())
    ax.set_xlim([np.percentile(xpts, 3), np.percentile(xpts, 97)])

    cax,kw = mcbar.make_axes(ax,location='bottom',pad=0.01,shrink=0.8)
    cb=fig.colorbar(im1,cax=cax,extend='both',**kw)

    cb.set_ticks(np.linspace(minval[i], maxval[i], tickint) )
    cb.set_label(cbarlabels[i], labelpad=3)
    #im1=hexbin(xpts_sections, ypts_sections, C=data[i], vmin=minval[i], vmax=maxval[i], gridsize=gridbins, cmap=cmaps[i], zorder=2, rasterized=True)
    
    ax.annotate(labels[i], xy=(0.02, 1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

#plt.tight_layout()
fig.savefig(figPath+'/comp_maps'+var+strs[0][0:3]+'-'+strs[0][0:3]+'bm.png', dpi=300)



