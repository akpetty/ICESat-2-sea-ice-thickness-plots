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


import time
import sys
sys.path.append('../')
import common_functions as cF
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
cF.reset_matplotlib()


relStr='rel004'
runStr='run4'


#figPath='../../Figures/'
figPath='./figures/'
if not os.path.exists(figPath):
    os.makedirs(figPath)

dataPath='/sea_ice_pso/aapetty/thickness_data/'+relStr+'/'+runStr+'/raw/'
concDataPath='/sea_ice_pso/aapetty/raw_data/ICECONC/CDR/monthly/'
iceTypePath='/sea_ice_pso/aapetty/raw_data/ICETYPE/OSISAF/'

beams=['bnum1', 'bnum3', 'bnum5']
dayStr='*'
month=4
monStr='%02d' %(month)

fNum=-1 #-1=all
yearStr='2019'


cols=['snow_depth_NPdist','snow_density_N','ice_thickness_NPdist']
cols.append('lat')
cols.append('lon')
cols.append('ssh_flag')
cols.append('seg_length')

IS2data = cF.getProcessedATL10ShotdataNCDF(dataPath, vars=cols, yearStr=yearStr, ssh_mask=0, monStr=monStr, dayStr=dayStr, fNum=fNum, minseg=4, maxseg=200, concat=1, beamStrs=beams)


m = Basemap(projection='npstere',boundinglat=60,lon_0=0, resolution='l' , round=False)

xpts, ypts=m(IS2data['lon'].values, IS2data['lat'].values)
xptsc, yptsc, iceconc=cF.get_cdr_conc(concDataPath, m, yearStr, monStr)

xptst, yptst, icetype=cF.getIceTypeRaw(iceTypePath, m, '15', monStr, yearStr, res=1)

xpts1, ypts1=m(176.16, 85.7)
xpts2, ypts2=m(87.34, 78.5)

vmins=[0, 200, 0]
vmaxs=[0.5, 350, 5]
cmaps=[cm.YlOrRd, cm.BuPu, cm.cubehelix_r]

unit='m'
labelStr=runStr+'-'+'-'+yearStr+monStr+dayStr+'_'+str(size(beams))+'bms'
monLabels=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
monLabel=monLabels[month-1]

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(7, 3))

#============ Freeboard ============
ax1=axs.flatten()[0]
sca(ax1)
im1=hexbin(xpts, ypts, C=IS2data[cols[0]].values, gridsize=5000,
        cmap=cmaps[0], vmin=vmins[0], vmax=vmaxs[0], zorder=2, rasterized=True)
m.drawcoastlines(linewidth=0.25, zorder=5)
m.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
m.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
m.fillcontinents(color='0.9',lake_color='grey', zorder=3)
im11 = m.contour(xptsc , yptsc, iceconc,levels=[0.5], colors='m', linewidths=0.5, zorder=5, alpha=1)
im21 = m.contour(xptst , yptst, icetype,levels=[0.5], colors='k', linewidths=0.4, zorder=5, alpha=0.7)

if (month==11):
	im12=plot(xpts1, ypts1, '*', color='k', markersize=5, zorder=5)
	im13=plot(xpts2, ypts2, '^', color='k', markersize=4, zorder=5)

cax1 = fig.add_axes([0.02, 0.14, 0.3, 0.035])
cbar=fig.colorbar(im1, cax=cax1, orientation='horizontal',extend='both')
cbar.set_label('snow depth (m)', labelpad=3)
cbar.set_ticks(np.arange(0, vmaxs[0]+0.1, 0.2))

ax1.annotate(' '+monLabel+' '+yearStr, xy=(0.02, 0.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

#============ Snow depth ============
sca(axs.flatten()[1])  

im2=hexbin(xpts, ypts, C=IS2data[cols[1]].values, gridsize=5000,
        cmap=cmaps[1], vmin=vmins[1], vmax=vmaxs[1], zorder=2, rasterized=True)
m.drawcoastlines(linewidth=0.25, zorder=5)
m.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
m.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
m.fillcontinents(color='0.9',lake_color='grey', zorder=3)
im11 = m.contour(xptsc , yptsc, iceconc,levels=[0.5], colors='m', linewidths=0.5, zorder=5, alpha=1)
im21 = m.contour(xptst , yptst, icetype,levels=[0.5], colors='k', linewidths=0.4, zorder=5, alpha=0.7)

if (month==11):
	im12=plot(xpts1, ypts1, '*', color='k', markersize=5, zorder=5)
	im13=plot(xpts2, ypts2, '^', color='k', markersize=4, zorder=5)

cax2 = fig.add_axes([0.35, 0.14, 0.3, 0.035])
cbar2 = colorbar(im2,cax=cax2, orientation='horizontal', extend='both', use_gridspec=True)
cbar2.set_label('NESOSIM snow density (m)', labelpad=3)
cbar2.set_ticks(np.arange(0, vmaxs[1]+0.1, 0.1))

#============ THICKNESS ============
sca(axs.flatten()[2])  
im3=hexbin(xpts, ypts, C=IS2data[cols[2]].values, gridsize=5000,
        cmap=cmaps[2], vmin=vmins[2], vmax=vmaxs[2], zorder=2, rasterized=True)
m.drawcoastlines(linewidth=0.25, zorder=5)
m.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
m.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
m.fillcontinents(color='0.9',lake_color='grey', zorder=3)
im11 = m.contour(xptsc , yptsc, iceconc,levels=[0.5], colors='m', linewidths=0.5, zorder=5, alpha=1)
im21 = m.contour(xptst , yptst, icetype,levels=[0.5], colors='k', linewidths=0.4, zorder=5, alpha=0.7)

if (month==11):
	im12=plot(xpts1, ypts1, '*', color='k', markersize=5, zorder=5)
	im13=plot(xpts2, ypts2, '^', color='k', markersize=4, zorder=5)

cax3 = fig.add_axes([0.68, 0.14, 0.3, 0.035])
cbar3 = colorbar(im3,cax=cax3, orientation='horizontal', extend='both', use_gridspec=True)
cbar3.set_label('Sea ice thickness (m)', labelpad=3)
cbar3.set_ticks(np.arange(0, vmaxs[2]+0.1, 1))

subplots_adjust(bottom=0.18, left=0.01, top = 0.98, right=0.99, wspace=0.02, hspace=0.2)

fig.savefig(figPath+'/map3_'+cols[-4]+labelStr+runStr+'.png', dpi=300)








