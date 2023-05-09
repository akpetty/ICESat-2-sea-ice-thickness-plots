
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

import time
import sys
sys.path.append('../../ICESat-2-sea-ice-thickness-v3/Code/')
import common_functions as cF

relStr='rel005'
runStr='run1'

figPath='/cooler/sea_ice_pso/aapetty/figures/IS2/'+relStr+'/'+runStr+'/'
if not os.path.exists(figPath):
    os.makedirs(figPath)
#figPath='./../../Figures/'
dataPath='/cooler/sea_ice_pso/aapetty/thickness_data/'+relStr+'/'+runStr+'/final_data_along_track/'


date='20201211'
beam='bnum3'
granule=0
lat1=70
lat2=90
atrack1=1
atrack2=3


#IS2data = cF.getATL10Shotdata(dataOutPathM, runStr, campaignStr, cols='all', yearStr=yearStr, monStr=monStr, dayStr=dayStr, fNum=fNum, beamStr=beams[beamNum])
IS2data = cF.get_IS2_along_track(dataPath, date, beam=beam, relStr=relStr[-3:], granule=granule)
print(IS2data['latitude'].values[0], IS2data['latitude'].values[-1])

# Crop using lat/lon first
IS2data=IS2data.where(((IS2data['latitude']>lat1)&(IS2data['latitude']<lat2)), drop=True)

# Reset along_track_distance to be from start of granule
IS2data['along_track_distance']=IS2data['along_track_distance']-IS2data['along_track_distance'].values[0]
IS2data=IS2data.where(((IS2data['along_track_distance']>atrack1)&(IS2data['along_track_distance']<atrack2)), drop=True)


ice_type=IS2data['ice_type'].values[0]
if (ice_type>0.5):
	iceTypeStr='MYI'
else:
	iceTypeStr='FYI'

region_label=cF.get_region_mask_sect_labels(int(IS2data['region_flag'].values[0]))
lon1=str(np.round(IS2data['longitude'].values[0], 2))
lat1=str(np.round(IS2data['latitude'].values[0], 2))

label_str=date+beam+'_gr'+str(granule)+region_label+lon1+lat1

density=str(int(floor(IS2data['snow_density'].values[0])))
#densityWarren=str(int(floor(IS2data['snow_density_W99'].values[0])))

titleStr=region_label+' ('+iceTypeStr+') '+lat1+'N, '+lon1+'E,  Snow density (NESOSIM): '+density+r' kg m$^{-3}$'


fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(10, 8))

sca(axs.flatten()[0])
ax1=gca()
ax1.plot(IS2data['along_track_distance'].values, IS2data['freeboard'].values, 'x-', color='k', alpha=0.8, label='ATL10 '+beam, markersize=5, linewidth=1)
#ax1.errorbar(IS2data['along_track_distance'].values, IS2data['freeboard'].values, yerr=IS2data['freeboard_sigma'].values, fmt='', linestyle='', marker='.', color='k', lw=0.5, capsize=0.5,)

ax1.set_ylabel('freeboard (m)')
#ax1.set_yticks(np.arange(0, 1.5, 0.5))
ax1.set_xticklabels([])

sca(axs.flatten()[1]) 
ax2=gca()

#ax2.plot(IS2data['along_track_distance'].values, IS2data['snow_depth_N'].values, '--', color='m', alpha=0.9,label='NSIM (100 km)', markersize=3, linewidth=1.5)
ax2.plot(IS2data['along_track_distance'].values, IS2data['snow_depth'].values, 'x-', color='k', alpha=0.9, markersize=5, linewidth=1)
#ax2.plot(IS2data['along_track_distance'].values, IS2data['snow_depth_W99mod5'].values, '.-', color='m', alpha=0.9,label='W99m5', markersize=3, linewidth=1.5)

#ax2.yaxis.set_label_position("right")
#ax2.yaxis.tick_right()
#ax2.set_yticks(np.arange(0, 0.4, 0.1))
#ax2.set_ylim([0, 0.35])
ax2.set_ylabel('Snow depth (m)')
ax2.set_xticklabels([])


sca(axs.flatten()[2]) 
ax3=gca()


ax3.plot(IS2data['along_track_distance'].values, IS2data['freeboard'].values-IS2data['snow_depth'].values, 'x-', color='k', markersize=5, linewidth=1)

#ax3.plot(dist, IS2data['ice_thickness_W99mod5dist'].values, '.-',color='r', alpha=0.9, label=r'W99m5$_{rd-pw}$', markersize=3, linewidth=1)
#ax3.plot(dist, IS2data['ice_thickness_NPdist'].values, '.-',  color='b', alpha=0.9, label=r'NSIM$_{rd-pw}$', markersize=3, linewidth=1)
#ax4.plot(IS2data['delta_time'].values-IS2data['delta_time'].iloc[0], IS2data['ice_thickness_W99mod7'].values, '-',color='b', alpha=0.5, label='m7W99', linewidth=1.5)

#ax4.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_N'].values-IS2data['ice_thickness_N_unc'].values,IS2data['ice_thickness_N'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
#ax4.set_yticks(np.arange(0, 5, 1))
ax3.set_ylabel('Eff ice freeboard (m)')
ax3.set_xticklabels([])


sca(axs.flatten()[3]) 
ax4=gca()

ax4.fill_between(IS2data['along_track_distance'].values, IS2data['ice_thickness'].values-IS2data['ice_thickness_unc'].values, IS2data['ice_thickness'].values + IS2data['ice_thickness_unc'].values, alpha=0.5, label='total uncertainity', edgecolor='none', facecolor='k', zorder=1)

ax4.plot(IS2data['along_track_distance'].values, IS2data['ice_thickness'].values, 'x', color='k', markersize=3, linewidth=1)

#ax3.plot(dist, IS2data['ice_thickness_W99mod5dist'].values, '.-',color='r', alpha=0.9, label=r'W99m5$_{rd-pw}$', markersize=3, linewidth=1)
#ax3.plot(dist, IS2data['ice_thickness_NPdist'].values, '.-',  color='b', alpha=0.9, label=r'NSIM$_{rd-pw}$', markersize=3, linewidth=1)
#ax4.plot(IS2data['delta_time'].values-IS2data['delta_time'].iloc[0], IS2data['ice_thickness_W99mod7'].values, '-',color='b', alpha=0.5, label='m7W99', linewidth=1.5)

#ax4.fill_between(IS2data['delta_time'].values-IS2data['delta_time'].values[0], IS2data['ice_thickness_N'].values-IS2data['ice_thickness_N_unc'].values,IS2data['ice_thickness_N'].values+IS2data['ice_thickness_NPdist_unc'].values, alpha=0.3, edgecolor='none', facecolor='b', label='NESOSIM_dist')
#ax4.set_yticks(np.arange(0, 5, 1))
ax4.set_ylabel('Ice thickness (m)')
ax4.set_xticklabels([])


sca(axs.flatten()[4]) 
ax5=gca()
#ax4.yaxis.set_label_position("right")
#ax4.yaxis.tick_right()
#ax.fill_between(days[x], ma.mean(snowBudgetRegionAll[x][n], axis=0)-ma.std(snowBudgetRegionAll[x][n], axis=0), ma.mean(snowBudgetRegionAll[x][n], axis=0)+ma.std(snowBudgetRegionAll[x][n], axis=0), alpha=0.3, edgecolor='none', facecolor=colors[x], zorder=1)
ax5.plot(IS2data['along_track_distance'].values, IS2data['ice_thickness_unc'].values, '.-', alpha=0.9, label='total', color='k', markersize=3, linewidth=1)
ax5.plot(IS2data['along_track_distance'].values, IS2data['ice_thickness_uncrandom'].values, '.-', alpha=0.9, label='random', color='c', markersize=3, linewidth=1)
ax5.plot(IS2data['along_track_distance'].values, IS2data['ice_thickness_uncsys'].values, '.-', alpha=0.9, label='systematic', color='m', markersize=3, linewidth=1)

ax5.set_ylabel('Uncertainity (m)')
ax5.set_xlabel('Along track distance (km)')

ax1.legend(loc=1, ncol=1, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
ax4.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
#ax3.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)
ax5.legend(loc=1, ncol=6, markerscale=3, bbox_to_anchor=[0.98, 1.15], labelspacing=0.45, frameon=False)

# Common subplot additions
for ax in axs.flatten():
    sca(ax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    #ax.set_xlim([0, 1000])
    ax.yaxis.set_ticks_position('left')
    ax.grid(axis='both', linestyle='--')
    ax.yaxis.set_major_formatter(FormatStrFormatter("%1.2f"))

ax1.annotate(titleStr, xy=(0.01, 1.01), xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom')

subplots_adjust(left = 0.07, right = 0.98, bottom=0.07, top = 0.96, hspace=0.18)
plt.savefig(figPath+'/IS2SITDAT4'+label_str+'.png', dpi=500)
#plt.savefig(figPathM+'/ts4'+campaignStr+'_F'+str(fileNum)+'shotData.pdf')
#plt.savefig(figPathM+'/ts3'+labelStr+runStr+'_F'+str(fNum)+'shotData.png', dpi=500)
#fig.show()

