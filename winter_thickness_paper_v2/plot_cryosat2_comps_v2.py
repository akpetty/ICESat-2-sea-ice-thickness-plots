
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
sys.path.append('../')
import common_functions as cF

import pyproj
import cartopy.crs as ccrs
import cartopy
cF.reset_matplotlib()

def getCS2gsfc(dataPathCS2, yearStr, mStr, proj):

	f = Dataset(dataPathCS2+'/GSFC/'+yearStr+'/RDEFT4_'+yearStr+mStr+'15.nc', 'r')
	thicknessCS = f.variables['sea_ice_thickness'][:]
	thicknessCS=ma.masked_where(thicknessCS<0, thicknessCS)
	thicknessCS=ma.masked_where(thicknessCS>15, thicknessCS)
	thicknessCS=ma.masked_where(np.isnan(thicknessCS), thicknessCS)
	#thicknessCS=thicknessCS.filled(0)
	#thicknessCS=ma.masked_where(region_mask!=1, thicknessCS)

	latsCS = f.variables['lat'][:]
	lonsCS = f.variables['lon'][:]

	xptsT, yptsT = proj(lonsCS, latsCS)

	#files = glob(dataPath+ystr+mstr+'*')
	return xptsT, yptsT, thicknessCS

def getCS2cpom(dataPathCS2, yearStr, mStr, proj):
	print(dataPathCS2+'/CPOM/thk_'+yearStr+'_'+mStr+'*')
	files=glob(dataPathCS2+'/CPOM/thk_'+yearStr+'_'+mStr+'*')
	
	print(files[0])
	f = Dataset(files[0], 'r')
	latsCS = f.variables['latitude'][::5]
	lonsCS = f.variables['longitude'][::5]
	
	thicknessCS = f.variables['thickness'][::5]

	xptsT, yptsT = proj(lonsCS, latsCS)

	#files = glob(dataPath+ystr+mstr+'*')
	return xptsT, yptsT, thicknessCS

def getCS2jpl(dataPathCS2, yearStr, mStr, proj):
	flat = open(dataPathCS2+'JPL/psn25lats_v3.dat', 'rb')
	flon = open(dataPathCS2+'JPL/psn25lons_v3.dat', 'rb')
	lats = fromfile(file=flat, dtype='<i4')/100000.0
	latsCS = reshape(lats, [448, 304])
	lons = fromfile(file=flon, dtype='<i4')/100000.0
	lonsCS = reshape(lons, [448, 304])
	
	lonsCS=lonsCS[135:-113, 53:-51]
	latsCS=latsCS[135:-113, 53:-51]

	mon1=loadtxt(dataPathCS2+'JPL/month_'+mStr+'.txt')
	thicknessCS=np.flip(mon1, axis=0) 

	xptsT, yptsT = proj(lonsCS, latsCS)

	#files = glob(dataPath+ystr+mstr+'*')
	return xptsT, yptsT, thicknessCS


def getCS2awi(dataPathCS2, yearStr, mStr, proj):
	print(dataPathCS2+'/AWI/*'+yearStr+mStr+'*')
	files=glob(dataPathCS2+'/AWI/*'+yearStr+mStr+'*')
	
	print(files[0])
	f = Dataset(files[0], 'r')
	latsCS = f.variables['lat'][:]
	lonsCS = f.variables['lon'][:]
	
	thicknessCS = f.variables['sea_ice_thickness'][:]

	xptsT, yptsT = proj(lonsCS, latsCS)

	#files = glob(dataPath+ystr+mstr+'*')
	return xptsT, yptsT, thicknessCS


def main(CS2PRODUCT, date):
	
	dataPathIS2='/sea_ice_pso/aapetty/thickness_data/'
	dataPathCS2='/sea_ice_pso/aapetty/raw_data/CS2/'

	mStr=date[4:]
	monLabel=cF.monLabels(int(mStr)-1)
	yearStr=str(date[0:4])
	print(yearStr, mStr)

	relStr='rel003'
	runStr='run1'
	figPath='../../Figures/'
	beamStr='bnum1'
	version='psn'

	savePath=dataPathIS2+'/'+relStr+'/'+runStr+'/products/'
	labelStr=runStr+'-'+beamStr+'-'+yearStr+monLabel

	#xptsIS2, yptsIS2, latG, lonG, proj = cF.create_grid(epsg_string='3413', dxRes=25000., lllat=30, llon=-90, urlat=30, urlon=90)
	xptsIS2, yptsIS2, latG, lonG, proj = cF.create_grid_nsidc()
	#proj = pyproj.Proj("+init=EPSG:3413")
	#mapProj = Basemap(epsg='3413',lon_0=0.,resolution='l',width=6000000,height=6000000)

	#mapProj = pyproj.Proj("+init=EPSG:3413")

	if (CS2PRODUCT=='AWI'):
		snowVar='awi'
		xptsCS2, yptsCS2,ice_thicknessCS2 = getCS2awi(dataPathCS2, yearStr, mStr, proj)
	elif (CS2PRODUCT=='CPOM'):
		snowVar='cpom'
		xptsCS2, yptsCS2,ice_thicknessCS2 = getCS2cpom(dataPathCS2, yearStr, mStr, proj)
	elif (CS2PRODUCT=='JPL'):
		snowVar='nasa'
		xptsCS2, yptsCS2,ice_thicknessCS2 = getCS2jpl(dataPathCS2, yearStr, mStr, proj)
	elif (CS2PRODUCT=='GSFC'):
		snowVar='nasa'
		xptsCS2, yptsCS2,ice_thicknessCS2 = getCS2gsfc(dataPathCS2, yearStr, mStr, proj)
	else:
		print('no data')

	resolution=25.
	ice_thicknessIS2=np.load(savePath+'gridded_data'+'ice_thickness_'+snowVar+labelStr+str(int(resolution))+version, allow_pickle=True)

	ice_thicknessCS2[where(ice_thicknessCS2<0.25)]=np.nan
	ice_thicknessIS2[where(ice_thicknessIS2<0.25)]=np.nan

	ice_thicknessCS2G=griddata((xptsCS2.flatten(), yptsCS2.flatten()), ice_thicknessCS2.flatten(), (xptsIS2, yptsIS2), method='nearest')
	ice_thicknessCS2G=ma.masked_where(~np.isfinite(ice_thicknessIS2), ice_thicknessCS2G)
	ice_thicknessCS2G=ma.masked_where(~np.isfinite(ice_thicknessCS2G), ice_thicknessCS2G)

	ice_thicknessIS2=ma.masked_where(~np.isfinite(ice_thicknessCS2G), ice_thicknessIS2)
	ice_thicknessIS2=ma.masked_where(~np.isfinite(ice_thicknessIS2), ice_thicknessIS2)

	region_mask, xptsI, yptsI = cF.get_region_mask_sect('../../AncData/', proj, xypts_return=1)
	region_maskG=griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsIS2, yptsIS2), method='nearest')
	regions=[10, 11, 12, 13, 15]
	ice_thicknessIS2=ma.masked_where(~np.isin(region_maskG, regions), ice_thicknessIS2)
	ice_thicknessCS2G=ma.masked_where(~np.isin(region_maskG, regions), ice_thicknessCS2G)


	ice_thicknessIS2M = ice_thicknessIS2.flatten()[ice_thicknessIS2.flatten().mask == False]
	ice_thicknessCS2M = ice_thicknessCS2G.flatten()[ice_thicknessCS2G.flatten().mask == False]



	trend, sig, r_a, intercept = cF.correlateVars(ice_thicknessCS2M, ice_thicknessIS2M)
	rStr = '%.2f' % r_a
	rmse=sqrt(mean((ice_thicknessIS2M-ice_thicknessCS2M)**2))
	rmsStr='%.2f' % rmse


	merr=mean(ice_thicknessCS2M-ice_thicknessIS2M)
	merrStr='%.2f' % merr	

	std=np.std(ice_thicknessIS2M+merr-ice_thicknessCS2M)
	stdStr='%.2f' % std

	minval=0
	maxval=5

	fig = plt.figure(figsize=(8, 2.6))
	plt.subplots_adjust(bottom=0.17, left=0.01, top = 0.9, right=0.99, hspace=0.22, wspace=0.22)

	ax1 = plt.subplot(141, projection=ccrs.NorthPolarStereo(central_longitude=-45))
	
	
	#ax1=axs.flatten()[0]
	sca(ax1)
	#im1 = imshow(ice_thicknessIS2, levels=np.arange(minval, maxval+0.1, 0.25), cmap=cm.viridis, vmin=minval, vmax=maxval, extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	im1 = ax1.pcolormesh(lonG , latG, ice_thicknessIS2, transform=ccrs.PlateCarree(), cmap=cm.viridis, vmin=minval, vmax=maxval, shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	#im1=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessIS2,
	#        cmap=cm.viridis, vmin=vmin, vmax=vmax, zorder=2, rasterized=True)
	#plt.clim(-0.3,5)
	#mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	#mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	#apProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	#mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)
	ax1.set_extent([-179, 179, 63, 90], ccrs.PlateCarree())
	ax1.coastlines('50m')
	ax1.add_feature(cartopy.feature.LAND, facecolor='0.95', edgecolor='0.5')
	ax1.gridlines(draw_labels=False,
         linewidth=0.22, color='gray', alpha=0.5, linestyle='--')

	#cbar.set_ticks(np.arange(0, vmaxs[0]+0.1, 0.2))
	ax1.annotate(monLabel+' '+yearStr+' '+relStr, xy=(0.02, 1.12), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')
	ax1.annotate('(a) IS-2 ('+CS2PRODUCT+' inputs)', xy=(0.02, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

	ax2 = plt.subplot(142, projection=ccrs.NorthPolarStereo(central_longitude=-45))
	sca(ax2)  

	im2 = ax2.pcolormesh(lonG , latG, ice_thicknessCS2G, transform=ccrs.PlateCarree(), cmap=cm.viridis, vmin=minval, vmax=maxval, shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
	ax2.set_extent([-179, 179, 63, 90], ccrs.PlateCarree())
	ax2.coastlines('50m')
	ax2.add_feature(cartopy.feature.LAND, facecolor='0.95', edgecolor='0.5')
	ax2.gridlines(draw_labels=False,
         linewidth=0.22, color='gray', alpha=0.5, linestyle='--')

         #im2=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessCS2G,
	#        cmap=cm.viridis, vmin=vmin, vmax=vmax, zorder=2, rasterized=True)
	#plt.clim(-0.3,5)
	#mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	#mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	#mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	#mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)

	#cax2 = fig.add_axes([0.28, 0.12, 0.2, 0.035])
	#cbar2 = colorbar(im2,cax=cax2, orientation='horizontal', extend='both', use_gridspec=True)
	#cbar2.set_label('CS2 thickness (m)', labelpad=3)
	ax2.annotate('(b) '+CS2PRODUCT+' CS-2', xy=(0.02, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')

	#cbar2.set_ticks(np.arange(0, vmaxs[1]+0.1, 0.1))

	ax3 = plt.subplot(143, projection=ccrs.NorthPolarStereo(central_longitude=-45))
	sca(ax3) 

	im3 = ax3.pcolormesh(lonG , latG, ice_thicknessCS2G-ice_thicknessIS2, transform=ccrs.PlateCarree(), cmap=cm.RdBu_r, vmin=-2, vmax=2, zorder=4, rasterized=True)
	ax3.set_extent([-179, 179, 63, 90], ccrs.PlateCarree())
	ax3.coastlines('50m')
	ax3.add_feature(cartopy.feature.LAND, facecolor='0.95', edgecolor='0.5')
	ax3.gridlines(draw_labels=False,
         linewidth=0.22, color='gray', alpha=0.5, linestyle='--')
	#im3=mapProj.pcolormesh(xptsIS2, yptsIS2,ice_thicknessCS2G-ice_thicknessIS2,
	#        cmap=cm.RdBu_r, vmin=-2, vmax=2, zorder=2, rasterized=True)
	#mapProj.drawcoastlines(linewidth=0.25, zorder=5)
	#mapProj.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	#mapProj.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	#mapProj.fillcontinents(color='0.9',lake_color='grey', zorder=3)
	ax3.annotate('(c) CS-2 minus IS-2', xy=(0.02, 1.02), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left',color='k')


	cax2 = fig.add_axes([0.14, 0.15, 0.2, 0.035])
	cbar2 = colorbar(im2,cax=cax2, orientation='horizontal', extend='both', use_gridspec=True)
	cbar2.set_label('Sea ice thickness (m)', labelpad=3)
	cbar2.set_ticks(np.arange(0, 5.1, 1))
	

	cax3 = fig.add_axes([0.53, 0.15, 0.2, 0.035])
	cbar3 = colorbar(im3,cax=cax3, orientation='horizontal', extend='both', use_gridspec=True)
	cbar3.set_label('difference (m)', labelpad=3)
	cbar3.set_ticks(np.arange(-2, 2.1, 1))
	


	ax4 = plt.subplot(144)
	sca(ax4)
	plt.scatter(ice_thicknessCS2G.flatten(), ice_thicknessIS2.flatten(), marker='x', color='0.2', s=4, alpha=.3)
	#nbins, binEdges, _=plt.hist(elevation, bins=30, linewidth=1.5, histtype='step', color='k', density=True, label='elevation')
	#histVals=binEdges+(binEdges[1]-binEdges[0])
	plt.plot(np.arange(0, 10, 0.1), np.arange(0, 10, 0.1), 'k', ls='-', alpha=.5)

	plt.plot(np.arange(0, 10, 0.1), trend*np.arange(0, 10, 0.1)+intercept, 'k', ls='--', alpha=.8)

	ax4.set_xlabel('CS2 thickness (m)', labelpad=1)
	ax4.set_ylabel('IS2 thickness (m)', labelpad=1)
	#ax4.annotate('(d)\n r: '+rStr+'\nBias: '+merrStr+' m'+'\nSD: '+stdStr+' m', xy=(0.02, 0.98), 
	#		xycoords='axes fraction', color='k', verticalalignment='top', horizontalalignment='left')
	ax4.annotate('(d)', xy=(0.02, 1.01), 
			xycoords='axes fraction', color='k', verticalalignment='bottom', horizontalalignment='left')

	ax4.annotate('r:'+rStr+'\nMB:'+merrStr+' m'+'\nSD:'+stdStr+' m', xy=(0.02, 0.98), 
			xycoords='axes fraction', color='k', fontsize=8, verticalalignment='top', horizontalalignment='left')

	#ax3.annotate('(c) Elevation distribution' , xy=(0., 1.02), xycoords='axes fraction', color='k', horizontalalignment='left', verticalalignment='bottom')
	ax4.set_xlim(0, 6)
	ax4.set_ylim(0, 6)

	fig.savefig(figPath+'/thicknessComp_'+labelStr+relStr+runStr+'CS2IS2corr4'+CS2PRODUCT+version+'3NP.png', dpi=300)


if __name__ == '__main__':
	dates=['201811', '201812', '201901', '201902', '201903', '201904']
	#for month in months:
	#	main('GSFC', month)
	products=[ 'CPOM', 'JPL', 'AWI', 'GSFC']
	for product in products:
		for date in dates:
			main(product, date)
	#months=[11, 12, 1, 2, 3]
	#for month in months:
	#	main('JPL', month)
	#months=[11, 12, 1, 2, 3]
	#for month in months:
	#	main('AWI', month)
	#months=[11]
	#or month in months:
		#main('GSFC', month)

