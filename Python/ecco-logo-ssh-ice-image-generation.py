#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 16:47:35 2020

PROGRAM GENERATES BACKGROUND IMAGES FOR ECCO LOGO

PROGRAM V4
2020-04-10

IAN FENTY
IAN.FENTY@JPL.NASA.GOV

"""

# ensure the following packages are installed
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import matplotlib as mp


######################################################
# define a directory that has ssh, sea ice area, and sea ice thickness fields
# interpolated to lat-lon grid. This is ECCO v4r4

# SSHDYN_1997_01.nc
# SIarea_1997_01.nc
# SIheff_1997_01.nc


data_dir = '/Users/ifenty/inSync Share/Projects/ECCOv4/logo/v14/SSH and sea-ice data/'
output_dir = '/Users/ifenty/tmp'


# open files

ssh_ds = xr.open_dataset(data_dir + '/SSHDYN_1997_01.nc')
ssh = ssh_ds.SSHDYN.where(ssh_ds.SSHDYN != 0)[0,:]
land_mask = np.ones(ssh.shape)
land_mask = np.where(np.isnan(ssh),1,0)


sia_ds = xr.open_dataset(data_dir + '/SIarea_1997_01.nc')
sia = sia_ds.SIarea.where(sia_ds.SIarea > .2)[0,:]

sih_ds = xr.open_dataset(data_dir + '/SIheff_1997_01.nc')
sih = sih_ds.SIheff.where(sih_ds.SIheff >  .2)[0,:]


# define levels for SSH contours
levels=np.linspace(-1.8,1.2,30)



# pull arrays from datasets, create land mask, 

lons_a = ssh.longitude.values
lats_a = ssh.latitude.values
ssh_a = ssh.values
sia_a = sia.values
sih_a = sih.values

ssh_a, lons_a = add_cyclic_point(ssh_a, lons_a)
sia_a = add_cyclic_point(sia_a)
sih_a = add_cyclic_point(sih_a)

land_mask_a = add_cyclic_point(land_mask)

land_mask_a2 = np.where(land_mask_a == 1, 1, np.nan)


######################################################
#
# PART 1:  COLOR FILL CONTOURS
#
#%%######################################################

plt.close('all')
fig = plt.figure(num=1, figsize=(14,7));
ax = plt.axes(projection=ccrs.Sinusoidal(central_longitude=-60))

vmax = 2.3 
vmax = 2.0

mp.rcParams['axes.linewidth'] = .05
mp.rcParams['lines.linewidth'] = .05
mp.rcParams['patch.linewidth'] = .05


cnt = plt.contourf(lons_a, lats_a, np.zeros(land_mask_a2.shape), levels=2,\
               transform=ccrs.PlateCarree(),vmin=-1,vmax=0,colors='black')

# This is the fix for the white lines between contour levels
for c in cnt.collections:
    c.set_edgecolor("face")


cnt = plt.contourf(lons_a, lats_a, ssh_a,\
               transform=ccrs.PlateCarree(),levels=levels,
               cmap=plt.get_cmap('nipy_spectral'),
               vmin=-2, vmax=vmax)

# This is the fix for the white lines between contour levels
for c in cnt.collections:
    c.set_edgecolor("face")

    
cnt = plt.contourf(lons_a, lats_a, sih_a/sia_a, \
             transform=ccrs.PlateCarree(),vmin=0,vmax=1.5, levels=40,
             cmap=cmocean.cm.ice, zorder=100)

# This is the fix for the white lines between contour levels
for c in cnt.collections:
    c.set_edgecolor("face")

cnt = plt.contourf(lons_a, lats_a, land_mask_a2, levels=20,\
             transform=ccrs.PlateCarree(),vmin=0,vmax=1,cmap='binary',
             linewidths=0.01)

# This is the fix for the white lines between contour levels
for c in cnt.collections:
    c.set_edgecolor("face")

#ax.coastlines('110m', linewidth=.1, zorder=101)
#ax.add_feature(cfeature.LAND, color='black')
ax.gridlines(color='black', linestyle='--',linewidth=0.5, zorder=110)


plt.savefig(output_dir + '/S5_color_fill_h.svg',format='svg')
plt.savefig(output_dir + '/S5_color_fill_h.pdf',format='pdf')
plt.savefig(output_dir + '/S5_color_fill_h.eps',format='eps')



######################################################
#
# PART 2:  COLOR CONTOURS
#
#%%######################################################

plt.close('all')
fig = plt.figure(num=2, figsize=(14,7));
ax = plt.axes(projection=ccrs.Sinusoidal(central_longitude=-60))

vmax = 2.3 # same as logo v 11
vmax = 2.0


mp.rcParams['axes.linewidth'] = 1
mp.rcParams['lines.linewidth'] = 1
mp.rcParams['patch.linewidth'] = 1


cnt = plt.contour(lons_a, lats_a, ssh_a,\
               transform=ccrs.PlateCarree(),levels=levels,
               cmap=plt.get_cmap('nipy_spectral'),
               vmin=-2, vmax=vmax,linewidths=.5)

ax.add_feature(cfeature.LAND, color='gray',zorder=100)
ax.coastlines('110m', linewidth=1, zorder=101, color='black')
ax.gridlines(color='black', linestyle='--',linewidth=0.5, zorder=110)

plt.savefig(output_dir + '/S5_colored_contours_h.svg',format='svg')
plt.savefig(output_dir + '/S5_colored_contours_h.pdf',format='pdf')
plt.savefig(output_dir + '/S5_colored_contours_h.eps',format='eps')



######################################################
#
# PART 3:  BLACK AND WHITE CONTOURS
#
#%%######################################################

plt.close('all')
fig = plt.figure(num=2, figsize=(14,7));
ax = plt.axes(projection=ccrs.Sinusoidal(central_longitude=-60))

vmax = 2.3 # same as logo v 11
vmax = 2.0


mp.rcParams['axes.linewidth'] = 1
mp.rcParams['lines.linewidth'] = 1
mp.rcParams['patch.linewidth'] = 1


cnt = plt.contour(lons_a, lats_a, ssh_a,\
               transform=ccrs.PlateCarree(),levels=levels,
               vmin=-2, vmax=vmax,linewidths=.5, colors='black')

ax.add_feature(cfeature.LAND, color='gray',zorder=100)
ax.coastlines('110m', linewidth=1, zorder=101, color='black')
ax.gridlines(color='black', linestyle='--',linewidth=0.5, zorder=110)

plt.savefig(output_dir + '/S5_bw_contours_h.svg',format='svg')
plt.savefig(output_dir + '/S5_bw_contours_h.pdf',format='pdf')
plt.savefig(output_dir + '/S5_bw_contours_h.eps',format='eps')

