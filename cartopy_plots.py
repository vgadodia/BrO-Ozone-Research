import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import numpy as np
import cv2
import matplotlib.path as mpath
import cartopy.feature
import matplotlib.colors as clr
from cartopy.util import add_cyclic_point

from pylab import *

# hexes = cm.get_cmap('nipy_spectral', 8)    # PiYG
# arr = []
# for i in range(hexes.N):
#     rgb = hexes(i)[:3] # will return rgba, we take only first 3 so we get rgb
#     arr.append(matplotlib.colors.rgb2hex(rgb))

# # print(arr)
def main(filename):
    rootgrps = []

    file = filename
    rootgrps.append(netCDF4.Dataset(file))
    rootgrp = rootgrps[0]

    # Function to Plot Nitrogen Dioxide
    def plot_ozone(ax, lat, lon, ozone, filename, vmin=None, vmax=None, title=None, cbar=True):
        # new_cmap = clr.LinearSegmentedColormap.from_list('blues', arr, N=256)
        plt.pcolormesh(lon, lat, ozone, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), cmap=plt.get_cmap("gist_rainbow_r"))
        plt.title(title)
        if cbar:
            cbar = plt.colorbar(shrink=0.7, orientation="horizontal", ticks=[100,200,300,400,500,600]);
            cbar.set_label('Column Amount O3 (Dobson Units)')
        # plt.show()
        plt.savefig(filename)

    # Function to Plot Ozone based on Column Amount filter
    def masked_ozone_func(ax, filter_val_less, filter_val_greater, ozone, lat, lon, vmin, vmax, filename, title=None, cbar=True):
        # masked_ozone = np.ma.masked_less(ozone, filter_val_less)
        # masked_ozone = np.ma.masked_greater(ozone, filter_val_greater)
        masked_ozone = np.ma.masked_outside(ozone, filter_val_less, filter_val_greater)
        print("Mean:", masked_ozone.mean())
        print("Median:", np.ma.median(masked_ozone))
        plot_ozone(ax, lat, lon, masked_ozone, title=title, vmin=vmin, vmax=vmax, filename=filename, cbar=cbar)

    lat = rootgrp.variables['latitude'][:]
    lon = rootgrp.variables['longitude'][:]
    ozone = rootgrp.groups["key_science_data"].variables["ColumnAmountO3"][:]

    # lat = np.ma.masked_less(lat, 30)

    ozone, lon = add_cyclic_point(ozone, coord=lon)

    # Number of pixels masked, loop through days
    # Date, mean, median, SD, maximum, minimum, percentiles, Range, IQR, area

    vmin = 150
    vmax = 550

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    # # PLOT GLOBAL OZONE
    # plt.figure(figsize=(10, 10))
    # ax = plt.axes(projection=ccrs.PlateCarree())
    # ax.coastlines(resolution='110m')
    # # ax.gridlines()
    # plot_ozone(ax, lat, lon, ozone, vmin=vmin, vmax=vmax, filename="ozone_global{}.png".format(file[20:24]), cbar=False)

    # PLOT ARCTIC OZONE HOLE
    plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.Stereographic(central_latitude=90.0, central_longitude=0.0, scale_factor=7))
    ax.coastlines(resolution='110m')
    # ax.gridlines()
    ax.set_boundary(circle, transform=ax.transAxes)
    plot_ozone(ax, lat, lon, ozone, vmin=vmin, vmax=vmax, filename="ozone_arctic{}.png".format(file[15:19] + "-" + file[20:24]), cbar=False)
    plt.close()
    
    # # PLOT MASKED ARCTIC OZONE HOLE
    # plt.figure(figsize=(10, 10))
    # ax = plt.axes(projection=ccrs.Stereographic(central_latitude=90.0, central_longitude=0.0, scale_factor=5))
    # ax.coastlines(resolution='110m')
    # # ax.gridlines()
    # ax.set_boundary(circle, transform=ax.transAxes)
    # low = 0
    # high = 220
    # masked_ozone_func(ax, low, high, ozone, lat, lon, vmin=vmin, vmax=vmax, filename="ozone_masked{}.png".format(file[15:19] + "-" + file[20:24]), cbar=False)

    # # PLOT BINARY IMAGE
    # B = ozone < 220
    # B = B.astype(np.int)
    # plt.figure(figsize=(10, 10))
    # ax = plt.axes(projection=ccrs.Stereographic(central_latitude=90.0, central_longitude=0.0, scale_factor=5))
    # ax.coastlines(resolution='110m')
    # # ax.gridlines()
    # ax.set_boundary(circle, transform=ax.transAxes)
    # plot_ozone(ax, lat=lat, lon=lon, ozone=B, filename="ozone_binary{}.png".format(file[20:24]), cbar=False)


    # image = cv2.imread('ozone_binary.png')
    # gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    # (thresh, blackAndWhiteImage) = cv2.threshold(gray, 127, 255, cv2.THRESH_BINARY)
    # cv2.imshow('gray', blackAndWhiteImage)
    # cv2.imwrite('ozone_bw.png', blackAndWhiteImage) 

    # ret, thresh = cv2.threshold(gray, 127, 255, 0)
    # contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    # cont = cv2.drawContours(image, contours, -1, (0,255,0), 3)
    # cv2.imshow('draw contours',cont)

    # cv2.waitKey()
    # cv2.destroyAllWindows()

# main("NPP_SAO_O3T_L3_2020m0101_v001-2020m0708t202411.nc")

# main("NPP_SAO_O3T_L3_2014m0127_v001-2020m0715t030633.nc")
main("NPP_SAO_O3T_L3_2020m0317_v001-2020m0715t023445.nc")
i = 0
import glob, os
for file in glob.glob("*.nc"):
    # main("SAO_NPP_O3T_L3_20190101.nc")
    # if "NPP_SAO_O3T_L3_2020" in file:
    #     main(file)
    pass
# print(2278/7)