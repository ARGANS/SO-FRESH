#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

# Import packages.
import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cv2
from datetime import datetime, date, time, timezone
import gc, glob
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cm
import matplotlib.ticker as mticker
import numpy as np
import os
import rasterio
import sys

parser = argparse.ArgumentParser(description="""
===================================================
### Geospatial density plot displayed through time. ###
Execution options:
# data-folder
    ~ for data access & saving
# start-date
    ~lower limit of date selection
# end-date 
    ~upper limit of date selection
# products
    ~products for execution
    (ensure all data is downloaded, pre-processed & (optional) fused & run through APIT)
===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Functions ###
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def info_builder(img):

    name = os.path.basename(img)
    
    roi = ((os.path.splitext(img)[0]).rsplit("_", 1)[1]).capitalize()
    date = (os.path.splitext(img)[0]).rsplit("_", 2)[1]
    day,month,year = (date[6:]), (date[4:-2]), (date[:-4])
 
    # Remove the rest of the elements from the list to select the products used.
    spl = name.rsplit("_")
    spl.remove(f"{(os.path.splitext(img)[0]).rsplit('_', 1)[1]}.tif")
    spl.remove(date)
    spl.remove("03")
    spl.remove("ESA")
    spl.remove("SOFRESH")
    spl.remove("APIT")
    

    prods = f"Products: {', '.join(spl)}"
    
    date_f = f"{day}/{month}/{year}"
    dt = datetime.strptime(f"{date_f}", "%d/%m/%Y")

    date = (dt.strftime('%d %B %Y'))
    return(date, prods, roi)


def graph_builder(img, outfile):
    
    date, prods, roi = info_builder(img)

    # Read image
    src = rasterio.open(img)
    arr = src.read(1)
    # Set no-data values to be 0.
    arr = np.where(arr == 1, arr, 0)
    # Mask values of 0.
    masked_data = np.ma.masked_where(arr == 0.0, arr)
 

    # Calculate the sum of each column and row.
    y_arr = np.sum(arr, axis=1)
    x_arr = np.sum(arr, axis=0)

    # Generate flattened arrays for line graph
    y_len = np.asarray([y for y in range(0, y_arr.shape[0])])
    x_len = np.asarray([x for x in range(0, x_arr.shape[0])])

    y_coords = []
    for y in y_len:
        y_coords.append(src.bounds[3]-(y/10))

    x_coords = []
    for x in x_len:
        x_cord = src.bounds[0]+(x/10)
        x_coords.append(src.bounds[0]+(x/10))


    resol = '50m'
    land = cfeature.NaturalEarthFeature('physical', 'land', scale=resol, edgecolor=None, facecolor=cfeature.COLORS['land_alt1'])
    #ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='none', facecolor=cfeature.COLORS['water'])

    fig = plt.figure(figsize=(10,8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[0.5,3.5], height_ratios=[3.5,0.5])

    # Add geospatial plot and format the subplot.
    ax_GEO = fig.add_subplot(gs[0,1], projection=crs.PlateCarree()) # top right
    #ax_GEO.set_extent([src.bounds[0],src.bounds[2],src.bounds[1],-50])
    ax_GEO.set_extent([-180,180,-90,-50])
    ax_GEO.add_feature(land)
    gl = ax_GEO.gridlines(crs=crs.PlateCarree(), draw_labels=True, linewidth=2, color='black', alpha=0.5, linestyle=':')
    gl.xlabels_bottom=False
    gl.ylabels_left = False
    gl.xlines = False
    #gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135, 180])
    gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
    gl.ylocator = mticker.FixedLocator([-55, -70, -85])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 16, 'color': 'black'}#, 'weight': 'bold'}
    gl.ylabel_style = {'size': 16, 'color': 'black'}#, 'weight': 'bold'}

    # Add density plots.
    ax_Y = fig.add_subplot(gs[0,0], sharey=ax_GEO) #top left
    ax_X = fig.add_subplot(gs[1,1], sharex=ax_GEO) # bottom right

    plt.figtext(0.25, 0.95, "The Automated Polynya Identification Tool", size=19, fontfamily="sans-serif", fontweight="bold")
    plt.figtext(0.63, 0.03, date, size=22, fontfamily="sans-serif", fontstyle="oblique", fontweight="bold")
    #plt.figtext(0.35, 0.95, roi, size=17, fontfamily="sans-serif", fontstyle="oblique")
    plt.figtext(0.01, 0.03, prods, size=17, fontfamily="sans-serif", fontstyle="oblique", fontweight="bold")

    # Remove density figure borders.
    ax_Y.spines['top'].set_visible(False)
    ax_Y.spines['right'].set_visible(False)
    ax_Y.spines['bottom'].set_visible(False)
    ax_Y.spines['left'].set_visible(False)
    ax_X.spines['top'].set_visible(False)
    ax_X.spines['right'].set_visible(False)
    ax_X.spines['bottom'].set_visible(False)
    ax_X.spines['left'].set_visible(False)
    # Remove axis labels
    ax_Y.xaxis.set_ticks([])
    ax_Y.yaxis.set_ticks([])
    ax_X.yaxis.set_ticks([])
    ax_X.xaxis.set_ticks([])
    # Remove space between all plots.
    fig.subplots_adjust(wspace=0.0015, hspace=0.0015)

    # Plot data and sort setup.
    ax_GEO.imshow(arr, cmap=cm.Greys, extent=(src.bounds[0],src.bounds[2],src.bounds[1],src.bounds[3]), origin="upper", aspect="auto")
    ax_GEO.imshow(masked_data, cmap="tab20", extent=(src.bounds[0],src.bounds[2],src.bounds[1],src.bounds[3]), origin="upper", aspect="auto")
    ax_X.bar(x_coords, x_arr)
    ax_X.invert_yaxis()
    ax_Y.barh(y_coords, y_arr)
    ax_Y.invert_xaxis()


    #plt.show()
    #sys.exit()
    plt.savefig(outfile, dpi=100)
    del fig, gs, ax_GEO, ax_X, ax_Y, src, arr
    gc.collect()
    #sys.exit()





imgs = sorted(glob.glob("/data/sofresh/03_APIT/SIC/2017/11/*/03*antarctica.tif"))
outputdir = "mp4/SIC/2017/07-11/"

#for img in imgs:
#    graph_builder(img, (f"{outputdir}{os.path.basename(img)[:-4]}.png"))
#sys.exit()


png_imgs = sorted(glob.glob("mp4/SIC/2017/07-11/*png"))
frame = cv2.imread(png_imgs[0])
height, width, layers = frame.shape
video = cv2.VideoWriter(f"{outputdir}APIT_SIC_07-11_2017.avi", 0, 1.75, (width,height))

for img in png_imgs:
    video.write(cv2.imread(img))

cv2.destroyAllWindows()
video.release()
sys.exit()