#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.morphology as morpho
from skimage import morphology
import gdal
import rasterio
import sys, os
from tqdm import tqdm

from rasterio.windows import Window

#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Iterate through an image and calculate the sum of the pixels in a 3x3 square and apply to a new image.
**************************************************************************
##Tasks:
- Read image.
- Calculate the sum of the pixels in a 3x3 box, which itterates across the whole image.
- Save the new classification image for areas scoring 6 and above.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def cumval(x):
    # Cumulative value.
    # Adds up all values.
    return x.sum()

def noise_removal(array, eroded_output):
    # Create an array of zeroes and another of ones.
    arrZero, arrOne = np.zeros(array.shape, dtype=bool), np.ones(array.shape, dtype=bool)
    # Where the array has a value more than 5 make the arrOne = True, and everywhere else (False) = arrZero.
    array_mod = np.where(array>5, arrOne, arrZero)
    array_erode = morphology.remove_small_objects(array_mod, 3, connectivity=2)

    outDataset = gdal.GetDriverByName("GTiff").Create(eroded_output, cols, rows, 1, gdal.GDT_Float32)
    outDataset.SetProjection(proj)
    outDataset.SetGeoTransform(geom)
    outBand = outDataset.GetRasterBand(1)
    outBand.WriteArray(array_erode)

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", nargs="+", required=True, help="Input threshold classification image.").completer = FilesCompleter(allowednames=(".tif"))
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------
        for img in args.input_img:
            if not os.path.exists(img):raise RuntimeError(f"{img} does not exist, please input a true filepath")
        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        print(f"Generating heatmap and eroded and dilated mask for {len(args.input_img)} images.")
        for input in tqdm(args.input_img):
            # Read image
            with rasterio.open(input, 'r') as threshold_class:
                threshold_array = threshold_class.read(1)
                profile = threshold_class.profile
            # Itterate through the image and calculate the sum of all pixels in the 3x3 box.
            # The 'constant' gives the values which lie on the size of the image a value of 0.
            window_filter = ndimage.generic_filter(threshold_array, cumval, size=(3,3), mode='constant')
            # Heatmap saved for future purposes. 
            heatmap_output = os.path.split(input)[0]+"/03b"+os.path.basename(os.path.splitext(input)[0])[3:].rsplit("classification")[0]+"heatmap.tif"
            with rasterio.open(heatmap_output, 'w', **profile) as h_output:
                h_output.write(window_filter, 1)
            img_open = gdal.Open(input)
            cols, rows, proj, geom = img_open.RasterXSize, img_open.RasterYSize, img_open.GetProjection(), img_open.GetGeoTransform()
            # Noise removal from the heatmap.
            nr_output = os.path.split(input)[0]+"/03c"+os.path.basename(os.path.splitext(input)[0])[3:].rsplit("classification")[0]+"heatmap_eroded.tif"
            noise_removal(window_filter, os.path.split(input)[0] + "/03c" + os.path.basename(os.path.splitext(input)[0])[3:].rsplit("classification")[0]+ "heatmap_eroded.tif")
        print("Process complete - heatmap and masks produced.")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
