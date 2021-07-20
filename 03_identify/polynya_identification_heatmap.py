#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.morphology as morpho
import gdal
import rasterio
import sys, os

#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Iterate through an image and calculate the sum of the pixels in a 3x3 square and apply to a new image.
**************************************************************************
##Tasks:
- Read image.
- Calculate the sum of the pixels in a 3x3 box, which itterates across the whole image.
- Save the new classification image.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def cumsum(x):
    return x.sum()


def erosion_dialation(heatmap, eroded_output):
    # Read image.
    array = rasterio.open(heatmap).read(1)
    # Create an array of zeroes and another of ones.
    arrZero, arrOne = np.zeros(array.shape), np.ones(array.shape)
    # Where the array has a value more than 5 make the arrOne = True, and everywhere else (False) = arrZero.
    array_mod = np.where(array>5, arrOne, arrZero)
    profile = rasterio.open(heatmap).profile
    profile.update({'dtype':'uint8'})
    profile.update({'nodata':'0'})
    # Set structure - on how the array is being checked.
    structure = ndimage.generate_binary_structure(2, 2)
    # Erode, dialate and erode again.
    print("Eroding & dilating...")
    array_mod = morpho.binary_erosion(array_mod, structure=structure, iterations=1, border_value=0)
    array_mod = morpho.binary_dilation(array_mod, structure=structure, iterations=2*1,border_value=0)
    array_mod = morpho.binary_erosion(array_mod, structure=structure,iterations=1, border_value=0)
    # Select the appropriate data.
    array_mod = np.where(array_mod,1,0)
    with rasterio.open(eroded_output, 'w', **profile) as imgout:
        imgout.write(array_mod.astype('uint8'), 1)

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
        for input in args.input_img:
            # Itterate through the image and calculate the sum of all pixels in the 3x3 box.
            # The 'constant' gives the values which lie on the size of the image a value of 0.
            window_filter = ndimage.generic_filter(gdal.Open(input).ReadAsArray(), cumsum, size=(3,3), mode='constant')
            #Set output information (stored from input image)
            driver = gdal.GetDriverByName("GTiff")
            outImg = driver.Create(os.path.split(input)[0] + "/04" + os.path.basename(os.path.splitext(input)[0])[2:] + ".heatmap.tif", gdal.Open(input).RasterYSize, gdal.Open(input).RasterXSize, 1, gdal.GDT_Byte)
            outBand = outImg.GetRasterBand(1)
            outBand.WriteArray(window_filter)
            outImg.SetGeoTransform(gdal.Open(input).GetGeoTransform())
            outImg.SetProjection(gdal.Open(input).GetProjection())
            del window_filter, driver, outImg, outBand

            # Erode and dialate the heatmap.
            erosion_dialation(os.path.split(input)[0] + "/04" + os.path.basename(os.path.splitext(input)[0])[2:] + ".heatmap.tif", os.path.split(input)[0] + "/05" + os.path.basename(os.path.splitext(input)[0])[2:] + ".heatmap.eroded.tif")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
