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
def cumval(x):
    # Cumulative value.
    # Adds up all values.
    return x.sum()

def erosion_dialation(array, eroded_output):
    # Create an array of zeroes and another of ones.
    arrZero, arrOne = np.zeros(array.shape), np.ones(array.shape)
    # Where the array has a value more than 5 make the arrOne = True, and everywhere else (False) = arrZero.
    array_mod = np.where(array>5, arrOne, arrZero)
    #profile.update({'dtype':'uint8'})
    profile.update({'nodata':'0'})
    # Set structure - on how the array is being checked.
    structure = ndimage.generate_binary_structure(2, 2)
    # Erode, dialate and erode again.
    print(f"Eroding & dilating {os.path.split(input)[1]}...")
    array_mod = morpho.binary_erosion(array_mod, structure=structure, iterations=1, border_value=0)
    array_mod = morpho.binary_dilation(array_mod, structure=structure, iterations=2, border_value=0)
    array_mod = morpho.binary_erosion(array_mod, structure=structure, iterations=1, border_value=0)

    # Select the appropriate data.
    array_mod = np.where(array_mod,1,0)
    array_mod = np.array(array_mod, dtype=np.uint8)
    with rasterio.open(eroded_output, 'w', **profile) as imgout:
        imgout.write(array_mod, 1)
        

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
            # Read image
            with rasterio.open(input, 'r') as threshold_class:
                threshold_array = threshold_class.read(1)
                profile = threshold_class.profile
            # Itterate through the image and calculate the sum of all pixels in the 3x3 box.
            # The 'constant' gives the values which lie on the size of the image a value of 0.
            window_filter = ndimage.generic_filter(threshold_array, cumval, size=(3,3), mode='constant')
            
            # Heatmap saved for future purposes. 
            heatmap_output = os.path.join(os.path.split(input)[0], "04" + os.path.basename(os.path.splitext(input)[0])[2:] + "_heatmap.tif")
            with rasterio.open(heatmap_output, 'w', **profile) as h_output:
                h_output.write(window_filter, 1)
            
            # Erode and dialate the heatmap.
            erosion_dialation(window_filter, os.path.split(input)[0] + "/05" + os.path.basename(os.path.splitext(input)[0])[2:] + "_heatmap_eroded.tif")
            
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
