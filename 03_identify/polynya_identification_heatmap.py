#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import numpy as np
import scipy.ndimage as ndimage
import gdal
import sys

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

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", required=True, help="Input classification image.").completer = FilesCompleter(allowednames=(".tif"))
        parser.add_argument("-o", "--output-img", required=True, help="Output heatmap image.").completer = FilesCompleter(allowednames=(".tif"))
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        # Read image as array.
        array = gdal.Open(args.input_img).ReadAsArray()
        # Itterate through the image and calculate the sum of all pixels in the 3x3 box.
        # The 'constant' gives the values which lie on the size of the image a value of 0.
        result = ndimage.generic_filter(array, cumsum, size=(3,3), mode='constant')
        #Set output information (stored from input image)
        driver = gdal.GetDriverByName("GTiff")
        outImg = driver.Create(args.output_img, gdal.Open(args.input_img).RasterYSize, gdal.Open(args.input_img).RasterXSize, 1, gdal.GDT_Byte)
        outBand = outImg.GetRasterBand(1)
        outBand.WriteArray(result)
        outImg.SetGeoTransform(gdal.Open(args.input_img).GetGeoTransform())
        outImg.SetProjection(gdal.Open(args.input_img).GetProjection())

        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
