#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import gdal
import glob
import sys, os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import osr
from tqdm import tqdm

#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="""
# Produces a projected tif file from the MODIS TCI images found on: https://e4ftl01.cr.usgs.gov/
**************************************************************************
##Tasks:
- Extract coordinates from the MODIS sinsusoidal tiles based on the file name.
- Project array of data using extracted coordinates.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def modis_extract_geom(input, coords):
    # Extract geometry from the sinusoidal tiles reference based on those reference in the filename. 
    if os.path.basename(input).rsplit(".", 7)[3][0] == 'h' and os.path.basename(input).rsplit(".", 7)[3][3] == 'v':
        h, v = os.path.basename(input)[28:-27], os.path.basename(input)[31:-24]
    else:
        raise RuntimeError("The filename does not match up.")

    df = pd.read_csv(coords)
    df_row = df[(df['iv'] == int(v)) & (df['ih'] == int(h))]

    if len(gdal.Open(input).ReadAsArray().shape) == 3:
        ny, nx = gdal.Open(input).ReadAsArray().shape[1], gdal.Open(input).ReadAsArray().shape[2]
    elif len(gdal.Open(input).ReadAsArray().shape) == 2:
        ny, nx = gdal.Open(input).ReadAsArray().shape[0], gdal.Open(input).ReadAsArray().shape[1]
    else:
        raise RuntimeError('Please input a 2D or 3D array.')

    xmin, ymin, xmax, ymax = float(df_row['lon_min']), float(df_row['lat_min']), float(df_row['lon_max']), float(df_row['lat_max'])
    xres = (xmax - xmin) / nx
    yres = (ymax - ymin) / ny

    return(xmin, xres, 0, ymax, 0, -yres)

def modis_jpg2tif(input, output, epsg, geom):
    # Convert image from a JPEG to a TIFF.
    img = gdal.Open(input).ReadAsArray()

    if len(gdal.Open(input).ReadAsArray().shape) == 3:
        shp, ny, nx = gdal.Open(input).RasterCount, gdal.Open(input).ReadAsArray().shape[1], gdal.Open(input).ReadAsArray().shape[2]
    elif len(gdal.Open(input).ReadAsArray().shape) == 2:
        shp, ny, nx = gdal.Open(input).RasterCount, gdal.Open(input).ReadAsArray().shape[0], gdal.Open(input).ReadAsArray().shape[1]
    else:
        raise RuntimeError('Please input a 2D or 3D array.')

    #Set output information (stored from input image)
    outDataset = gdal.GetDriverByName("GTiff").Create(output, ny, nx, shp, gdal.GDT_Byte)
    outDataset.SetGeoTransform(geom)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    outDataset.SetProjection(srs.ExportToWkt())

    count = 0
    for arr in img:
        count += 1
        outBand = outDataset.GetRasterBand(count)
        outBand.WriteArray(arr)
        outBand = None

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", required=True, nargs = "+", help="The input jpg file. The output will be output to this directory with the same filename with a 'tif' extension.").completer = FilesCompleter(allowednames=(".jpg"))
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        print(f"Reprojecting {len(args.input_img)} images.")
        for img in tqdm(args.input_img):
            modis_jpg2tif(img, os.path.split(img)[0] + "/02" + os.path.basename(os.path.splitext(img)[0])[2:] + ".tif", 4326, modis_extract_geom(img, os.path.split(__file__)[0]+"/modis_sinusoidal_tiles.csv"))
        print("Process complete - all images projected.")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
