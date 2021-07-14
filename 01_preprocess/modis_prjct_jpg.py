#!/usr/bin/env python

import argparse
import gdal
import sys, os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import osr
#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Produces a projected tif file from the MODIS TCI images found on: https://e4ftl01.cr.usgs.gov/
**************************************************************************
##Tasks:
- Extract coordinates from the MODIS sinsusoidal tiles based on the file name.
- Project array of data using extracted coordinates.
-
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def modis_extract_geom(input, coords):
    if os.path.basename(input)[24] == 'h' and os.path.basename(input)[27] == 'v':
        h, v = os.path.basename(input)[25:-27], os.path.basename(input)[28:-24]
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
    xres = (xmax - xmin) / float(nx)
    yres = (ymax - ymin) / float(ny)
    return(xmin, xres, 0, ymin, 0, -yres)

def modis_jpg2tif(input, output, epsg, geom):
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

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", required=True, help="The input jpg file. The output will be output to this directory with the same filename with a 'tif' extension.")
        parser.add_argument("-c", "--coords-csv", required=True, help="MODIS Sinusodial tiles CSV file.")
        parser.add_argument("-e", "--epsg", required=False, default=4326, help="The EPSG to set the image to, default is 4326.")
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        modis_jpg2tif(args.input_img, os.path.splitext(args.input_img)[0] + ".tif", args.epsg, modis_extract_geom(args.input_img, args.coords_csv))

        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
