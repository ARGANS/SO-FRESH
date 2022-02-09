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
    if input.endswith(".jpg"):
        if os.path.basename(input).rsplit(".")[3][0] == "h" and os.path.basename(input).rsplit(".")[3][3] == "v":
            h, v = os.path.basename(input).rsplit(".")[3][1:3], os.path.basename(input).rsplit(".")[3][4:6]
        else:
            raise RuntimeError("The filename does not match up.")
    elif input.endswith(".tif"):
        if os.path.basename(input).rsplit(".")[2][0] == "h" and os.path.basename(input).rsplit(".")[2][3] == "v":
            h, v = os.path.basename(input).rsplit(".")[2][1:3], os.path.basename(input).rsplit(".")[2][4:6]

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

def modis_extract_hdf(img):
    # Extract bands from MODIS MYDTGBA products.
    # Preprocessed directory output.
    pp_dir = os.path.split(os.path.abspath(img))[0] + ("/02a_" + os.path.basename(os.path.splitext(img)[0])[3:])
    if not os.path.isdir(pp_dir):
        os.mkdir(pp_dir)
    cmd = "gdal_translate -q -of GTIFF -sds %s %s"%(img, (pp_dir + "/02b_" + os.path.basename(os.path.splitext(img)[0])[3:] + ".tif"))
    os.system(cmd)
    return (pp_dir + "/02b_" + os.path.basename(os.path.splitext(img)[0])[3:])

def normalise_reprjct(img, epsg, geom):
    # Normalise MODIS imagery from storage values and assign projection. (Specifically MODIS Band 32 - Thermal Infrared (TIR) 11.77-12.27 um).
    cmd = (gdal.Open(img).ReadAsArray()) * 0.01
    output = os.path.dirname(os.path.dirname(img)) + ("/02_" + os.path.basename(img)[4:-6] + ".tif")
    shp, ny, nx = gdal.Open(img).RasterCount, gdal.Open(img).ReadAsArray().shape[0], gdal.Open(img).ReadAsArray().shape[1]
    outDataset = gdal.GetDriverByName("GTiff").Create(output, ny, nx, shp, gdal.GDT_Byte)
    outDataset.SetGeoTransform(geom)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    outDataset.SetProjection(srs.ExportToWkt())
    outBand = outDataset.GetRasterBand(1)
    outBand.SetDescription("Band_32")
    outBand.WriteArray(cmd)
    outband = None

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", required=True, nargs = "+", help="The input file. The output will be output to this directory with the same filename with a 'tif' extension.").completer = FilesCompleter(allowednames=(".jpg", ".hdf"))
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        print("Processing all entries...")
        for img in tqdm(args.input_img):
            if img.endswith(".jpg"):
                modis_jpg2tif(img, os.path.split(img)[0] + "/02" + os.path.basename(os.path.splitext(img)[0])[2:] + ".tif", 4326, modis_extract_geom(img, os.path.split(__file__)[0]+"/modis_sinusoidal_tiles.csv"))
            elif img.endswith(".hdf"):
                extracted = modis_extract_hdf(img)
                for e in sorted(glob.glob(extracted+"*")):
                    if e.endswith( "_4.tif"):
                        normalise_reprjct(e, 4326, modis_extract_geom(e, os.path.split(__file__)[0]+"/modis_sinusoidal_tiles.csv"))
        print("Processing complete.")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
