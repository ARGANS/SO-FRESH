#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import gdal
import geopandas as gpd
from osgeo import ogr
import os, sys
#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
#
**************************************************************************
##Tasks:
-
-
-
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def vectorize(mask):
    os.system("gdal_polygonize.py %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".shp", "/06" + os.path.basename(os.path.splitext(mask)[0])[2:]))
    return(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".shp")

def area_calculator(shapefile):
    gdf = gpd.read_file(shapefile)
    gdf = gdf["geometry"].to_crs({"proj":"cea"})
    area = gdf.area / 10 ** 6
    area.to_file(shapefile)
    #return(area)

def selector(shapefile):
    gdf = gpd.read_file(shapefile)
    return(gdf)

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", nargs="+", required=True, help="").completer = FilesCompleter(allowednames=(".tif"))
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        for img in args.input_img:
            vector = vectorize(img)
            area = area_calculator(vector)
            #print(area)
            select = selector(vector)
            print(select)
        #--------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
