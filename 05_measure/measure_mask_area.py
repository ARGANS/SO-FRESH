#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import csv
from csv import DictWriter 
from datetime import datetime
import difflib
import gdal
import geopandas as gpd
import os, sys
import pandas as pd
from pathlib import Path
#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
#
**************************************************************************
##Tasks:
- Vectorize mask and calculate the area of the vectors.
- If they meet a condition keep their filename / information in a csv.
-
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def vectorize(mask):
    # See if the file already exists, if not, create it and if so, remove it and crete it.
    if os.path.exists(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson"):
        os.remove(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson")
        os.system("gdal_polygonize.py %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson", "/06" + os.path.basename(os.path.splitext(mask)[0])[2:]))
    else:
        os.system("gdal_polygonize.py %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson", "/06" + os.path.basename(os.path.splitext(mask)[0])[2:]))

    return(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson")

def area_calculator(shapefile):
    gdf = gpd.read_file(shapefile)
    # Take the geometry of the polygons and measure their area, using "cylindrical equal area" as this is what we need to preserve. 
    cea = gdf["geometry"].to_crs({"proj":"cea"})
    # Calculate area and get it in km2.
    gdf['Area'] = cea.area / 10 ** 6
    # Write to shapefile.
    gdf.to_file(shapefile, driver='GeoJSON')

def selector(shapefile):
    # List of files
    files = os.path.split(shapefile)[0].rsplit('/')
    # The product in the file path - where to save the csv to.
    product = ''.join(difflib.get_close_matches(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[1], files))
    # Filepath reconstructed.
    filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/"
    # Read shapefile
    gdf = gpd.read_file(shapefile)
    outputlist = []
    # If area is more than X km2, collect file information.
    if len(gdf['Area'] > 200) >= 1:
        # This itterates through if there are > 1 polygon that meets the criteria.
        for i in range(len(gdf.bounds)):
            # Outputs:
            date = datetime.strptime(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[2][1:], "%Y%j").date()
            date.strftime("%Y-%m-%d")
            tile = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[4]
            area = gdf['Area'][i]
            filename = os.path.basename(shapefile).rsplit('_', 4)[0][3:]
            minx = gdf.bounds.iloc[i][0]
            miny = gdf.bounds.iloc[i][1]
            maxx = gdf.bounds.iloc[i][2]
            maxy = gdf.bounds.iloc[i][3]
            outputlist.append([date, tile, area, filename, minx, miny, maxx, maxy])
    return outputlist
    
def append_data(img, info):
    for i in info:
        files = os.path.split(img)[0].rsplit('/')
        product = ''.join(difflib.get_close_matches(os.path.split(img)[1].rsplit('_', 4)[0].rsplit('.', 7)[1], files))
        filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/"
        if not os.path.exists(filepath + "01_csv/"):
            os.mkdir(filepath + "01_csv/")
        else:
            pass
        if not os.path.exists(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):
            headers = ["Date", "Tile", "Area", "Filename", "minx", "miny", "maxx", "maxy"]
            Path(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv").touch()
            with open(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "w+") as f:
                writer = csv.DictWriter(f, fieldnames=headers)
                writer.writeheader() 
                writer.writerow({"Date":i[0], "Tile":i[1], "Area":i[2], "Filename":i[3], "minx":i[4], "miny":i[5], "maxx":i[6], "maxy":i[7]})

        elif os.path.exists(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):            
            with open(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "a") as infile:
                headers = ["Date", "Tile", "Area", "Filename", "minx", "miny", "maxx", "maxy"]
                writer = csv.DictWriter(infile, fieldnames=headers)
                writer.writerow({"Date":i[0], "Tile":i[1], "Area":i[2], "Filename":i[3], "minx":i[4], "miny":i[5], "maxx":i[6], "maxy":i[7]})
        else:
            pass


#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", nargs="+", required=True, help="Eroded and dilated heatmap image.").completer = FilesCompleter(allowednames=(".tif"))
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        for img in args.input_img:
            # Vectorize image.
            vector = vectorize(img)
            # Calculate area of polygons in shapefile.
            area = area_calculator(vector)
            # Selects those which fit in the criteria (i.e. area <= 200km2).
            select = selector(vector)
            append = append_data(img, select)        
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
