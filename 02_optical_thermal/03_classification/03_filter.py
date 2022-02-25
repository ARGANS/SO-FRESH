#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import csv
from csv import DictWriter
from datetime import datetime
import difflib
import gdal, ogr
import geopandas as gpd
import os, sys
from pathlib import Path
from shapely.ops import unary_union
from tqdm import tqdm
#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
#
**************************************************************************
##Tasks:
- Vectorize mask and calculate the area of the vectors.
- If they meet specific condition (specific area (km2) and are >200km from land) return information including:date, version, tile, total area, full filename and bounding box values.
- For all files which pass save the information to a csv.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Hard arguments:
#----------------------------------------------------------------------------------------------------
# Antarctica land mask buffered by 2 degrees ~ 222 km. 1 degree buffer also available.
antarctica_mask = gpd.read_file(os.path.split(__file__)[0] + "/03_landmask/antartica_landmask_2deg_buffer.geojson")["geometry"][0]
#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def vectorize(mask):
    # See if the file already exists, if not, create it and if so, remove it and create it.
    if os.path.exists(os.path.split(mask)[0] + "/03d" + os.path.basename(os.path.splitext(mask)[0])[3:] + ".geojson"):
        os.remove(os.path.split(mask)[0] + "/03d" + os.path.basename(os.path.splitext(mask)[0])[3:] + ".geojson")
        os.system("gdal_polygonize.py -q -f 'GeoJSON' %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/03d" + os.path.basename(os.path.splitext(mask)[0])[3:] + ".geojson", "/03d" + os.path.basename(os.path.splitext(mask)[0])[3:]))
    else:
        os.system("gdal_polygonize.py -q -f 'GeoJSON' %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/03d" + os.path.basename(os.path.splitext(mask)[0])[3:] + ".geojson", "/03d" + os.path.basename(os.path.splitext(mask)[0])[3:]))

    return(os.path.split(mask)[0] + "/03d" + os.path.basename(os.path.splitext(mask)[0])[3:] + ".geojson")

def area_calculator(polygon):
    gdf = gpd.read_file(polygon)
    if not gdf.empty:
        # Take the geometry of the polygons and measure their area, using "cylindrical equal area" as this is what we need to preserve.
        cea = gdf["geometry"].to_crs({"proj":"cea"})
        # Calculate area and get it in km2.
        gdf["Area"] = round(cea.area / 10 ** 6, 2)
        # Write to shapefile.
        gdf.to_file(polygon, driver='GeoJSON')

def filter(shapefile, mask):
    files = os.path.split(shapefile)[0].rsplit('/') # List of files
    #product = ''.join(difflib.get_close_matches(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[1], files)) # The product in the file path - where to save the csv to.
    #filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/" # Filepath reconstructed.
    gdf = gpd.read_file(shapefile) # Read shapefile
    if not gdf.empty:
        #If the polygon intersects with the land mask - ignore.
        # If polygon is smaller than 200 km2 - ignore.
        for index, row in gdf.iterrows():
            if row["geometry"].intersects(mask) or row["Area"] > 85000 or not row["Area"] > 200:
                gdf.drop(index, inplace=True)
        if not gdf.empty:
            trimmed_shp = os.path.split(shapefile)[0] + "/03e" + os.path.basename(os.path.splitext(shapefile)[0])[3:] + "_criteria.geojson"
            if os.path.exists(trimmed_shp):
                os.remove(trimmed_shp)
                gdf.to_file(trimmed_shp, driver="GeoJSON")
            else:
                gdf.to_file(trimmed_shp, driver="GeoJSON")

            return trimmed_shp
        else:
            pass
    else:
        pass

def mask(shapefile, img):
    xmin, xpixel, _, ymax, _, ypixel = gdal.Open(img).GetGeoTransform()
    xmax = xmin + xpixel * gdal.Open(img).RasterXSize
    ymin = ymax + ypixel * gdal.Open(img).RasterYSize
    os.system("gdalwarp -q -overwrite --config GDALWARP_IGNORE_BAD_CUTLINE YES -of GTiff -te %s %s %s %s -tr %s %s -tap -cutline %s -cl %s %s %s"%(xmin, ymin, xmax, ymax, xpixel, ypixel, shapefile, os.path.splitext(os.path.basename(shapefile))[0], img, os.path.split(img)[0] + "/03f" + os.path.basename(os.path.splitext(shapefile)[0])[3:] + "_masked.tif"))

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
        for img in tqdm(args.input_img):
            # Vectorize image.
            vector = vectorize(img)
            # Calculate area of polygons in shapefile.
            area = area_calculator(vector)
            # Selects those which fit in the criteria (i.e. area <= 200km2 and do not intersect land mask).
            if os.path.isfile(vector):
                criteria = filter(vector, antarctica_mask)
            # Remove areas which have not fit the criteria by masking them out of the vector files.
            if not criteria == None:masking = mask(criteria, img)
            else:pass
        print("Process complete.")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
