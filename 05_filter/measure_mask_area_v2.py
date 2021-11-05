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
antarctica_mask = gpd.read_file(os.path.split(__file__)[0] + "/antartica_landmask_2deg_buffer.geojson")["geometry"][0]
#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def vectorize(mask):
    # See if the file already exists, if not, create it and if so, remove it and crete it.
    if os.path.exists(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson"):
        os.remove(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson")
        os.system("gdal_polygonize.py -q -f 'GeoJSON' %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson", "/06" + os.path.basename(os.path.splitext(mask)[0])[2:]))
    else:
        os.system("gdal_polygonize.py -q -f 'GeoJSON' %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson", "/06" + os.path.basename(os.path.splitext(mask)[0])[2:]))

    return(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson")

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
    product = ''.join(difflib.get_close_matches(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[1], files)) # The product in the file path - where to save the csv to.
    filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/" # Filepath reconstructed.
    gdf = gpd.read_file(shapefile) # Read shapefile
   
    #If the polygon intersects with the land mask - ignore.
    # If polygon is smaller than 200 km2 - ignore.
    for index, row in gdf.iterrows():
        if row["geometry"].intersects(mask) or row["Area"] > 85000 or not row["Area"] > 200:
            gdf.drop(index, inplace=True)

    trimmed_shp = os.path.split(shapefile)[0] + "/07" + os.path.basename(os.path.splitext(shapefile)[0])[2:] + "_criteria.geojson"
    gdf.to_file(trimmed_shp, driver="GeoJSON")
    
    outputlist = []
    # Extract required information. 
    for area, geometry in zip(gdf["Area"], gdf["geometry"]):
        date = datetime.strptime(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[2][1:], "%Y%j").date()
        date.strftime("%Y-%m-%d")
        doy = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[2][1:][4:]
        year = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[2][1:][:-3]
        tile = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[3]
        
        outputlist.append([date, doy, year, geometry, area, tile])
    else:
        pass
    return outputlist, trimmed_shp

def mask(shapefile, img):
    xmin, xpixel, _, ymax, _, ypixel = gdal.Open(img).GetGeoTransform() 
    os.system("gdalwarp -q -overwrite --config GDALWARP_IGNORE_BAD_CUTLINE YES -of GTiff -tr %s %s -tap -cutline %s -cl %s -crop_to_cutline %s %s"%(xpixel, ypixel, shapefile, os.path.splitext(os.path.basename(shapefile))[0], img, os.path.split(img)[0] + "/08" + os.path.basename(os.path.splitext(shapefile)[0])[2:] + "_masked.tif"))



def append_data_by_year(img, info):
    # This function appends information in to csv files based on years.
    if not len(info) == 0:
        for i in info:
            # Extract necesary file information to where the csv is to be saved.
            files = os.path.split(img)[0].rsplit('/')
            product = ''.join(difflib.get_close_matches(os.path.split(img)[1].rsplit('_', 4)[0].rsplit('.', 7)[1], files))
            filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/"
            # Create output file if it does not exist.
            if not os.path.exists(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3] + "/"):
                os.mkdir(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3] + "/")
            else:
                pass
            # Create a new csv with specified headers and insert a row.
            if not os.path.exists(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):
                headers = ["Date", "Version", "Tile", "Filename", "Number of Polygons", "minx", "miny", "maxx", "maxy", "geometry"]
                Path(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv").touch()
                with open(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "w+") as f:
                    writer = csv.DictWriter(f, fieldnames=headers)
                    writer.writeheader()
                    writer.writerow({"Date":i[0], "Version":i[1], "Tile":i[2], "Filename":i[3], "Number of Polygons":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8], "geometry":i[9]})
            # If the file exists, insert the following data in a new row.
            elif os.path.exists(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):
                with open(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "a") as infile:
                    headers = ["Date", "Version", "Tile", "Filename", "Number of Polygons", "minx", "miny", "maxx", "maxy", "geometry"]
                    writer = csv.DictWriter(infile, fieldnames=headers)
                    writer.writerow({"Date":i[0], "Version":i[1], "Tile":i[2], "Filename":i[3], "Number of Polygons":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8], "geometry":i[9]})

            else:
                pass
    else:
        pass

def append_data_by_month(img, info):
    # This function appends information in to csv files based on years.
    if not len(info) == 0:
        for i in info:
            # Extract necesary file information to where the csv is to be saved.
            files = os.path.split(img)[0].rsplit('/')
            product = ''.join(difflib.get_close_matches(os.path.split(img)[1].rsplit('_', 4)[0].rsplit('.', 7)[1], files))
            filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/"
            # Create output file if it does not exist.
            if not os.path.exists(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3] + "/"):
                os.mkdir(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3] + "/")
            else:
                pass
            # Create a new csv with specified headers and insert a row.
            if not os.path.exists(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + files[:-1][-1] + "_" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):
                headers = ["Date", "Version", "Tile", "Filename", "Number of Polygons", "minx", "miny", "maxx", "maxy", "geometry"]
                Path(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + files[:-1][-1] + "_" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv").touch()
                with open(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + files[:-1][-1] + "_" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "w+") as f:
                    writer = csv.DictWriter(f, fieldnames=headers)
                    writer.writeheader()
                    writer.writerow({"Date":i[0], "Version":i[1], "Tile":i[2], "Filename":i[3], "Number of Polygons":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8], "geometry":i[9]})
            # If the file exists, insert the following data in a new row.
            elif os.path.exists(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + files[:-1][-1] + "_" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):
                with open(str(filepath + "01_csv/" +  img.rsplit(".", 6)[1][1:-3] + "/" + files[:-1][-1] + "_" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "a") as infile:
                    headers = ["Date", "Version", "Tile", "Filename", "Number of Polygons", "minx", "miny", "maxx", "maxy", "geometry"]
                    writer = csv.DictWriter(infile, fieldnames=headers)
                    writer.writerow({"Date":i[0], "Version":i[1], "Tile":i[2], "Filename":i[3], "Number of Polygons":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8], "geometry":i[9]})

            else:
                pass
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
        print(f"Examining the criteria of {len(args.input_img)} masks and appending those which fit to CSV.")
        for img in tqdm(args.input_img):
            # Vectorize image.
            vector = vectorize(img)
            # Calculate area of polygons in shapefile.
            area = area_calculator(vector)
            # Selects those which fit in the criteria (i.e. area <= 200km2 and do not intersect land mask).
            criteria = filter(vector, antarctica_mask)
            # Remove areas which have not fit the criteria by masking them out of the vector files.
            masking = mask(criteria[1], img)
            # For files which pass - save them to the csv by year.
            #append_by_year = append_data_by_year(img, select)
            # For files which pass - save them to the csv by month in corresponding year folders.
            #append_by_month = append_data_by_month(img, select)
        print("Process complete.")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
