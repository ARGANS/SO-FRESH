#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import csv
from csv import DictWriter
from datetime import datetime
import difflib
import fiona
import gdal, ogr
import geopandas as gpd
import os, sys
import pandas as pd
from pathlib import Path
import shapely.geometry
from shapely.geometry import Point, Polygon, LineString, mapping
from shapely.ops import nearest_points
import matplotlib.pyplot as plt
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
antartica_mask = gpd.read_file(os.path.split(__file__)[0] + "/antartica_landmask_2deg_buffer.geojson")["geometry"][0]
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

def area_calculator(polygon):
    gdf = gpd.read_file(polygon)
    if not gdf.empty:
        # Take the geometry of the polygons and measure their area, using "cylindrical equal area" as this is what we need to preserve.
        cea = gdf["geometry"].to_crs({"proj":"cea"})
        # Calculate area and get it in km2.
        gdf["Area"] = round(cea.area / 10 ** 6, 2)
        # Write to shapefile.
        gdf.to_file(polygon, driver='GeoJSON')

# This function is not currently in use. 
def generate_polygon_centroid(polygon):
    # Generate a point vector for every polygon present. 
    gdf = gpd.read_file(polygon)
    if not gdf.empty:
        gdf.geometry = gdf.representative_point()
        gdf.to_file(os.path.split(polygon)[0] + "/07" + os.path.basename(os.path.splitext(polygon)[0])[2:] + "_centroid.geojson", driver='GeoJSON')

    return(os.path.split(polygon)[0] + "/07" + os.path.basename(os.path.splitext(polygon)[0])[2:] + "_centroid.geojson")

# This function is not currently in use. 
def distance_calculator(point, landmask):
    # Calculate the distance between a point vector and the nearest point of a specified polygon. This will print the distance in "km".
    gdf_pnt = gpd.read_file(point)
    if not gdf_pnt.empty:
        polygon_lst = []
        centroid_lst = []
        for point_geom in gdf_pnt["geometry"]:
            # This returns the lat lon of the point vector and polygon point which is nearest. 
            polygon, centroid = nearest_points(landmask, point_geom)
            polygon_lst.append(polygon)
            centroid_lst.append(centroid)
        distance_lst = []
        for poly, pnt in zip(polygon_lst, centroid_lst):
            # As the projection for the data is in EPSG:4326 WGS 84 (global projection), the distance is in degrees, where 1 degree == 111 km. 
            distance = round(poly.distance(pnt)*111, 2)
            distance_lst.append(distance)
        gdf_pnt["Distance(KM)"] = distance_lst
        gdf_pnt.to_file(point, driver='GeoJSON')


def selector(shapefile, mask):
    # List of files
    files = os.path.split(shapefile)[0].rsplit('/')
    # The product in the file path - where to save the csv to.
    product = ''.join(difflib.get_close_matches(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[1], files))
    # Filepath reconstructed.
    filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/"
    # Read shapefile
    gdf = gpd.read_file(shapefile)
    outputlist = []
    # If the polygon intersects with the land mask - ignore.
    for i, geom_intersect in enumerate(gdf["geometry"].intersects(mask)):
        if geom_intersect == False:
            # If area is more than X km2, collect file information.
            if len(gdf["Area"] > 200) >= 1:
                # Outputs:
                date = datetime.strptime(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[2][1:], "%Y%j").date()
                date.strftime("%Y-%m-%d")
                version = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[4]
                tile = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[3]
                area = gdf['Area'][i]
                filename = os.path.basename(shapefile).rsplit('_', 4)[0][3:]
                minx = gdf.bounds.iloc[i][0]
                miny = gdf.bounds.iloc[i][1]
                maxx = gdf.bounds.iloc[i][2]
                maxy = gdf.bounds.iloc[i][3]
                outputlist.append([date, version, tile, area, filename, minx, miny, maxx, maxy])
                # Buffer the bounding box by 1.5 (this increases the lat lon by 1.5)
                geometry_gpd = gpd.GeoSeries(Polygon([(minx, miny), (maxx, miny), (maxx, maxy), (minx,maxy)]))
                s = geometry_gpd.buffer(1.5, join_style=2)
                s.to_file(filename='polygon.geojson', driver='GeoJSON')
                sys.exit()
                #fix, axs = plt.subplots(3, 2, figsize=(12, 12), sharex=True, sharey=True)
                #for ax in axs.flatten():
                #    geometry_gpd.plot(ax=ax)
                #    ax.set(xticks=[], yticks=[])
                #s.plot(ax=axs[1, 0], alpha=0.6)
                #axs[1, 0].set_title("s.buffer(0.2, cap_style=2)")
                #plt.show()
                sys.exit()
                
            else:
                pass
        else:
            pass
    return outputlist

def append_data(img, info):
    for i in info:
        files = os.path.split(img)[0].rsplit('/')
        product = ''.join(difflib.get_close_matches(os.path.split(img)[1].rsplit('_', 4)[0].rsplit('.', 7)[1], files))
        filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/"
        # Create output file if it does not exist.
        if not os.path.exists(filepath + "01_csv/"):
            os.mkdir(filepath + "01_csv/")
        else:
            pass
        # Create a new csv with specified headers and insert a row.
        if not os.path.exists(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):
            headers = ["Date", "Version", "Tile", "Area", "Filename", "minx", "miny", "maxx", "maxy"]
            Path(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv").touch()
            with open(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "w+") as f:
                writer = csv.DictWriter(f, fieldnames=headers)
                writer.writeheader()
                writer.writerow({"Date":i[0], "Version":i[1], "Tile": i[2], "Area":i[3], "Filename":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8]})
        # If the file exists, insert the following data in a new row.
        elif os.path.exists(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):
            with open(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "a") as infile:
                headers = ["Date", "Version", "Tile", "Area", "Filename", "minx", "miny", "maxx", "maxy"]
                writer = csv.DictWriter(infile, fieldnames=headers)
                writer.writerow({"Date":i[0], "Version":i[1], "Tile": i[2], "Area":i[3], "Filename":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8]})
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
            ### Create polygon centroids
            ###centroid = generate_polygon_centroid(vector)
            ### Calculate distance of polygons in shapefile from land mask.
            ###distance = distance_calculator(centroid, antartica_mask)
            # Selects those which fit in the criteria (i.e. area <= 200km2 and out of the land mask).
            select = selector(vector, antartica_mask)
            # For files which pass - save them to the csv.
            # Check if the csv exists prior and whether it contains information with similar dates - as this function applies them to the csv as new rows. 
            append = append_data(img, select)
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
