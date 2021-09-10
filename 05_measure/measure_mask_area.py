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
import itertools
import pandas as pd
import numpy as np
from pathlib import Path
import shapely.geometry
from shapely.geometry import Point, Polygon, LineString, mapping
from shapely.ops import nearest_points, unary_union
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
antarctica_mask = gpd.read_file(os.path.split(__file__)[0] + "/antartica_landmask_2deg_buffer.geojson")["geometry"][0]
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


def area_distance_criteria(shapefile, mask):
    # Check to see whether the polygons reach the initial criteria (area and distance from land).
    # Read shapefile
    gdf = gpd.read_file(shapefile)
    shp_list = []
    gdf_list = []
    # If the polygon intersects with the land mask - pass
    for i, geom_intersect in enumerate(gdf["geometry"].intersects(mask)):
        if geom_intersect == False:
            # If area is more than X km2, collect file information.
            if len(gdf["Area"] > 200) >= 1:
                if shapefile not in shp_list: # Avoids appending multiple of the same filename. 
                    shp_list.append(shapefile)
                    gdf_list.append(gdf)
            elif len(gdf["Area"] < 200) >=1: # if the area is below 200 km2 - remove polygon
                gdf = gdf.drop(labels=i, axis=0)
            else:
                pass
        
        elif geom_intersect == True:
            gdf = gdf.drop(labels=i, axis=0) # If there is an intersection - remove polygon
        else:
            pass
        
    return(shp_list, gdf_list)

def bbox_extract(geodataframe):
    # Extract the bounding boxes of all polygons that have passed the initial criteria. 
    bbox_geom_lst = []
    for i, val in enumerate(geodataframe["geometry"]):
        minx = geodataframe.bounds.iloc[i][0]
        miny = geodataframe.bounds.iloc[i][1]
        maxx = geodataframe.bounds.iloc[i][2]
        maxy = geodataframe.bounds.iloc[i][3]
        bbox_geom_lst.append([minx, miny, maxx, maxy]) 
    
    return(bbox_geom_lst)


def selector(shapefile, mask):
    # List of files
    files = os.path.split(shapefile)[0].rsplit('/')
    # The product in the file path - where to save the csv to.
    product = ''.join(difflib.get_close_matches(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[1], files))
    # Filepath reconstructed.
    filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/"
    # Read shapefile
    gdf = gpd.read_file(shapefile)
    
    #If the polygon intersects with the land mask - ignore.
    # If polygon is smaller than 200 km2 - ignore.
    for index, row in gdf.iterrows():
        if row["geometry"].intersects(mask) or not row["Area"] > 200:
            gdf.drop(index, inplace=True)
    
    # Buffer the bounding box of the polygons.
    gdf_buffered = gdf.envelope.buffer(0.75, join_style=2)
    # Union all those that overlap.
    union = gdf_buffered.unary_union
    # Turn the union in to a GeoSeries.
    shapes_series = gpd.GeoSeries(union, crs=gdf.crs)
    # If the file is a 'MultiPolygon' the explode function splits it into individual polygons.
    exploded = shapes_series.explode()
    #exploded.to_file(filename='test_v1.geojson', driver='GeoJSON')
    outputlist = []
    for index, row in exploded.bounds.iterrows():
        minx, miny, maxx, maxy = row[0], row[1], row[2], row[3]
        date = datetime.strptime(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[2][1:], "%Y%j").date()
        date.strftime("%Y-%m-%d")
        version = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[4]
        tile = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[3]
        filename = os.path.basename(shapefile).rsplit('_', 4)[0][3:]
        quantity_poly = len(exploded.bounds)
        outputlist.append([date, version, tile, filename, quantity_poly, row[0], row[1], row[2], row[3]])
    else:
        pass

    return outputlist



    ##################################################
    #- Append all items in shapes_series bounding boxes to csv.
    #################################################


    # Order: minx, miny, maxx, maxy
    if len(overlapping_bbox_lst) >= 1:
        full_bbox = [(min(list(list(zip(*overlapping_bbox_lst))[0]))), (min(list(list(zip(*overlapping_bbox_lst))[1]))), (max(list(list(zip(*overlapping_bbox_lst))[2]))), (max(list(list(zip(*overlapping_bbox_lst))[3])))]
        # The next two lines give the functionality for the bounding box to be saved as geojson. 
        #full_bbox_gdf = gpd.GeoSeries(Polygon([(full_bbox[0], full_bbox[1]), (full_bbox[2], full_bbox[1]), (full_bbox[2], full_bbox[3]), (full_bbox[0],full_bbox[3])]))
        #full_bbox_gdf.to_file(filename='polygon.geojson', driver='GeoJSON')
        
        #############
        date = datetime.strptime(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[2][1:], "%Y%j").date()
        date.strftime("%Y-%m-%d")
        version = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[4]
        tile = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[3]
        filename = os.path.basename(shapefile).rsplit('_', 4)[0][3:]
        quantity_poly = len(overlapping_bbox_lst)
        outputlist.append([date, version, tile, filename, quantity_poly, full_bbox[0], full_bbox[1], full_bbox[2], full_bbox[3]])

    else:
        pass
    
    return outputlist
   
def append_data(img, info):
    # If wanting to search for solo polygons.
    '''
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
            headers = ["Date", "Version", "Tile", "Filename", "minx", "miny", "maxx", "maxy"]
            Path(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv").touch()
            with open(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "w+") as f:
                writer = csv.DictWriter(f, fieldnames=headers)
                writer.writeheader()
                writer.writerow({"Date":i[0], "Version":i[1], "Tile": i[2], "Area": i[3], "Filename":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8]})
        # If the file exists, insert the following data in a new row.
        elif os.path.exists(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):
            with open(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "a") as infile:
                headers = ["Date", "Version", "Tile", "Filename", "minx", "miny", "maxx", "maxy"]
                writer = csv.DictWriter(infile, fieldnames=headers)
                writer.writerow({"Date":i[0], "Version":i[1], "Tile": i[2], "Area": i[3], "Filename":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8]})
        else:
            pass
    '''
    if not len(info) == 0:
        for i in info:
            # Extract necesary file information to where 
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
                headers = ["Date", "Version", "Tile", "Filename", "Number of Polygons", "minx", "miny", "maxx", "maxy"]
                Path(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv").touch()
                with open(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "w+") as f:
                    writer = csv.DictWriter(f, fieldnames=headers)
                    writer.writeheader()
                    writer.writerow({"Date":i[0], "Version":i[1], "Tile":i[2], "Filename":i[3], "Number of Polygons":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8]})
            # If the file exists, insert the following data in a new row.
            elif os.path.exists(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv"):
                with open(str(filepath + "01_csv/" + img.rsplit(".", 6)[1][1:-3]) + "_imgs_in_criteria.csv", "a") as infile:
                    headers = ["Date", "Version", "Tile", "Filename", "Number of Polygons", "minx", "miny", "maxx", "maxy"]
                    writer = csv.DictWriter(infile, fieldnames=headers)
                    writer.writerow({"Date":i[0], "Version":i[1], "Tile":i[2], "Filename":i[3], "Number of Polygons":i[4], "minx":i[5], "miny":i[6], "maxx":i[7], "maxy":i[8]})

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
        for img in args.input_img:
            # Vectorize image.
            vector = vectorize(img)
            # Calculate area of polygons in shapefile.
            area = area_calculator(vector)
            '''
            ### Create polygon centroids
            ###centroid = generate_polygon_centroid(vector)
            ### Calculate distance of polygons in shapefile from land mask.
            ###distance = distance_calculator(centroid, antartica_mask)
            '''
            # Selects those which fit in the criteria (i.e. area <= 200km2 and do not intersect land mask).
            select = selector(vector, antarctica_mask)
            # For files which pass - save them to the csv.
            # Check if the csv exists prior and whether it contains information with similar dates - as this function applies them to the csv as new rows. 
            append = append_data(img, select)
            
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
