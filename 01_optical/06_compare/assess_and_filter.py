#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import cartopy.crs as ccrs
from collections import Counter
import csv
from datetime import datetime  
from datetime import timedelta
from functools import reduce
import geopandas as gpd
import geoplot as gplt
import geoplot.crs as gcrs
import itertools
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
from pyproj import CRS
from rtree import index
import seaborn as sns
from shapely import wkt
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
def csv_2_gdf(input_csv):
    # Read csv to pandas dataframe, and applying the geometry column to give coordinates in the WKT format. 
    df = pd.read_csv(input_csv)
    df['geometry'] = df['geometry'].apply(wkt.loads)
    gdf = gpd.GeoDataFrame(df, geometry="geometry", crs=4326)
    return(gdf)

def centroid(gdf):
    return(gdf.geometry.centroid)
    

def overlap_check(gdf):
    max_gdf_date = datetime.strptime(max(gdf["Date"]), "%Y-%m-%d")
    delta = max_gdf_date - datetime.strptime(min(gdf["Date"]), "%Y-%m-%d")
    i = 0
    j = 13
    intersect_gdf = gpd.GeoDataFrame(columns=["Start Date", "End Date", "geometry"], crs=gdf.crs)
    for days in range(delta.days)[:3]:
        date_start = datetime.strptime(min(gdf["Date"]), "%Y-%m-%d") + timedelta(days=i)
        date_end = datetime.strptime(min(gdf["Date"]), "%Y-%m-%d") + timedelta(days=j)
        if not date_end == max_gdf_date:
            i = i + 1
            j = j + 1
            date_range_gdf = gdf[gdf["Date"].between(str(date_start), str(date_end))]
            '''
            data = []
            for index1, row1 in date_range_gdf.iterrows():
                for index2, row2 in date_range_gdf.iterrows():
                    if row1["geometry"].intersects(row2["geometry"]):
                        #owdspd=row2['id']
                        data.append({'geometry':row1['geometry'].intersection(row2['geometry'])})
            df = gpd.GeoDataFrame(data, columns=["geometry"])
            df.to_file("intersection_test.json", driver="GeoJSON")
            sys.exit()
            '''
            
            overlay = gpd.overlay(date_range_gdf, date_range_gdf, how="intersection", keep_geom_type=False)
            
            #matching_filename_indexes = overlay[overlay["Filename_1"] == overlay["Filename_2"]].index
            

            matching_minmax_xy = overlay[(overlay["Filename_1"] == overlay["Filename_2"]) & 
                                        (overlay["minx_1"] == overlay["minx_2"]) & 
                                        (overlay["miny_1"] == overlay["miny_2"]) & 
                                        (overlay["maxx_1"] == overlay["maxx_2"]) & 
                                        (overlay["maxy_1"] == overlay["maxy_2"])].index
            overlay.drop(matching_minmax_xy, inplace=True)
            overlay.to_file(f"test_overlay_{str(i)}.json", driver="GeoJSON")

            union = overlay.unary_union
            shapes_series = gpd.GeoSeries(union, crs=gdf.crs)
            print(shapes_series)

            shapes_series.to_file(f"test_union_{str(i)}.json", driver="GeoJSON")

            #intersect_gdf.append([date_start, date_end, mp], ignore_index=True)
    
    sys.exit()
    
    
    '''
    i = 0
    j = 0
    for iterr in combinations(gdf.iterrows(), 2):
        print(iterr)
        if iterr[0][1][9].intersects(iterr[1][1][9]):
            i = i +1
            print(i)
        else:
            j = j +1
            print(j)
    sys.exit()
        
    crs =  pyproj.CRS.from_user_input("EPSG:4326")
    series = gpd.GeoSeries(iterr, crs=crs)

    
    for index, row in gdf.iterrows():
        print(row.geometry)
        sys.exit()

        #test = np.where(row.geometry.overlaps(row.geometry, align=False), True, False)
        test = (row.geometry).overlaps(row.geometry, align=False)
        #print(row.geometry)
        print(test)
        sys.exit()
    '''
    '''
    data_overlaps = gpd.GeoDataFrame(crs=gdf.crs)
    for index, row in gdf.iterrows():
        gdf1=gdf.loc[gdf.Filename!=row.Filename,]
        overlaps = gdf1[gdf1.geometry.overlaps(row.geometry)]["Filename"].tolist()
        if len(overlaps) > 0:
            for y in overlaps:
                # Tests for the intersection between individual rows.
                intersection_area=gpd.overlay(gdf.loc[gdf.Filename==y,],gdf.loc[gdf.Filename==row.Filename,],how="intersection")
                # Produces dataframe with one row, containing all information of both geometries which intersect. 
                intersection_area=intersection_area.loc[intersection_area.geometry.area>=9e-9]
                if intersection_area.shape[0] > 0: # Check to see if there is more than one row. 
                    data_overlaps=gpd.GeoDataFrame(pd.concat([intersection_area,data_overlaps],ignore_index=True),crs=gdf.crs)
                    data_overlaps['sorted']=data_overlaps.apply(lambda y: sorted([y["Filename_1"],y["Filename_2"]]),axis=1)
                    data_overlaps['sorted']=data_overlaps.sorted.apply(lambda y: ''.join(y))
                    data_overlaps=data_overlaps.drop_duplicates('sorted')
                    data_overlaps=data_overlaps.reset_index()[['Filename_1','Filename_2','geometry']] # Final GeoDataFrame with all possible combinations of intersections.
    data_overlaps.to_csv("test.csv")
    print(data_overlaps)
    print(type(data_overlaps))
    #sys.exit()
    '''

def heatmap(gdf):
    antarctica_mask = gpd.read_file("github_jhickson/SO-FRESH/05_filter/antartica_landmask.geojson")[['geometry']]
    
    # Add date range.
    #####################################################################
    gdf = gdf[gdf["Date"].between("2017-09-01", "2017-11-30")]
    #print(gdf)
    #gdf = gdf[gdf["Tile"].between("h14v15", "h24v15")]
    #print(gdf)
    
    #####################################################################
    # Rather than a KDE plot, could a sns.jointplot be used, where the localities are specified on a map, and the frequency is displayed on the X and Y axis


    # This produces a plot of kernel density across the whole region.
    ax = gplt.polyplot(antarctica_mask, projection=gcrs.PlateCarree(), facecolor="White", edgecolor="Black", linewidth=0.1, zorder=1)
    polynyas = centroid(gdf)
    gplt.kdeplot(polynyas, ax=ax, cmap='cividis', shade=True, thresh=0, extent=gdf['geometry'].total_bounds)
    #gplt.pointplot(polynyas, ax=ax, extent=antarctica_mask.total_bounds)
    ax.set(xlabel="Longitude (degrees)", ylabel="Latitude (degrees)")
    gridliner = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gridliner.ylabels_left = False
    gridliner.xlabels_top = False
    plt.show()
    
    tiles = gdf["Tile"].unique()
    print(tiles)
    centroid_list = []
    


def month_selector(csvfile):
    month_range = []
    with open(csvfile) as file:
            reader = csv.reader(file)
            next(reader, None)
            for index, elem in enumerate(reader):
                month_range.append(elem[0][5:-3])
            date_range = sorted(set(month_range))
    return date_range

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-csv", required=True, help="Input csv file containg all details of areas identified as polynyas.").completer = FilesCompleter(allowednames=(".csv"))
        #parser.add_argument("-m", "--input-mask", required=True, help="").completer = FilesCompleter(allowednames=(".jpg"))
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        gdf = csv_2_gdf(args.input_csv)
        #centroid(gdf)
        heatmap(gdf)
        #overlap_check(gdf)
        #print(gdf)
        


        '''
        month_range = month_selector(args.input_csv)
        #Things to be done:
        # Smash today - you're awesome.
        # Itterate through loop, where if the date is the same i.e. '09' and '09', then keep going. But figure out the loop order, and is this even the right way. 
        # Once thats figured - itterate over and see which bounding boxes intersect. You should be able to figure out which is which based on their location in the list. 
        # Could always save the information aswell, such as date, for easy reference. 

        with open(args.input_csv) as file:
            reader = csv.reader(file)
            next(reader, None)
            i = 1
            test = []
            #for row, month in zip(reader, month_range):
            #for row in reader:
            for month in month_range:
                for row in reader:
                    print(month)
                    print(row)
                    sys.exit()
                    if month == row[0][5:-3]:
                        minx = float(row[5])
                        miny = float(row[6])
                        maxx = float(row[7])
                        maxy = float(row[8])
                        bbox = (minx, miny, maxx, maxx)
                        idx = index.Index()
                        idx.insert(i, bbox)
                        i = i+1
                        results = []
                        for internal in reader:
                            iminx = float(internal[5])
                            iminy = float(internal[6])
                            imaxx = float(internal[7])
                            imaxy = float(internal[8])
                            ibbox = (iminx, iminy, imaxx, imaxy)
                            #print(bbox)
                            #print(ibbox)
                            #ibbox = ', '.join((iminx, iminy, imaxx, imaxy))
                            results.append(idx.count(ibbox))
                        test.append(results)
                    else:
                        continue
            print(results)
            print(len(results))
            print(test)
            print(len(test))
            sys.exit()
        '''

                   

        '''
                if month == row[0][5:-3]:
                    print('yay')
                    
                    minx = float(row[5])
                    miny = float(row[6])
                    maxx = float(row[7])
                    maxy = float(row[8])
                    bbox = (minx, miny, maxx, maxx)
                    #print(bbox)
                    #bbox = ', '.join((minx, miny, maxx, maxy))
                    idx = index.Index()
                    idx.insert(i, bbox)
                    i = i+1
                    results = []
                    for internal in reader:
                        iminx = float(internal[5])
                        iminy = float(internal[6])
                        imaxx = float(internal[7])
                        imaxy = float(internal[8])
                        ibbox = (iminx, iminy, imaxx, imaxy)
                        #print(bbox)
                        #print(ibbox)
                        #ibbox = ', '.join((iminx, iminy, imaxx, imaxy))
                        results.append(idx.count(ibbox))
                    test.append(results)
                    
                else:
                    print('nay')
                '''
            #print(results)
            #print(len(results))
            #print(test)
            #print(len(test))
        sys.exit()
                    
                    
        # Read csv information.
        # Take monthly intervals - that should be possible to be set. Reading dates
        # Take bounding box - assess the overlap of bounding boxes.
        # This will likely be a loop in a loop, where you iterate through each bounding box and compare to every other one.
        # (If dates / BB / all other info match = skip)


        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
