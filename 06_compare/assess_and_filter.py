#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
from collections import Counter
import csv
import geopandas as gpd
import pandas as pd
from rtree import index
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
    gdf = gpd.GeoDataFrame(df, geometry = 'geometry')
    return(gdf)


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
        print(gdf)
        


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
