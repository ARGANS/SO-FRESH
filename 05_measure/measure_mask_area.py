#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import csv
from csv import DictWriter
import difflib
import gdal
import geopandas as gpd
import os, sys
import pandas as pd
#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
#
**************************************************************************
##Tasks:
- Vectorize mask and calculate the area of the vectors.
- If they meet a condition keep their filename / information.
-
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def vectorize(mask):
    if os.path.exists(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".shp"):
        os.remove(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".cpg")
        os.remove(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".dbf")
        os.remove(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".prj")
        os.remove(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".shp")
        os.remove(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".shx")
        os.system("gdal_polygonize.py %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".shp", "/06" + os.path.basename(os.path.splitext(mask)[0])[2:]))
    else:
        os.system("gdal_polygonize.py %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".shp", "/06" + os.path.basename(os.path.splitext(mask)[0])[2:]))
    return(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".shp")

def area_calculator(shapefile):
    gdf = gpd.read_file(shapefile)
    cea = gdf["geometry"].to_crs({"proj":"cea"})
    gdf['Area'] = cea.area / 10 ** 6
    gdf.to_file(shapefile)

def selector(shapefile):
    # List of files
    files = os.path.split(shapefile)[0].rsplit('/')
    # The product in the file path - where to save the csv to.
    product = ''.join(difflib.get_close_matches(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[2], files))
    # Filepath reconstructed.
    filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/"



    gdf = gpd.read_file(shapefile)
    ds = gdf['Area'] > 200
    if len(ds) >= 1:
        filename = os.path.split(shapefile)[1].rsplit('_', 4)[0][3:]
        doy = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[3][5:]
        year = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[3][1:-3]
        tile = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[4]
        version = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[5]
        minarea = min(gdf['Area'])
        maxarea = max(gdf['Area'])

    return([filename, doy, year, tile, version, minarea, maxarea])



def append_data(data, dataframe):
    df = pd.DataFrame(data=data).T

    test = dataframe.append(df)
    print(df)
    print(test)
    # writing to csv could probably be easier as a column name is not necesarilly required



    """
    #with open(filepath + product[3:] + "_datadownload.csv", "w") as csv:
    with open(filepath + product[3:] + "_datadownload.csv", 'a+', newline='') as write_obj:
        field_names = ["filename", "DOY", "Year", "Tile", "Version", "Minimum Area (km2)", "Maximum Area (km2)"]
        content_dict = {"filename":filename, "DOY":doy, "Year":year, "Tile":tile, "Version":version, "Minimum Area (km2)":minarea, "Maximum Area (km2)":maxarea}

        # Create a writer object from csv module
        dict_writer = DictWriter(write_obj, fieldnames=field_names)
        dict_writer.writeheader()
        # Add dictionary as wor in the csv
        dict_writer.writerow(content_dict)
        '''
        # Create a writer object from csv module
        csv_writer = writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(content)
        '''

    sys.exit()


    r = csv.reader(csv)
    #row = r.next()
    row.append("filename", "DOY", "Year", "Tile", "Version", "Minimum Area (km2)", "Maximum Area (km2)")
    print(row)
    sys.exit()


    print(filename)
    print(doy)
    print(year)
    print(tile)
    print(version)
    print(minarea)
    print(maxarea)
    sys.exit()
    """
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
        df = pd.DataFrame(columns=("filename", "DOY", "Year", "Tile", "Version", "Minimum Area (km2)", "Maximum Area (km2)"))
        for img in args.input_img:
            vector = vectorize(img)
            area = area_calculator(vector)
            #print(area)
            select = selector(vector)
            save = append_data(select, df)
        print(df)


        #print(select)
        #--------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
