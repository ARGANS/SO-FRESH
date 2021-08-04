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
    if os.path.exists(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson"):
        os.remove(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson")
        os.system("gdal_polygonize.py %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson", "/06" + os.path.basename(os.path.splitext(mask)[0])[2:]))
    else:
        os.system("gdal_polygonize.py %s %s -b 1 %s"%(mask, os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson", "/06" + os.path.basename(os.path.splitext(mask)[0])[2:]))

    return(os.path.split(mask)[0] + "/06" + os.path.basename(os.path.splitext(mask)[0])[2:] + ".geojson")

def area_calculator(shapefile):
    gdf = gpd.read_file(shapefile)
    cea = gdf["geometry"].to_crs({"proj":"cea"})
    gdf['Area'] = cea.area / 10 ** 6
    gdf.to_file(shapefile, driver='GeoJSON')

def selector(shapefile):
    # List of files
    files = os.path.split(shapefile)[0].rsplit('/')
    # The product in the file path - where to save the csv to.
    product = ''.join(difflib.get_close_matches(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[2], files))
    # Filepath reconstructed.
    filepath = '/'.join(files[0:files.index(product)]) + "/" + product + "/"
    # Read shapefile
    gdf = gpd.read_file(shapefile)
    outputlist = []
    # If area is more than X km2, collect file information.
    if len(gdf['Area'] > 200) >= 1:
        # This itterates through if there are > 1 polygon that meets the criteria.
        for i in range(len(gdf.bounds)):
            date = datetime.strptime(os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[3][1:-3] + "-" + os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[3][5:], "%Y-%j").strftime("%Y-%m-%d")
            tile = os.path.split(shapefile)[1].rsplit('_', 4)[0].rsplit('.', 7)[4]
            area = gdf['Area'][i]
            filename = os.path.basename(shapefile).rsplit('_', 4)[0][3:]
            minx = gdf.bounds.iloc[i][0]
            miny = gdf.bounds.iloc[i][1]
            maxx = gdf.bounds.iloc[i][2]
            maxy = gdf.bounds.iloc[i][3]
            outputlist.append([date, tile, area, filename, minx, miny, maxx, maxy])

    # Here I have written the following outputs to be inputted into a list. This essentially itterates through polygon and extract specific information and saves each one as a new list. (so list in a list).
    # if you do a print(outputlist) it'll be very clear.
    # Things to do:
    # - Code so the csv is being saved in a suitable location - lines 47 to 51 basically do that, where it originally looks for the product name in the filename and looks for it in the path.
    # - This may need to be done out of the function and done with the code below. (This should be done on the VM, so it works and is solid.)
    # - Where there are lists in lists, each element should be written to a new column every time and when it is complete, it should then go to a new line.
    # - Header for CSV are on line 148.



    '''
    with open('test_csv.csv', 'w') as f1:
        writer=csv.writer(f1, delimiter='\t',lineterminator='\n',)
        for elem in outputlist:
            for i in elem:
                writer.writerow(map(lambda x: [x], i))

        sys.exit()
    '''
    #return([date, tile, area, filename, minx, miny, maxx, maxy])



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
        df = pd.DataFrame(columns=("Date", "Tile", "Area", "Filename", "minx", "miny", "maxx", "maxy"))
        df.to_csv('test_csv.csv')
        for img in args.input_img:
            vector = vectorize(img)
            area = area_calculator(vector)
            select = selector(vector)
            print(select)
            #save = append_data(select, df)
        #print(df)


        #print(select)
        #--------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
