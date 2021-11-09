#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
from datetime import datetime, timedelta
import gdal
import numpy as np
from PIL import Image
import sys, os, glob

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
def img_to_array(img):
    return (gdal.Open(img).ReadAsArray())
    

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", nargs="+", help="").completer = FilesCompleter(allowednames=(".tif"))
        parser.add_argument("-d", "--days", type=int, help="Number of days to look at from input image.")
        parser.add_argument("-t", "--tile", help="Filepath to MODIS tile of interest.")
        parser.add_argument("-s", "--time-start", help="Time lower bound (YYYY/MM/DD)")
        parser.add_argument("-e", "--time-end", help="Time upper bound (YYYY/MM/DD)")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------
        #if args.tile ends with '/' remove it!
        # time-end or days must be specified!
        #make sure date is split with "/" or if it is as '-' then change it?
        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        tile = args.tile
        sdate = datetime.strptime(os.path.join(args.time_start), "%Y/%m/%d").date()
        # End date based on specified or calculated from the start date from the number of days. 
        if args.time_end:
            edate = datetime.strptime(os.path.join(args.time_end), "%Y/%m/%d").date()
        elif args.days:
            edate = sdate + timedelta(days=args.days)
        #else:
        #    raise RuntimeError("Please specify either an end-date using '-e' or a number of days from the start date using '-d'.")

        # Acquire list of all the dates.
        date_list = [sdate + timedelta(days=x) for x in range((edate-sdate).days + 1)]
        # Full directory (including tile & date).
        directory_dates = ["/".join((tile, str(d.year), str('%02d' %d.month), str('%02d' %d.day))) for d in date_list if os.path.exists("/".join((tile, str(d.year), str('%02d' %d.month), str('%02d' %d.day))))]
        # List of directory and filenames.
        full_filepath = [files for dir in directory_dates for files in glob.glob(os.path.join(dir, "08*.tif"))]
        # Extract array of images within the date range that are available.
        array = np.array([np.array(Image.open(fimg)) for fimg in full_filepath])
        cumulative_array = np.sum(np.stack(array), axis=0)
        
        print(array.shape)
        sys.exit()

        
        
        
        print(tile)
        print(time_start)
        print(time_end)
        print("".join((tile, time_start)))
        sys.exit()
        
        for img in args.input_img:
            # Split path to get date.
            path_split = os.path.split(img)[0].rsplit('/')
            # Extract components of date to be fed into datetime. 
            day, month, year = path_split[-1], path_split[-2], path_split[-3]
            # Filepath to the tile which is being investigated.
            tile_filepath = "/".join(path_split[0:-3])
            
            # Acquire date range.            
            start_date = datetime.strptime(os.path.join(year,month,day), "%Y/%m/%d").date()
            end_date = start_date + timedelta(days=args.days)
            # Acquire list of all the dates.
            date_list = [start_date + timedelta(days=x) for x in range((end_date-start_date).days + 1)]
            # Generate filepaths for all the dates.
            dir_date_list = ["/".join((tile_filepath, str(d.year), str('%02d' %d.month), str('%02d' %d.day))) for d in date_list if os.path.exists("/".join((tile_filepath, str(d.year), str('%02d' %d.month), str('%02d' %d.day))))]
            # Pull files for all dates.     
            files_date_list = [files for dir in dir_date_list for files in glob.glob(os.path.join(dir, "08*.tif"))]
            
            array = np.array([np.array(Image.open(fimg)) for fimg in files_date_list])
            cumulative_array = np.sum(np.stack(array), axis=0)
            print(array.shape)
            sys.exit()
            



            # write to image for viewing purposes
            cols = gdal.Open(img).RasterXSize
            rows = gdal.Open(img).RasterYSize
            proj = gdal.Open(img).GetProjection()
            geom = gdal.Open(img).GetGeoTransform()

            outDataset = gdal.GetDriverByName("GTiff").Create("test.tif", cols, rows, 1, gdal.GDT_Float32)
            outDataset.SetProjection(proj)
            outDataset.SetGeoTransform(geom)
            outBand = outDataset.GetRasterBand(1)
            outBand.WriteArray(cumulative_array) 
            
            sys.exit()
            
                    
            
            sys.exit()
            #array = img_to_array(img)
            #print(array)
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
