#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
from datetime import datetime, timedelta
import gdal
import numpy as np
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
        parser.add_argument("-i", "--input-img", nargs="+", required=True, help="").completer = FilesCompleter(allowednames=(".tif"))
        parser.add_argument("-d", "--days", required=True, type=int, help="Number of days to look at from input image.")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
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
            files_date_list = ["/".join((fp, i)) for fp in dir_date_list for i in glob.glob(os.path.join(fp, "08*.tif"))]
            
            # PRINT ALL VALUES FOR THE DATE RANGE
            # CHECK 05 FILTER script with "writing to an empty dataframe"
            print(files_date_list)
            sys.exit()
            
            for fp in file_date_list:
                for i in glob.glob(os.path.join(fp, "08*")):
                    print(i)
                    sys.exit()
                   
                
            
            
            
            
            sys.exit()
            #array = img_to_array(img)
            #print(array)
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
