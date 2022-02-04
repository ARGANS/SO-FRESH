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
# Specified date range enables the assessment of open-water appearance frequency over the same location.
**************************************************************************
##Tasks:
- Open images from within specified date range where criterias have previously been passed.
- Stack arrays of images and generate a cumulative sum layer of all layers.
- Upon script completion, information is printed including date range, number of days, available images and the maximum value acquired from the stacked array.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def array_to_img(in_img, output_img):
    cols = gdal.Open(in_img).RasterXSize
    rows = gdal.Open(in_img).RasterYSize
    proj = gdal.Open(in_img).GetProjection()
    geom = gdal.Open(in_img).GetGeoTransform()
    outDataset = gdal.GetDriverByName("GTiff").Create(output_img, cols, rows, 1, gdal.GDT_Float32)
    outDataset.SetProjection(proj)
    outDataset.SetGeoTransform(geom)
    outBand = outDataset.GetRasterBand(1)
    outBand.WriteArray(cumulative_array)
#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-s", "--time-start", required=True, help="Time lower bound (YYYY/MM/DD)")
        parser.add_argument("-t", "--tile", required=True, help="Filepath to MODIS tile of interest.")
        parser.add_argument("-d", "--days", type=int, help="Number of days to look at from input image.")
        parser.add_argument("-e", "--time-end", help="Time upper bound (YYYY/MM/DD)")
        parser.add_argument("-save", "--save", action="store_true", help="Include if you would like to save the image to current directory with outputfile named 'h**v**_YYYYMMDD_YYYYMMDD.tif'.")
        parser.add_argument("-o", "--output-img", help="Save image to specified files path with specified name.")

        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------
        # Ensure tile filepath ends correctly.
        if args.tile[-1] == "/":args.tile = args.tile[:-1]
        else:pass
        # Make sure either number of days or end date is specified. 
        if args.days == None and args.time_end == None:raise RuntimeError("Either the number of days ('-d') or an end date ('-time-end') must be specified.")
        # Make sure date format is all coorectly laid out. 
        if not (args.time_start[4] == "/" and args.time_start[-3] == "/"):raise RuntimeError("Start date format must be 'YYYY/MM/DD'.")
        if not args.time_end == None:
            if not (args.time_end[4] == "/" and args.time_end[-3] == "/"):raise RuntimeError("End date format must be 'YYYY/MM/DD'.")
        
        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        tile = args.tile
        sdate = datetime.strptime(os.path.join(args.time_start), "%Y/%m/%d").date()
        # End date based on specified or calculated from the start date from the number of days. 
        if args.time_end:edate = datetime.strptime(os.path.join(args.time_end), "%Y/%m/%d").date()
        elif args.days:edate = sdate + timedelta(days=args.days)
        # Acquire list of all the dates.
        date_list = [sdate + timedelta(days=x) for x in range((edate-sdate).days + 1)]
        # Full directory (including tile & date).
        directory_dates = ["/".join((tile, str(d.year), str('%02d' %d.month), str('%02d' %d.day))) for d in date_list if os.path.exists("/".join((tile, str(d.year), str('%02d' %d.month), str('%02d' %d.day))))]
        # List of directory and filenames.
        full_filepath = [files for dir in directory_dates for files in glob.glob(os.path.join(dir, "08*.tif"))]
        # Extract array of images within the date range that are available.
        array = np.array([np.array(Image.open(fimg)) for fimg in full_filepath])
        cumulative_array = np.sum(np.stack(array), axis=0)
        

        # If speecified, save image.
        if args.save == True and bool(args.output_img) == True:
            raise RuntimeError("Please include either '-save' or '-o' based on where you would like to save the output.")
        elif args.save == True:
            outfile = "_".join((os.path.split(tile)[1], "".join((str(sdate.year), str('%02d' %sdate.month), str('%02d' %sdate.day))), "".join((str(edate.year), str('%02d' %edate.month), str('%02d' %edate.day))))) + ".tif"
            array_to_img(full_filepath[0], outfile)  
        elif bool(args.output_img) == True:
            array_to_img(full_filepath[0], args.output_img)

        print("Based on specified criteria:")
        print(f"    Date range:                 {sdate} - {edate}")
        if bool(args.days) == True:
            print(f"    Range between dates:        {args.days} days")
        elif bool(args.time_end) == True:
            print(f"    Range between dates:        {(edate-sdate).days} days")
        print(f"    Number of available images: {len(full_filepath)}")
        print(f"    Maximum pixel overlap:      {int(np.amax(cumulative_array))}")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
