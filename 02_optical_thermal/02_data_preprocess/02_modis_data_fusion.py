#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import gdal
import glob
import sys, os
from datetime import datetime, timedelta

#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# 
**************************************************************************
##Tasks:
- 
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------


#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-dir", required=True, help="Input directory containing MODIS products.").completer = FilesCompleter(allowednames=(".tif"))
        parser.add_argument("-t", "--tile", nargs="+", help="Opportunity to specify tiles of interest.")
        parser.add_argument("-s", "--time-start", required=True, help="Time lower bound (YYYY/MM/DD).")
        parser.add_argument("-e", "--time-end", required=True, help="Time upper bound (YYYY/MM/DD).")
        parser.add_argument("-p", "--products", required=True, nargs="+", help="Products to fuse.")
        parser.add_argument("-v", "--version", required=True, nargs="+", help="Version of the product to fuse. Please specify for each input product")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------
        # Check for matching length of products and version.
        if not len(args.products) == len(args.version):
            raise RuntimeError("Number of products and version must be equal.")
        # Check that the full filepaths exist. 
        products = [(args.input_dir+p+"_"+v+"/") for p, v in zip(args.products, args.version)]
        if not all([os.path.isdir(p) for p in products]):
            new_line = '\n'
            raise RuntimeError(f"There is an issue with the product name and/or version, please check spelling. Generated filepaths include the following: {new_line}{new_line.join(map(str,products))}")

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        # Check if tiles are specified, if not take all available tiles.
        #print(args.input_dir)
        #print(args.tile)
        
        if args.tile == None:
            filepath = [t for t in sorted(glob.glob(os.path.join(args.input_dir, "*")))]
            # rather than the *, could I put the 'products' variable, to specify and then have the * to select all tiles and thn I can play with the date stuff after it.
            print(filepath)
            sys.exit()
        elif not args.tile == None:
            filepath = [os.path.join(args.filepath_tile, t) for t in args.tile if os.path.isdir(os.path.join(args.input_dir, t)) == True]




        print(args.input_dir)
        sdate = datetime.strptime(os.path.join(args.time_start), "%Y/%m/%d").date()
        edate = datetime.strptime(os.path.join(args.time_end), "%Y/%m/%d").date()
        dates = [sdate + timedelta(days=x) for x in range((edate - sdate).days + 1)]
        
        ffp = ["/".join((p[0], str(p[1].year), str('%02d' %p[1].month), str('%02d' %p[1].day))) for p in itertools.product(filepath,dates) if os.path.isdir("/".join((p[0], str(p[1].year), str('%02d' %p[1].month), str('%02d' %p[1].day))))]
        print(sdate)
        print(edate)
        print('---')
        print(dates)

        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
