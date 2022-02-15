#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import gdal
import glob
import itertools
from operator import itemgetter
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
        # Extract all dates in datetime format.
        sdate = datetime.strptime(os.path.join(args.time_start), "%Y/%m/%d").date()
        edate = datetime.strptime(os.path.join(args.time_end), "%Y/%m/%d").date()
        dates = [sdate + timedelta(days=x) for x in range((edate - sdate).days + 1)]
       
        # Extract all tiles.
        tiles = set()
        if args.tile == None:
            all_tiles = [os.listdir(p) for p in products]
            for at in all_tiles:
                for t in at:
                    if t not in tiles:
                        tiles.add(t)
        elif not args.tile == None:
            for t in args.tile:
                tiles.add(t)
        print('date')
        print(dates)
        print('tiles')
        






        for p in products:
            tiles = [os.listdir(p)]
        print(tiles)
        sys.exit()

        # Check if tiles are specified, if not take all available tiles.
        if args.tile == None:
            filepath = [sorted(glob.glob(os.path.join(t+"/*"))) for t in products]
        elif not args.tile == None:
            filepath = [sorted(glob.glob(os.path.join(p+t))) for p in products for t in args.tile]
        

        # Need the length of the amount of lists to look at,
        # make sure those lists are the same lengths
        # the based on that, select element 0 - #. 
        # Rather than range could this be just len

        #length = (len(f) for f in filepath)
        length = []
        for f in filepath:
            lgth = len(f)
            length.append(lgth)
        if all(length):
            length = length[0]

        for f in filepath:
            print(f)
        sys.exit()
        
        for fp in filepath:
            for rng in range(len(fp)):
                print(fp[rng])
        sys.exit()





        print(list(map(itemgetter(0), filepath)))
        sys.exit()

        for f in filepath:
            print(itemgetter(f))
            sys.exit()





        for fp in filepath:
            print(fp)
            sys.exit()
            #print('--')
            for rang in range(len(fp)):
                #print(fp)
                print(rang)
                #sys.exit()






        count = 0
        for r in range(len(filepath)):
            for f in filepath:
                print(f[r])
            

            sys.exit()
        # Select files which fit in date range.
        sdate = datetime.strptime(os.path.join(args.time_start), "%Y/%m/%d").date()
        edate = datetime.strptime(os.path.join(args.time_end), "%Y/%m/%d").date()
        dates = [sdate + timedelta(days=x) for x in range((edate - sdate).days + 1)]
        
        # Full filepaths. (i.e. "product/tile/year/month/day")
        ffp = ["/".join((ftp, str(p[1].year), str('%02d' %p[1].month), str('%02d' %p[1].day))) for p in itertools.product(filepath,dates) for ftp in p[0] if os.path.isdir("/".join((ftp, str(p[1].year), str('%02d' %p[1].month), str('%02d' %p[1].day))))]

        


        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
