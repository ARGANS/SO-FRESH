#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
from datetime import datetime, timedelta
import numpy as np
import itertools
from PIL import Image
import pprint, gdal
import glob, sys, os

#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Produce visuals based on results.
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
def extract_coords(product):
    ds = gdal.Open(product)
    xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
    xmax = xmin+xpixel*ds.RasterXSize
    ymin = ymax+ypixel*ds.RasterYSize
    return xmin, xmax, xpixel, ymin, ymax, ypixel

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-f", "--filepath-tile", required=True, help="Filepath to the folder containing all tiles.")
        parser.add_argument("-t", "--tile-folder", nargs="+", help="Opportunity to specify tiles of interest.")
        parser.add_argument("-s", "--time-start", required=True, help="Time lower bound (YYYY/MM/DD)")
        parser.add_argument("-e", "--time-end", required=True, help="Time upper bound (YYYY/MM/DD)")
        #parser.add_argument("-d", "--days", type=int, help="Time upper bound (YYYY/MM/DD)")

        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        # Check if tiles are specified, if not take all available tiles.
        if args.tile_folder == None:
            filepath = [t for t in sorted(glob.glob(os.path.join(args.filepath_tile, "*")))]
        elif not args.tile_folder == None:
            filepath = [os.path.join(args.filepath_tile, t) for t in args.tile_folder if os.path.isdir(os.path.join(args.filepath_tile, t)) == True]
        
        # Add date to the filepath.
        sdate = datetime.strptime(os.path.join(args.time_start), "%Y/%m/%d").date()
        edate = datetime.strptime(os.path.join(args.time_end), "%Y/%m/%d").date()
        dates = [sdate + timedelta(days=x) for x in range((edate - sdate).days + 1)]
        # Full filepath with tile/year
        ffp = ["/".join((p[0], str(p[1].year), str('%02d' %p[1].month), str('%02d' %p[1].day))) for p in itertools.product(filepath,dates) if os.path.isdir("/".join((p[0], str(p[1].year), str('%02d' %p[1].month), str('%02d' %p[1].day))))]
        
        # Create nested lists for each tile.
        ffp_split = []
        for tile in filepath:
            ffp_by_tile = []
            for full in ffp:
                if tile in full:
                    ffp_by_tile.append(full)
                else:
                    pass
            if ffp_by_tile:
                ffp_split.append(ffp_by_tile)
        imagery = []
        for fp in ffp_split:
            img_by_tile = []
            for f in fp:
                for i in glob.glob(os.path.join(f, "08*.tif")):
                    if i:
                        imagery.append(i)
            #if img_by_tile:
                #imagery.append(img_by_tile)
        #pprint.pprint(imagery)
        #sys.exit()
        img_coords = []
        for img in imagery:
            (xmin, xmax, xpixel, ymin, ymax, ypixel) = extract_coords(img)
            img_coords.append((xmin, xmax, xpixel, ymin, ymax, ypixel))
            
        #sys.exit()
        array = np.array([np.array(Image.open(fimg)) for fimg in imagery])

        print(len(array))
        print('--------')
        print(img_coords)
        print(len(img_coords))
        sys.exit()
        cumu_array = np.sum(np.stack(array), axis=0)

        import matplotlib.pyplot as plt
        plt.imshow(cumu_array)
        plt.colorbar()
        plt.show()
        #print(cumu_array)
        sys.exit()



        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
