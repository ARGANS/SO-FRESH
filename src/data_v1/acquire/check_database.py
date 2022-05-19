#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import toolkit as tk
import os, sys
#sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+"/adjust")
#import toolkit_ as tk

parser = argparse.ArgumentParser(description="""
===================================================

===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

if __name__ == "__main__":
    try:
        parser.add_argument("-d", "--data-folder", help="Filepath to the data folder. i.e. 'fp/02_data'")
        parser.add_argument("-s", "--start-date", help="Lower-limit of date selection in YYYY/MM/DD format.")
        parser.add_argument("-e", "--end-date", help="Upper-limit of date selection in YYYY/MM/DD format.")
        parser.add_argument("-p", "--product", help="Product of interest for acquisition. Setup for: MYD09GA / MYDTBGA.")
        parser.add_argument("-aoi", "--area-of-interest", default=None, help="Specify area of interest, i.e. 'antarctica'.")
        parser.add_argument("-hmin", "--horizontal-minimum", default=None, help="Minimum horizontal MODIS sinusoidal tile.")
        parser.add_argument("-hmax", "--horizontal-maximum", default=None, help="Maximum horizontal MODIS sinusoidal tile.")
        parser.add_argument("-vmin", "--vertical-minimum", default=None, help="Minimum vertical MODIS sinusoidal tile.")
        parser.add_argument("-vmax", "--vertical-maximum", default=None, help="Maximum vertical MODIS sinusoidal tile.")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #==================================================================
        # Check for errors:
        
        #==================================================================
        # Execute acquisition:
        # hmin, hmax, vmin, vmax
        # args.horizontal_minimum, args.horizontal_maximum, args.vertical_minimum, args.vertical_maximum
        if (args.horizontal_minimum and args.horizontal_maximum and args.vertical_minimum and args.vertical_maximum) == None:
            if args.area_of_interest == "antarctica":
                hmin, hmax, vmin, vmax = 14, 24, 15, 17
        else: hmin, hmax, vmin, vmax = int(args.horizontal_minimum), int(args.horizontal_maximum), int(args.vertical_minimum), int(args.vertical_maximum)
        data_folder, sdate, edate, product = args.data_folder, args.start_date, args.end_date, args.product
        if product == "MYDTBGA":
                version = "006"
        elif product == "MYD09GA":
            version = input("What MYD09GA version would you like, 006 or 061?:\n")
        scdb = tk.scan_database(data_folder, sdate, edate, product, hmin, hmax, vmin, vmax)
        dates = scdb.dates_builder()
        tiles = scdb.tiles_builder()
        root = data_folder+"MODIS/"+product+"_"+version+"/"+"01_tiles/"
        if not os.path.isdir(root):
            raise RuntimeError(f"This directory does not exist:\n {root}")
        output = scdb.file_iterator(dates=dates,tiles=tiles,root=root)
        
       #==================================================================
        # Run and errors:
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
