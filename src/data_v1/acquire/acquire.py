#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import toolkit as tk
#import os, sys # If required - import the adjust toolkit.
#sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+"/adjust")
#import toolkit_ as tk

parser = argparse.ArgumentParser(description="""
===================================================
When an archive or series of data is required, provide the following inputs for data download.
##################################
# data-folder
    ~for data saving filepath building
# start-date
    ~lower limit of date selection
# end-date 
    ~upper limit of date selection
# product
    ~product for download
# horizontal-minimum
# horizontal-maximum
# vertical-minimum
# vertical-maximum
    ~aoi specification
##################################
===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

if __name__ == "__main__":
    try:
        parser.add_argument("-d", "--data-folder", help="Filepath to the data folder. i.e. 'fp/02_data'")
        parser.add_argument("-s", "--start-date", help="Lower-limit of date selection in YYYY/MM/DD format.")
        parser.add_argument("-e", "--end-date", help="Upper-limit of date selection in YYYY/MM/DD format.")
        parser.add_argument("-p", "--product", help="Product of interest for acquisition. Setup for: MYD09GA / MYDTBGA.")
        parser.add_argument("-hmin", "--horizontal-minimum", help="Minimum horizontal MODIS sinusoidal tile.")
        parser.add_argument("-hmax", "--horizontal-maximum", help="Maximum horizontal MODIS sinusoidal tile.")
        parser.add_argument("-vmin", "--vertical-minimum", help="Minimum vertical MODIS sinusoidal tile.")
        parser.add_argument("-vmax", "--vertical-maximum", help="Maximum vertical MODIS sinusoidal tile.")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #==================================================================
        # Check for errors:
        
        #==================================================================
        # Execute acquisition:
        data_folder, sdate, edate, product, hmin, hmax, vmin, vmax = args.data_folder, args.start_date, args.end_date, args.product, args.horizontal_minimum, args.horizontal_maximum, args.vertical_minimum, args.vertical_maximum
        links, outdirs = tk.lpdaac_download(data_folder, sdate, edate, product, hmin, hmax, vmin, vmax).url_extractor()
        tk.lpdaac_download(data_folder, sdate, edate, product, hmin, hmax, vmin, vmax).download(links, outdirs)
        #==================================================================
        # Run and errors:
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
