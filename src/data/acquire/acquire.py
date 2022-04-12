#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import toolkit as tk
#import os, sys
#sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+"/adjust")
#import toolkit_ as tk

parser = argparse.ArgumentParser(description="""
===================================================

===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

if __name__ == "__main__":
    try:
        parser.add_argument("-d", "--data-folder")
        parser.add_argument("-s", "--start-date")
        parser.add_argument("-e", "--end-date")
        parser.add_argument("-p", "--product")
        parser.add_argument("-hmin", "--horizontal-minimum")
        parser.add_argument("-hmax", "--horizontal-maximum")
        parser.add_argument("-vmin", "--vertical-minimum")
        parser.add_argument("-vmax", "--vertical-maximum")
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
