#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import toolkit as tk
import os, sys

parser = argparse.ArgumentParser(description="""
===================================================
### Imagery fusion ###
Fusing of listed products - mosaic MODIS imagery before.

# data-folder
    ~for data access & saving
# start-date
    ~lower limit of date selection
# end-date 
    ~upper limit of date selection
# products
    ~list of products for fusing
===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

if __name__ == "__main__":
    try:
        parser.add_argument("-d", "--data-folder", help="Filepath to the folder containing all data.)
        parser.add_argument("-s", "--start-date", help="Time lower bound (YYYY/MM/DD).")
        parser.add_argument("-e", "--end-date", help="Time upper bound (YYYY/MM/DD).")
        parser.add_argument("-p", "--products", nargs="+", help="Product(s) to execute APIT on.")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #==================================================================
        # Check for errors:
        
        #==================================================================
        # Execute fusion:
        tk.fusion(args.data_folder, args.start_date, args.end_date, args.products)
        #==================================================================
        # Run and errors:
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
