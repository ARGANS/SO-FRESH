#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import toolkit as tk
import os, sys

parser = argparse.ArgumentParser(description="""
===================================================
### Automated Polynya Identification Tool ###
Execution options:
# data-folder
    ~ for data access & saving
# start-date
    ~lower limit of date selection
# end-date 
    ~upper limit of date selection
# products
    ~products for execution
    (ensure all data is downloaded, pre-processed & (optional) fused)
===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

if __name__ == "__main__":
    try:
        parser.add_argument("-d", "--data-folder", help="Filepath to the folder containing all data.")
        parser.add_argument("-s", "--start-date", help="Time lower bound (YYYY/MM/DD).")
        parser.add_argument("-e", "--end-date", help="Time upper bound (YYYY/MM/DD).")
        parser.add_argument("-p", "--products", nargs="+", help="Product(s) to execute APIT on.")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #==================================================================
        # Check for errors:

        #==================================================================
        # Execute the Automated Polynya Identification Tool:
        data_folder, sdate, edate, products=args.data_folder, args.start_date, args.end_date, args.products
        data = tk.data_selector(data_folder, sdate, edate, products).collector()
        run = tk.execute_APIT(data, products).classification(data_folder)
        tk.execute_APIT(data, products).generate_netcdf(run, data_folder, sdate, edate)
        #==================================================================
        # Run and errors:
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
