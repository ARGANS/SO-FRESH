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
===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

if __name__ == "__main__":
    try:
        parser.add_argument("-d", "--data-folder")
        parser.add_argument("-s", "--start-date")
        parser.add_argument("-e", "--end-date")
        parser.add_argument("-p", "--products", nargs="+")
        #parser.add_argument()
        #parser.add_argument()
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
