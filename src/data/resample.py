#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import toolkit as tk
import os, sys

parser = argparse.ArgumentParser(description="""
===================================================
# Resample imagery to a given resolution.
===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

if __name__ == "__main__":
    try:
        parser.add_argument("-i", "--input", required=True, nargs="+", help="Provide filepath to the source of imagery you would like to process i.e. -i file/path/*/*/*.jpg")
        #parser.add_argument()
        #parser.add_argument()
        #parser.add_argument()
        #parser.add_argument()
        #parser.add_argument()
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #==================================================================
        # Check for errors:



        #==================================================================
        # Execute adjustments:
        for img in args.input:
            filename = os.path.basename(img)
        #==================================================================
        # Run and errors:
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
