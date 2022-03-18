#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import toolkit as tk
import os, sys

parser = argparse.ArgumentParser(description="""
===================================================
# Imagery pre-processing and adjustment.
Pre-processing and adjustment options:
- File renaming
    - input: Add 01 at the begining of the raw downloaded files.
- Project MYD09GA (optical).
    - input: 01_MYD09GA.jpg file.
- Extract and normalise MYDTBGA (thermal).
    - input: 01_MYDTBGA.hdf file.
- Extract and project AMSR2 SIC (Sea Ice Concentration).
    - input: 01_
===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

if __name__ == "__main__":
    try:
        parser.add_argument("-i", "--input", required=False, nargs="+", help="Provide filepath to the source of imagery you would like to process i.e. -i file/path/*/*/*.jpg")
        parser.add_argument("-s", "--start-date", required=False, help"")
        parser.add_argument("-e", "--end-date", required=False, help"")
        parser.add_argument("-p", "--product", required=False, help"")
        #parser.add_argument()
        #parser.add_argument()
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #==================================================================
        # Check for errors:



        #==================================================================
        # Execute adjustments:
        if args.input == True:
            for img in args.input:
                filename = os.path.basename(img)
                # Check if file requires to be renamed.
                if filename.startswith(("BROWSE", "MYDTBGA", "ESACCI")) and filename.endswith((".jpg", ".hdf", ".nc")):
                    walk_fp = os.path.split(img)[0]
                    tk.modis_rename(walk_fp)
                # Check if thermal data requires to be extracted from hdf.
                elif filename.startswith("01_") and filename.endswith(".hdf"):
                    extract = tk.modis_preprocess(img).extract_MYDTBGA()
                    thermal_bnd = (extract+"_4.tif")
                    tk.modis_preprocess(thermal_bnd).normalise()
                # Check if optical data requires to be projected.
                elif filename.startswith("01_") and filename.endswith(".jpg"):
                    tk.modis_preprocess(img).assign_geometry()
                # Check if SIC requires to be extracted from netCDF.
                elif filename.startswith("01_") and filename.endswith(".nc"):
                    sic_img = tk.amsr2_preprocess(img).extract_from_netcdf()
                    tk.amsr2_preprocess(sic_img).reproject()
        elif args.start_date==True and args.end_date==True and args.product==True:
            # set up mosaic building
            #   -dates & products
            # resample
            print('here')
        #==================================================================
        # Run and errors:
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
