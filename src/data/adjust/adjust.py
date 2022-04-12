#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import toolkit as tk
from tqdm import tqdm
import os, sys

parser = argparse.ArgumentParser(description="""
===================================================
### Imagery pre-processing and adjustment. ###
Pre-processing and adjustment options:
- File renaming (add 01 at the begining of the file)
    - input: raw downloaded file.
- Project MYD09GA (optical).
    - input: 01_MYD09GA.jpg file.
- Extract and normalise MYDTBGA (thermal).
    - input: 01_MYDTBGA.hdf file.
- Extract and project AMSR2 SIC (Sea Ice Concentration).
    - input: 01_ESACCI_AMSR2_SIC.nc file.
    - include '-r' for data to be resampled*
===================================================
### Imagery mosaic and resampling. ###
For AMSR2 Ice Concentration - check step above.
For MYD09GA & MYDTBGA:
    - data_folder: filepath to folder containing all data in.
    - start/end_date: lower and upper time threshold (each date will be an individual image).
    - product: MYD09GA & MYDTBGA
    - include '-r' for the data to be resampled.
===================================================""",
formatter_class=argparse.RawDescriptionHelpFormatter)

if __name__ == "__main__":
    try:
        parser.add_argument("-i", "--input", required=False, nargs="+", help="Provide filepath to the source of imagery you would like to process i.e. -i file/path/*/*/*.jpg")
        parser.add_argument("-d", "--data-folder", required=False, help="Filepath to the folder containing all data.")
        parser.add_argument("-s", "--start-date", required=False, help="Time lower bound (YYYY/MM/DD).")
        parser.add_argument("-e", "--end-date", required=False, help="Time upper bound (YYYY/MM/DD).")
        parser.add_argument("-p", "--product", required=False, help="Product to merge [MYD09GA, MYDTBGA].")
        parser.add_argument("-r", "--resample", required=False, action="store_true", help="For data to be resampled too, include '-r' in the command. This will currently be done to 0.1 (x & y) degrees resolution.")
        #parser.add_argument()
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #==================================================================
        # Check for errors:
        #if not all([args.data_folder, args.start_date, args.end_date, args.product]):
            #raise RuntimeError("For mosaic creationg required inputs are missing.")
        #if args.product=="MYD09GA" or args.product=="MYDTBGA":pass
        #else: raise RuntimeError("This product does not exist, consider using either MYD09GA or MYDTBGA")
        #==================================================================
        # Execute adjustments:
        if bool(args.input) == True:
            for img in tqdm(args.input):
                filename = os.path.basename(img)
                # Check if file requires to be renamed.
                if filename.startswith(("BROWSE", "MYDTBGA", "ESACCI")) and filename.endswith((".jpg", ".hdf", ".nc")):
                    walk_fp = os.path.split(img)[0]
                    tk.aux_func(walk_fp).rename()
                # Check if thermal data requires to be extracted from hdf.
                elif filename.startswith("01_") and filename.endswith(".hdf"):
                    extract = tk.MYDTBGA_preprocess(img).extract_MYDTBGA()
                    thermal_bnd = (extract+"_4.tif")
                    tk.MYDTBGA_preprocess(thermal_bnd).normalise()
                # Check if optical data requires to be projected.
                elif filename.startswith("01_") and filename.endswith(".jpg"):
                    tk.MYD09GA_preprocess(img).assign_geometry()
                # Check if SIC requires to be extracted from netCDF.
                elif filename.startswith("01_") and filename.endswith(".nc"):
                    sic_img = tk.amsr2_preprocess(img).extract_from_netcdf()
                    rpjct_sic_img = tk.amsr2_preprocess(sic_img).reproject(args.resample)
                    if args.resample == True:
                        resampled_sic_img = tk.amsr2_preprocess(rpjct_sic_img).resample()
                        tk.amsr2_preprocess(resampled_sic_img).mask()
        elif bool(args.data_folder)==True and bool(args.start_date)==True and bool(args.end_date)==True and bool(args.product)==True:
            data_folder, sdate, edate, product=args.data_folder, args.start_date, args.end_date, args.product
            images, dates= tk.MODIS(data_folder, product).build_filepath(sdate, edate)
            for i, d in zip(images, dates):
                mosaic_img = tk.MODIS(data_folder, product).build_mosaic(i, d, args.resample)
                if args.resample == True:
                    tk.MODIS(data_folder, product).resample(mosaic_img)
        #==================================================================
        # Run and errors:
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
