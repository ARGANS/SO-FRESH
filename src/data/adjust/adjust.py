#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
# @James Hickson | Argans UK | jhickson@argans.co.uk

import argparse, argcomplete
from argcomplete.completers import ChoicesCompleter, FilesCompleter
from pathlib import Path
import toolkit_ as tk
from tqdm import tqdm
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+"/acquire")
import toolkit as tk_db

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
            if product == "MYD09GA":
                version = input("What MYD09GA version would you like, 006 or 061?:\n")
            images, dates, missing= tk.MODIS(data_folder, product).build_filepath(sdate, edate)

            '''
            if not len(missing) == 0:
                # Data inconsitencies means the following two tiles require removal from the "missing" list.
                for m in missing[:]:
                    if (os.path.split(m[0])[1]) == "h14v17" or (os.path.split(m[0])[1]) == "h21v17":
                        missing.remove(m)
            '''
            if not len(missing) == 0:
                # Offer the opportunity to download the missing data. 
                dwnld_missing = input("Based on entries, there is missing data. Would you like to see this list? (Y or N)\n")
                if dwnld_missing == "Y":
                    print("The following dates have missing data:")
                    print(*missing,sep='\n')
                    dwnld = input("Would you like to download the missing data? (Y or N)\n")
                elif dwnld_missing == "N":
                    dwnld = input("Would you like to download the missing data? (Y or N)\n")
                if dwnld == "Y":
                    mdate = list(set([m[1]+"/" for m in missing]))
                    mtile = list(set([os.path.split(m[0])[1] for m in missing]))
                    mroot = list(set([os.path.split(m[0])[0]+"/" for m in missing]))[0]
                    output = tk_db.scan_database(data_folder, None, None, product, None, None, None, None).file_iterator(dates=mdate,tiles=mtile,root=mroot)
                    if not output == None:
                        for imgs in images:
                            if all(("/".join(os.path.split(output)[0].rsplit("/", 3)[1:])) in i for i in imgs):
                                imgs.append(output)
                if dwnld == "N":
                    # Remove the date to which data is missing to avoid creating mosaic.
                    mdates = list(set([m[1] for m in missing]))
                    if any(x in dates for x in mdates):
                        for md in mdates:
                            for i, d in zip(images, dates):                                
                                if all(md in im for im in i) and md in d:
                                    images.remove(i)
                                    dates.remove(d)

            for i, d in zip(images, dates): 
                if product == "MYDTBGA": version="006"
                mosaic_img = tk.MODIS(data_folder, product).build_mosaic(i, d, version, args.resample)
                if args.resample == True:
                    tk.MODIS(data_folder, product).resample(mosaic_img)
                for img in i:
                    os.remove(img)
        #==================================================================
        # Run and errors:
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
