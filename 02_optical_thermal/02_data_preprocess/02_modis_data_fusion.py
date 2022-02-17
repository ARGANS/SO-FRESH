#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import gdal
import glob
import itertools
from operator import itemgetter
import sys, os
from datetime import datetime, timedelta
import pprint

#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# 
**************************************************************************
##Tasks:
- 
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def imagery_fusion(imgs_fp, outdir, products, versions):
    # Fuse imagery from list. 
    # Output name
    files = [sorted(glob.glob((i+"/02_*.tif"))) for i in imgs_fp]
    files = [item for sublist in files for item in sublist]
    if all(len(list(files)) == len(l) for l in list((products, versions))):
        prod_ver = "_".join([(p+"_"+v) for p, v in zip(args.products, args.version)]) # Product & version
        fp = []
        for i in imgs_fp:
            day, month, year, tile = os.path.split(i)[1], os.path.split(os.path.split(i)[0])[1], os.path.split(os.path.split(os.path.split(i)[0])[0])[1], os.path.split(os.path.split(os.path.split(os.path.split(i)[0])[0])[0])[1]
            fp_indv = "/".join((prod_ver, tile, year, month, day))
            fp.append(fp_indv)
        if all(fp): fp = fp[0]
    else:
        raiseRuntimeError("The number of original inputs (products and versions) do not match the number of selected images.")

    output_dir = (outdir+prod_ver+"/"+tile+"/"+year+"/"+month+"/"+day+"/")
    output_name = "_".join(("02",prod_ver,tile,(year+month+day+".tif")))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    cmd = "gdal_merge.py -q -of GTIFF -seperate -o %s %s"%((output_dir+output_name), " ".join(files))
    os.system(cmd)
    return((output_dir+output_name))

#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-dir", required=True, help="Input directory containing MODIS products.").completer = FilesCompleter(allowednames=(".tif"))
        parser.add_argument("-t", "--tile", nargs="+", help="Opportunity to specify tiles of interest.")
        parser.add_argument("-s", "--time-start", required=True, help="Time lower bound (YYYY/MM/DD).")
        parser.add_argument("-e", "--time-end", required=True, help="Time upper bound (YYYY/MM/DD).")
        parser.add_argument("-p", "--products", required=True, nargs="+", help="Products to fuse.")
        parser.add_argument("-v", "--version", required=True, nargs="+", help="Version of the product to fuse. Please specify for each input product")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------
        # Check for matching length of products and version.
        if not len(args.products) == len(args.version):
            raise RuntimeError("Number of products and version must be equal.")
        # Check that the full filepaths exist. 
        products = [(args.input_dir+p+"_"+v+"/") for p, v in zip(args.products, args.version)]
        if not all([os.path.isdir(p) for p in products]):
            new_line = '\n'
            raise RuntimeError(f"There is an issue with the product name and/or version, please check spelling. Generated filepaths include the following: {new_line}{new_line.join(map(str,products))}")

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        # Extract all dates in datetime format.
        sdate, edate  = datetime.strptime(os.path.join(args.time_start), "%Y/%m/%d").date(), datetime.strptime(os.path.join(args.time_end), "%Y/%m/%d").date()
        dates = [sdate + timedelta(days=x) for x in range((edate - sdate).days + 1)]
        # Check if tiles are specified, if not take all available tiles.
        if args.tile == None:
            tile_filepath = [sorted(glob.glob(os.path.join(t+"/*"))) for t in products]
        elif not args.tile == None:
            tile_filepath = [sorted(glob.glob(os.path.join(p+t))) for p in products for t in args.tile]
        # Generate filepaths with tile & dates. 
        filepath = []
        for sub_f in tile_filepath:
            for f in itertools.product(sub_f, dates):
                filepath.append("/".join((f[0], str(f[1].year), str('%02d' %f[1].month), str('%02d' %f[1].day))))
        # Create dictionary to store product filepaths in.
        filepath_dict = {key:[] for key in args.products}
        for f in filepath:
            match = [prods for prods in args.products if (prods in f)]
            if len(match) > 1: raiseRuntimeError("Multiple matches were found in string, please name filepaths more appropriately.")
            filepath_dict[match[0]].append(f) 
        
        keys, values = list(filepath_dict.keys()), list(filepath_dict.values())
        for i in range(len(filepath_dict[keys[0]])):
            items = [lst[i] for lst in values]
            imagery_fusion(items, args.input_dir, args.products, args.version)
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
