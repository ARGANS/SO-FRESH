#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
from datetime import datetime, timedelta
import numpy as np
import itertools
from PIL import Image
import rasterio
import gdal
import glob, sys, os
import netCDF4 as nc

#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Produce a netCDF of results from a specific date range, displaying occurence 
**************************************************************************
##Tasks:
- Select all data based on input specifications (i.e. date range, tiles of interest (if not all))
- Stack arrays of those from the same tiles and turn in to one image for each tile. Extract co-ordinates for each tile and make sure they are all the same.
- Cumulative sum all stacked arrays in to raster files.
- Merge raster files and convert in to a sigular NetCDF file.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def array_to_img(array, coords):
    cols, rows, proj, geom = coords[0], coords[1], coords[2], coords[3]
    outDataset = gdal.GetDriverByName("GTiff").Create(out_file, cols, rows, 1, gdal.GDT_Float32)
    outDataset.SetProjection(proj)
    outDataset.SetGeoTransform(geom)
    outBand = outDataset.GetRasterBand(1)
    outBand.WriteArray(array)

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-f", "--filepath-tile", required=True, help="Filepath to the folder containing all tiles.")
        parser.add_argument("-t", "--tile-folder", nargs="+", help="Opportunity to specify tiles of interest.")
        parser.add_argument("-s", "--time-start", required=True, help="Time lower bound (YYYY/MM/DD)")
        parser.add_argument("-e", "--time-end", required=True, help="Time upper bound (YYYY/MM/DD)")
        parser.add_argument("-o", "--output-dir", required=True, help="Output directory.")
        parser.add_argument("-r", "--remove", required=False, action="store_true", help="Include to remove created 'tif' files")
        parser.add_argument("-c", "--clip", required=False, action="store_true", help="Time upper bound (YYYY/MM/DD)")
        parser.add_argument("-ymin", "--y-min", required=False, type=int, help="Minimum latitude to clip the data, if left empty, this will take the data extent.")
        parser.add_argument("-ymax", "--y-max", required=False, type=int, help="Maximum latitude to clip the data, if left empty, this will take the data extent.")
        parser.add_argument("-xmin", "--x-min", required=False, type=int, help="Minimum longitude to clip the data, if left empty, this will take the data extent.")
        parser.add_argument("-xmax", "--x-max", required=False, type=int, help="Maximum longitude to clip the data, if left empty, this will take the data extent.")
        parser.add_argument("-z", "--zero-nan", required=False, action="store_true", help="Set all values equal to 0 to NaN.")
        parser.add_argument("-n", "--nan", required=False, help="Set values equal to or more than this value to zero")

        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------
        if args.clip == True:
            if args.y_min == None and args.y_max == None and args.x_min == None and args.x_max == None:
                raise RuntimeError("If you wish to clip the data to a specific extent, please specify at least one extent you would like to clip to.")

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        # Check if tiles are specified, if not take all available tiles.
        if args.tile_folder == None:
            filepath = [t for t in sorted(glob.glob(os.path.join(args.filepath_tile, "*")))]
        elif not args.tile_folder == None:
            filepath = [os.path.join(args.filepath_tile, t) for t in args.tile_folder if os.path.isdir(os.path.join(args.filepath_tile, t)) == True]

        # Add date to the filepath.
        print("Collecting data based on entries...")
        sdate = datetime.strptime(os.path.join(args.time_start), "%Y/%m/%d").date()
        edate = datetime.strptime(os.path.join(args.time_end), "%Y/%m/%d").date()
        dates = [sdate + timedelta(days=x) for x in range((edate - sdate).days + 1)]
        # Full filepath with tile/year
        ffp = ["/".join((p[0], str(p[1].year), str('%02d' %p[1].month), str('%02d' %p[1].day))) for p in itertools.product(filepath,dates) if os.path.isdir("/".join((p[0], str(p[1].year), str('%02d' %p[1].month), str('%02d' %p[1].day))))]

        # Create nested lists for each tile.
        ffp_split = []
        for tile in filepath:
            ffp_by_tile = []
            for full in ffp:
                if tile in full:
                    ffp_by_tile.append(full)
                else:
                    pass
            if ffp_by_tile:
                ffp_split.append(ffp_by_tile)
        # Add image to the filepath in a dictionary where each key is a tile.
        img_dict = {}
        for fp in ffp_split:
            img_by_tile = []
            for f in fp:
                tile = [t for t in (f.split(os.sep)) if "h" and "v" in t][0]
                if not tile in img_dict:
                    img_dict[tile] = []
                # Search for images and add them to the appropriate tile key.
                for i in glob.glob(os.path.join(f, "03f*.tif")):
                    if i:
                        img_dict[tile].append(i)
        # Remove key which do not any associated values.
        img_dict = {k : v for k, v in img_dict.items() if v}
        img_array = {}
        img_coords = {}
        for key in img_dict.keys():
            # Add keys to the new dictionaries.
            if not key in img_array:
                img_array[key] = []
            if not key in img_coords:
                img_coords[key] = []
            # From images extract them to arrays and their co-ordinates.
            for img in img_dict[key]:
                data = gdal.Open(img)
                img_array[key].append(data.ReadAsArray())
                proj_info = data.RasterXSize, data.RasterYSize, data.GetProjection(), data.GetGeoTransform()
                img_coords[key].append(proj_info)
                del data
        # If all co-ordinates match in the list, select the first one.
        for key in img_coords.keys():
            if all(img_coords[key]):
                img_coords[key] = [img_coords[key][0]]
            else:
                raiseRuntimeError(f"The key ({key}) co-ordinates do not all match")
        print("Processing...")
        # Stack all available arrays and calculate their cumulative sum.
        for key in img_array.keys():
            sum_arr = np.sum(np.stack(img_array[key]), axis=0)
            img_array[key] = sum_arr

        # Pull product name for file names. If statement incase the filepath ends with a '/'.
        if bool(args.filepath_tile.rsplit("/",2)[-1]) == True:product = args.filepath_tile.rsplit("/",2)[-1]
        else:product = args.filepath_tile.rsplit("/", 2)[-2]

        # Export stacked tile arrays to tifs files.
        print("Exporting images...")
        imgs_4_netcdf = []
        for key in img_array.keys():
            output_fp = args.output_dir + (args.time_start.replace("/", "") + "_" + args.time_end.replace("/", "")+"/")
            out_file = os.path.join(output_fp + "tmp/" + "SOFRESH_MODIS_POLYNYA_" + product + "_" + args.time_start.replace("/", "") + "_" + args.time_end.replace("/", "") + "_" + key + ".tif")
            if not os.path.exists(output_fp):
                os.mkdir(output_fp)
            if not os.path.exists(os.path.join(output_fp + "tmp/")):
                os.mkdir(os.path.join(output_fp + "tmp/"))
            array_to_img(img_array[key], img_coords[key][0])
            imgs_4_netcdf.append(out_file)

        if not os.path.exists(os.path.join(output_fp + "00_netcdf/")):
            os.mkdir(os.path.join(output_fp + "00_netcdf/"))

        # Mosaic list of tif tiles to one netCDF file.
        out_tif= os.path.join(output_fp + "tmp/" + "SOFRESH_MODIS_POLYNYA_" + product + "_" + args.time_start.replace("/", "") + "_" + args.time_end.replace("/", "") + ".tif")
        if args.clip == True:clip_tif = os.path.join(output_fp + "tmp/" + "SOFRESH_MODIS_POLYNYA_" + product + "_" + args.time_start.replace("/", "") + "_" + args.time_end.replace("/", "") + "_clipped.tif")
        out_nc = os.path.join(output_fp + "00_netcdf/" + "SOFRESH_MODIS_POLYNYA_" + product + "_" + args.time_start.replace("/", "") + "_" + args.time_end.replace("/", "") + ".nc")
        # Merge all tif files in to one.
        print("Merging images...")
        os.system("gdal_merge.py -q -of GTIFF -o %s %s"%(out_tif, ' '.join(imgs_4_netcdf)))
        # If clip is specified, clip the tif before converting to netCDF file.
        if args.clip == True:
            img = gdal.Open(out_tif)
            array = img.ReadAsArray()
            ypixels, xpixels = array.shape[0], array.shape[1]
            # Extract image geometry
            gt = img.GetGeoTransform()
            xmin, xres = gt[0], gt[1]
            ymax, yres = gt[3], gt[5]
            xmax = (xres*xpixels)+xmin
            ymin = (yres*ypixels)+ymax
            # Check if x or y limits have been set.
            if bool(args.y_min) == True:ymin=args.y_min
            if bool(args.y_max) == True:ymax=args.y_max
            if bool(args.x_min) == True:xmin=args.x_min
            if bool(args.x_max) == True:xmax=args.x_max
            latbounds = [ymin, ymax]
            lonbounds = [xmin, xmax]
            ulx, uly, lrx, lry = lonbounds[0], latbounds[1], lonbounds[1], latbounds[0]
            os.system("gdal_translate -q -of GTiff -projwin %s %s %s %s %s %s"%(ulx, uly, lrx, lry, out_tif, clip_tif))
            # If required, set values to NaN.
            if args.zero_nan == True or args.nan == True:
                in_tif = clip_tif
                with rasterio.open(in_tif) as src:
                    array=src.read(1)
                    profile = src.profile
                if args.zero_nan == True:array[array==0] = np.nan
                if args.nan == True:array[array>=args.nan] = np.nan
                with rasterio.open(clip_tif, "w", **profile) as output:
                    output.write(array, 1)
            out_tif = clip_tif
        # Convert tif to netCDF.
        os.system("gdal_translate -q -of netCDF %s %s"%(out_tif, out_nc))
        print("NetCDF file successfully created.")
        # Change variable long name in netCDF file to "Count Frequency".
        os.system("ncatted -O -a long_name,Band1,o,c,Count_Frequency %s"%(out_nc))

        # Remove all '.tif' files created if specified.
        [os.remove(i) for i in imgs_4_netcdf if args.remove == True]
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
