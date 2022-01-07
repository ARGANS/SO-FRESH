#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
from datetime import datetime, timedelta
import numpy as np
import itertools
from PIL import Image
import rasterio
import xarray as xr
import pprint, gdal
import glob, sys, os

#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Produce visuals based on results.
**************************************************************************
##Tasks:
-
-
-
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def extract_coords(product):
    ds = gdal.Open(product)
    xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
    xmax = xmin+xpixel*ds.RasterXSize
    ymin = ymax+ypixel*ds.RasterYSize
    return [xmin, xmax, xpixel, ymin, ymax, ypixel]

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
        #parser.add_argument("-d", "--days", type=int, help="Time upper bound (YYYY/MM/DD)")

        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        # Check if tiles are specified, if not take all available tiles.
        if args.tile_folder == None:
            filepath = [t for t in sorted(glob.glob(os.path.join(args.filepath_tile, "*")))]
        elif not args.tile_folder == None:
            filepath = [os.path.join(args.filepath_tile, t) for t in args.tile_folder if os.path.isdir(os.path.join(args.filepath_tile, t)) == True]

        # Add date to the filepath.
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

        img_dict = {}
        for fp in ffp_split:
            img_by_tile = []
            for f in fp:
                tile = [t for t in (f.split(os.sep)) if "h" and "v" in t][0]
                if not tile in img_dict:
                    img_dict[tile] = []
                # Search for images and add them to the appropriate tile key.
                for i in glob.glob(os.path.join(f, "08*.tif")):
                    if i:
                        img_dict[tile].append(i)
        # Remove key which do not any associated values.
        img_dict = {k : v for k, v in img_dict.items() if v}

        #pprint.pprint(img_dict)
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
                '''
                ############################################################
                #import rioxarray as rxr
                #data = rxr.open_rasterio(img)
                #data = xr.open_dataset(img, engine="rasterio")
                ### ARRAY ###
                #img_array[key].append(data.data)
                #print(data.data)
                ### COORDS ###
                #img_coords[key].append(data.coords)
                ############################################################
                #XARRAY ONLY
                #dA = xr.DataArray(data.data, dims=['band', 'x', 'y'], coords={'band': data.coords["band"], 'x': data.coords["x"], 'y': data.coords["y"]})
                #img_array[key].append(dA)
                ###########################################################
                # RASTERIO ONLY
                data = gdal.Open(img)
                #data = rasterio.open(img)
                #print(data.ReadAsArray())
                ### ARRAY ###
                img_array[key].append(data.ReadAsArray())
                #img_array[key].append(data.read())
                ### COORDS ###
                # Order: cols, rows, proj, geom
                proj_info = data.RasterXSize, data.RasterYSize, data.GetProjection(), data.GetGeoTransform()
                #cols = data.RasterXSize
                #rows = data.RasterYSize
                #proj = data.GetProjection()
                #geom = data.GetGeoTransform()
                #print(cols, rows, proj, geom)
                ##### Is it worth making this in to a function?! and arr to img function too
                img_coords[key].append(proj_info)
                #img_coords[key].append(data.bounds)
                del data
                ############################################################
                '''
        # If all co-ordinates match in the list, select the first one.
        for key in img_coords.keys():
            if all(img_coords[key]):
                img_coords[key] = [img_coords[key][0]]

        # Stack all available arrays and calculate their cumulative sum.
        for key in img_array.keys():
            sum_arr = np.sum(np.stack(img_array[key]), axis=0)
            img_array[key] = sum_arr

            #cum_array = xr.DataArray(np.sum(np.stack(img_array[key]), axis=0), dims=["band", "x", "y"], coords={'band': img_coords[key][0]["band"], 'x': img_coords[key][0]["x"], 'y': img_coords[key][0]["y"]})

        # Export stacked tile arrays to tifs files.
        for key in img_array.keys():
            output_fp = "/".join((args.filepath_tile, "99_outputs/"))
            out_file = os.path.join(output_fp + args.time_start.replace("/", "") + "_" + args.time_end.replace("/", "") + "_" + key + ".tif")
            if not os.path.exists(output_fp):
                os.mkdir(output_fp)
            array_to_img(img_array[key], img_coords[key][0])

            sys.exit()









        for arr, coords in zip(img_array['h18v15'], img_coords['h18v15']):
            lon = [coords[0], coords[2]]
            lat = [coords[1], coords[3]]
            xtest = xr.DataArray(data=img_array['h18v15'], coords=dict(lon=(["xmin", "xmax"], lon), lat=(["ymin", "ymax"], lat)))
        print(xtest)
        sys.exit()

        ############
        # Write to netCDF4 file for each tile and extracting the appropriate coords for that image.
        ############










        sys.exit()

        ####################
        img_coords = {}
        for key, value in img_dict.items():
            if not key in img_coords:
                img_coords[key] = []
            if isinstance(value, list):
                for v in value:
                    #List order: xmin, xmax, xpixel, ymin, ymax, ypixel
                    xycoords = extract_coords(v)
                    img_coords[key].append(xycoords)

        for ik, ck in zip(img_dict.keys(), img_coords.keys()):
            if ik == ck:
                for iv, cv in zip(img_dict[ik], img_coords[ck]):
                    # function that takes img and coordinates and outputs an georeferenced array
                    import rasterio
                    data = rasterio.open(iv)
                    print(data.bounds)
                    print(cv)

                    sys.exit()
                    array = np.array([np.array(Image.open(iv))])
                    print(array)
                    ### gdal info the img before opening to an array ###
                    sys.exit()
            else:
                raise RuntimeError("Keys do not match.")

            
            
            print(d)
            print(c)
        sys.exit()



        pprint.pprint(img_dict.keys())
        pprint.pprint(img_coords.keys())
        sys.exit()
        array = np.array([np.array(Image.open(fimg)) for fimg in imagery])

        print(len(array))
        print('--------')
        print(img_coords)
        print(len(img_coords))
        sys.exit()
        cumu_array = np.sum(np.stack(array), axis=0)

        import matplotlib.pyplot as plt
        plt.imshow(cumu_array)
        plt.colorbar()
        plt.show()
        #print(cumu_array)
        sys.exit()



        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
