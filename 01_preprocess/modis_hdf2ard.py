#!/usr/bin/env python
#Author: James Hickson
#Date: 02/07/21

import argparse
import os, sys
import gdal
import glob
import shutil
#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Convert a MODIS MYD09GA HDF (AQUA) file to a image with bands 1-7 from the sensor.
**************************************************************************
##Tasks:
- Extract data from HDF files and select the desired files which include bands 1-7.
- Normalise the data - this is done by multiplying the data by 0.0001.
- From the normalised data, reconstruct the image to a multiband vrt or tif.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
def extract_hdf(hdf):
    if not os.path.exists(os.path.split(hdf)[0] + "/tmp/"):
        os.makedirs(os.path.split(hdf)[0] + "/tmp/")
    output = os.path.split(hdf)[0] + "/tmp/" + os.path.basename(os.path.splitext(hdf)[0]) + ".tif"
    cmd = "gdal_translate -of GTiff -sds %s %s"%(hdf, output)
    os.system(cmd)
    return(output)

def normalise(input, len):
    cols, rows, proj, geom = gdal.Open(input).RasterXSize, gdal.Open(input).RasterYSize, gdal.Open(input).GetProjection(), gdal.Open(input).GetGeoTransform()
    cmd = (gdal.Open(input).ReadAsArray()) * 0.0001
    if not os.path.exists(os.path.split(input)[0] + "/normalised/"):
        os.makedirs(os.path.split(input)[0] + "/normalised/")
    output = str(os.path.split(input)[0] + "/normalised/" + os.path.basename(os.path.splitext(input)[0])[:7] + "_" + os.path.basename(os.path.splitext(input)[0])[9:16] + "_band" + len + "_" + os.path.basename(os.path.splitext(input)[0])[17:23] + ".tif")
    out_img = gdal.GetDriverByName("GTiff").Create(output, cols, rows, 1, gdal.GDT_Float32)
    out_img.SetGeoTransform(geom)
    out_img.SetProjection(proj)
    outBand = out_img.GetRasterBand(1)
    outBand.SetDescription("Band " + len)
    outBand.WriteArray(cmd)
    return(output)

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-hdf", required=True, nargs="+", help="Input MODIS HDF files.")
        parser.add_argument("-o", "--output-img", required=True, nargs="+", help="Output image (Bands 1 - 7) filepath and name, ensure the order is the same as the inputs. ")
        parser.add_argument("-r", "--remove-files", required=False, action="store_true", help="Include this command to remove the 'tmp' folder created.")
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------
        for inHDF in args.input_hdf:
            if not os.path.exists(inHDF):raise RuntimeError(f"Input HDF:{args.input_hdf} not found. Check filepath and try again.")
        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        for hdf, output in zip(args.input_hdf, args.output_img):
            # Extract components of the hdf file into a temporary folder.
            # Select tif files which end with number between 12 - 18, these are bands 1 - 7.
            band_imgs = []
            for imgs in glob.glob(os.path.splitext(extract_hdf(hdf))[0] + "*.tif"):
                if (int((imgs).rsplit("_", 1)[1][:-4]) == 12) or (int((imgs).rsplit("_", 1)[1][:-4]) == 13) or (int((imgs).rsplit("_", 1)[1][:-4]) == 14) or (int((imgs).rsplit("_", 1)[1][:-4]) == 15) or (int((imgs).rsplit("_", 1)[1][:-4]) == 16) or (int((imgs).rsplit("_", 1)[1][:-4]) == 17) or (int((imgs).rsplit("_", 1)[1][:-4]) == 18):
                    band_imgs.append(imgs)
                else:
                    continue
            # Normalise the MODIS data by multiplying the data by 0.0001.
            img_list = []
            for img, i in zip(band_imgs, range(len(band_imgs))):
                img_list.append(normalise(img, str(i+1)))
            # Build the output in either vrt or tif format.
            if os.path.splitext((os.path.basename(output)))[1] == ".vrt":
                os.system("gdalbuildvrt -separate -overwrite -resolution highest -srcnodata 0 -vrtnodata 0 %s %s"%(output, " ".join(img_list)))
            elif os.path.splitext((os.path.basename(output)))[1] == ".tif":
                input = os.path.split(output)[0] + "/" + os.path.splitext((os.path.basename(output)))[0] + ".vrt"
                os.system("gdalbuildvrt -separate -overwrite -resolution highest -srcnodata 0 -vrtnodata 0 %s %s"%(input, " ".join(img_list)))
                os.system("gdal_translate -of GTiff %s %s"%(input, output))
                os.remove(input)
            else:
                raise RuntimeError("Please use '.vrt' or '.tif' as the output file extension.")

            # Remove the tmp files.
            if args.remove_files == False:
                pass
            else:
                shutil.rmtree(os.path.split(hdf)[0] + "/tmp/")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
