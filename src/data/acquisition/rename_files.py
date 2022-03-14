#!/usr/bin/env python
import argparse
import os, sys

#--------------------------------------------------------------------------------
# Description of script:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Rename all files within given directory
**************************************************************************
##Tasks:
- User supllies a directory containing raw modis data downloaded from LPDAAC.
- Full filepaths or list of filepaths can be provided to the directory containing the data.
- Renames those files to start with '01_'.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Hard commands:
#----------------------------------------------------------------------------------------------------
def rename(path):
    # Take MODIS raw product (optical and thermal) and add 01 to the begining for easy file distinguishability.
    name = "01_"
    for root, dirs, files in os.walk(path):
        for i in files:
            if i.startswith("BROWSE") and i.endswith(".jpg"):
                os.rename(os.path.join(root, i), os.path.join(root,name+i))
            elif i.startswith("MYDTBGA") and i.endswith(".hdf"):
                os.rename(os.path.join(root, i), os.path.join(root,name+i))

#==========================================================
# main
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-d", "--inDir", required=True, nargs="+", help="Filepath to the directory containing images to be renamed. Can also be a filepath to multiple files (i.e. 'file/path/*/*/*/')")
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        for file in args.inDir:
            rename(file)
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\nERROR: ", msg)
