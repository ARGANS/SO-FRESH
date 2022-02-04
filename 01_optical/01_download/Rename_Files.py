#!/usr/bin/env python
import argparse
import os

#--------------------------------------------------------------------------------
# Description of script
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Rename all files within given directory
**************************************************************************
##Tasks:
User supllies a directory containing RAW modis .jpg
Script will iterate thorugh all sub directories to find files starting with 'BROWSE' endsing with '.jpg'
Renames those files to start with '01_'
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Hard commands
#----------------------------------------------------------------------------------------------------
def rename(path):
    name = '01_'
    for root, dirs, files in os.walk(path):
            for i in files:
                if i.startswith('BROWSE') and i.endswith('.jpg'):
                    os.rename(os.path.join(root, i), os.path.join(root,name+i))

#==========================================================
#main
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        print('')

        #----------------------------------------------------------------------------------------------------
        # Arguments used
        #----------------------------------------------------------------------------------------------------
        #Required arguments
        parser.add_argument('-d', '--inDir', required=True, help='Dirctory containing files that need to be renamed')
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # scripting
        #----------------------------------------------------------------------------------------------------
        path = args.inDir
        rename(path)
        
        #----------------------------------------------------------------------------------------------------
        # Run and errors
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\nERROR: ", msg)