#!/usr/bin/env python
import argparse
import requests
from bs4 import BeautifulSoup
import sys
from tqdm import tqdm
import os
from datetime import date, timedelta

#--------------------------------------------------------------------------------
# Description of script
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Create a .txt file with all file paths from LP DAAC which fit within specified criteria. 
**************************************************************************
##Tasks:
- User selects year, product type and tiles.
- Scrape all the file paths meeting the conditions.
- Save all file paths into a .txt.
  NOTE: Products: optical=MYD09GA & thermal=MYDTBGA
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Hard commands
#----------------------------------------------------------------------------------------------------
# URLs of LPDAAC website:

# MYD09GA - Daily MODIS Aqua surface spectral reflectance.
## Version 6 and 6.1. (Version 6 soon to be retired (07/02/2022).
### https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.006/
### https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/

# MYDTBGA - Daily MODIS Aqua brightness temperature from thermal bands.
## Version 6.
### https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/

MYD09GA_urlV006 = "https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.006/"
MYD09GA_urlV061 = "https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/"
MYDTBGA_urlV006 = "https://e4ftl01.cr.usgs.gov/MOLA/MYDTBGA.006/"

#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-s","--time-start", required=True, help="Start of date range. Format: YYYY-MM-DD")
        parser.add_argument("-e","--time-end", required=True, help="End of date range. Format: YYYY-MM-DD")
        parser.add_argument("-o","--outpath", required=True, help="Output file for text file")
        parser.add_argument("-v","--version", required=True, help="Select version of MODIS product [006 or 061]")
        parser.add_argument("-t", "--data-type", required=True, type=str, help="Required data type (i.e. optical or thermal)")
        parser.add_argument("-hmin", "--horizontal-minimum", required=True, help="Minimum horizontal MODIS sinusoidal tile")
        parser.add_argument("-hmax", "--horizontal-maximum", required=True, help="Maximum horizontal MODIS sinusoidal tile")
        parser.add_argument("-vmin", "--vertical-minimum", required=True, help="Minimum vertical MODIS sinusoidal tile")
        parser.add_argument("-vmax", "--vertical-maximum", required=True, help="Maximum vertical MODIS sinusoidal tile")
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for errors:
        #----------------------------------------------------------------------------------------------------
        if args.data_type == "thermal" and args.version == "061":
            raise RuntimeError("MODIS product MYDTBGA does not have a V061, please use V006.")
        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        # Empty list for all dates wanted
        datesSelected = []

        # Format start date.
        startYear = int(args.time_start[0:4])
        startMonth = int(args.time_start[5:7])
        startDay = int(args.time_start[8:10])
        d0 = date(startYear, startMonth, startDay)

        # Format end date.
        endYear = int(args.time_end[0:4])
        endMonth = int(args.time_end[5:7])
        endDay = int(args.time_end[8:10])
        d1 = date(endYear, endMonth, endDay)

        # Find each date in between start and end and append to list
        delta = d1 - d0
        for i in range(delta.days + 1):
            day = d0 + timedelta(days=i)
            datesSelected.append(str(day).replace('-','.')+'/')

        if args.data_type == "optical":
            # Read the URL text
            if args.version == '006':
                response = requests.get(MYD09GA_urlV006).text
                mainUrl = MYD09GA_urlV006
            elif args.version == '061':
                response = requests.get(MYD09GA_urlV061).text
                mainUrl = MYD09GA_urlV061
        elif args.data_type == "thermal":
            response = requests.get(MYDTBGA_urlV006).text
            mainUrl = MYDTBGA_urlV006

        # Output txt name.
        if args.data_type == "optical":
            product = "MYD09GA"
        elif args.data_type == "thermal":
            product = "MYDTBGA"

        textFileName =os.path.join(args.outpath,'modis_%s_%s_%s_V%s.txt'%(product, args.time_start, args.time_end, args.version))
        if not os.path.exists(textFileName):
            os.system("touch %s"%(textFileName))
        textFile = open(textFileName, 'w')

        #Appened the date to the URL to access every dateFolder
        for date in tqdm(datesSelected):
            datePage = requests.get(mainUrl + date).text
            soup = BeautifulSoup(datePage, "lxml")
            if args.data_type == "optical":
                for link in soup.select("a[href$='.jpg']"):
                    jpgPath = (mainUrl + date + link.get('href'))
                    h = jpgPath.rsplit("/", 1)[1].rsplit(".", 7)[3][1:3]
                    v = jpgPath.rsplit("/", 1)[1].rsplit(".", 7)[3][4:]
                    if h >= args.horizontal_minimum and h <= args.horizontal_maximum:
                        if v >= args.vertical_minimum and v <= args.vertical_maximum:
                            textFile.write("%s\n"%jpgPath)
            elif args.data_type == "thermal":
                for link in soup.select("a[href$='.hdf']"):
                    hdfPath = (mainUrl + date + link.get('href'))
                    h = hdfPath.rsplit("/", 1)[1].rsplit(".", 7)[2][1:3]
                    v = hdfPath.rsplit("/", 1)[1].rsplit(".", 7)[2][4:]
                    if h >= args.horizontal_minimum and h <= args.horizontal_maximum:
                        if v >= args.vertical_minimum and v <= args.vertical_maximum:
                            textFile.write("%s\n"%hdfPath)
        textFile.close()
        #----------------------------------------------------------------------------------------------------
        # Run and errors
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\nERROR: ", msg)
