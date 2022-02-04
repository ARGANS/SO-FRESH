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
# Create a .txt file with all .jpg file paths from LPDAAC
**************************************************************************
##Tasks:
- User selects year (saves on memory and processing time)
- Scrape all the file paths meeting the conditions select
- Save all file paths into a .txt
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Hard commands
#----------------------------------------------------------------------------------------------------
# URL of LPDAAC website
# I have hard coded this in to reduce errors:
# https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/
# https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.006/

urlV006 = "https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.006/"
urlV061 = "https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/"

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

        parser.add_argument('-s','--startDate', required=True, help='Start YYYY-MM-DD')
        parser.add_argument('-e','--endDate', required=True, help='End YYYY-MM-DD')
        parser.add_argument('-o','--outpath', required=True, help='Out put for text file')
        parser.add_argument('-v','--version', required=True, help='Select version of modis product [006 or 061]')
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # scripting
        #----------------------------------------------------------------------------------------------------

        #Empty list for all dates wanted
        datesSelected = []

        #Turn start date into datetime format
        startYear = int(args.startDate[0:4])
        startMonth = int(args.startDate[5:7])
        startDay = int(args.startDate[8:10])
        d0 = date(startYear, startMonth, startDay)

        #Turn end date into datetime format
        endYear = int(args.endDate[0:4])
        endMonth = int(args.endDate[5:7])
        endDay = int(args.endDate[8:10])
        d1 = date(endYear, endMonth, endDay)

        #Find each date in between start and end and append to list
        delta = d1 - d0
        for i in range(delta.days + 1):
            day = d0 + timedelta(days=i)
            datesSelected.append(str(day).replace('-','.')+'/')

        # Read the URL text
        if args.version == '006':
            response = requests.get(urlV006).text
            mainUrl = urlV006
        elif args.version == '061':
            response = requests.get(urlV061).text
            mainUrl = urlV061

        # text file output name
        textFileName =os.path.join(args.outpath,'modis_%s_%s_V%s.txt'%(args.startDate,args.endDate,args.version))

        textFile = open(textFileName, 'w')
        #Appened the date to the URL to access every dateFolder
        for date in tqdm(datesSelected):
            #print(mainUrl + date)
            datePage = requests.get(mainUrl + date).text
            #print(mainUrl + date)
            soup = BeautifulSoup(datePage,"lxml")
            for link in soup.select("a[href$='.jpg']"):
                jpgPath = (mainUrl + date + link.get('href'))
                if int(jpgPath.rsplit("/", 1)[1].rsplit(".", 7)[3][4:]) >=14:
                    textFile.write("%s\n"%jpgPath)

        textFile.close()
        #----------------------------------------------------------------------------------------------------
        # Run and errors
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\nERROR: ", msg)