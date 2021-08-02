#!/usr/bin/env python
import argparse
import requests
from bs4 import BeautifulSoup
import sys
from tqdm import tqdm
import os

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
# I have hard coded this in to reduce errors but remember to change the different versions between .061 or .006:
# https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.061/
# https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.006/

mainUrl = "https://e4ftl01.cr.usgs.gov/MOLA/MYD09GA.006/"
version = mainUrl[-4:-1]
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
        parser.add_argument('year',help='Year of imagery required ex. 2002')
        parser.add_argument('txtOutPath',help='Out put for text file')
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # scripting
        #----------------------------------------------------------------------------------------------------
        # Read the URL text
        response = requests.get(mainUrl).text

        #Empty list for dates to be stored
        dateFolder = []

        #Iterate thorugh first page to get all dates
        soup = BeautifulSoup(response,"lxml")
        for link in soup.select("a[href$='/']"):
            if link.get('href')[0].isdigit() and link.get('href').startswith(args.year):
                dateFolder.append(link.get('href'))

        # text file output name
        textFileName =os.path.join(args.txtOutPath,'%s_Modis_Img_V%s.txt'%(args.year,version))

        textFile = open(textFileName, 'w')
        #Appened the date to the URL to access every dateFolder
        print('Year being processed: ',args.year)
        for date in tqdm(dateFolder):
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
