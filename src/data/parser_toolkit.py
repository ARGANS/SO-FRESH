#!/usr/bin/env/python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Argument parser for the Automated Polynya Identification Tool.
"""

# Package loader # 
import os, sys
from datetime import date, datetime, timedelta

# Toolkit loader #
from data.acquire import acquire_toolkit as acq_tk
from data.adjust import adjust_toolkit as adj_tk


class argument_receiver:
    def __init__(self, process, data_directory, product, aoi, start_date, end_date, resample):
        self.process = process
        self.data_directory = data_directory
        self.product = product
        self.aoi = aoi
        self.start_date = start_date
        self.end_date = end_date
        self.resample = resample

    def acquire_adjust_parser(self):
        # Format inputs.
        PF = parameter_formatting(self.aoi, self.product, self.start_date, self.end_date)
        aoi = None
        if all("MYD" in p for p in self.product):aoi = PF.modis_aoi_formater()
        elif all("SIC" in p for p in self.product):
            raise RuntimeError("The download setup for AMSR-2 sea-ice concentration is not setup, please head over to the European Space Angency's Climate Change Initiative (CCI) toolbox.")
        acquire_dates = PF.acquire_date_formater()
        adjust_dates = PF.adjust_date_formater()

        # Acquire data links & generate output directories.
        AM = acq_tk.acquire_modis(self.data_directory, acquire_dates, self.product[0], aoi)
        LPDAAC = acq_tk.lpdaac()
        product_url = AM.product_url()
        items, output_dirs = AM.url_extractor(product_url)
        # Login & setup workspace for LPDAAC data download.
        netrcdir, urs = LPDAAC.earthdata_authentication()
        # Execute download and pre-processing.
        LPDAAC.execute_download(self.data_directory, items, output_dirs, adjust=True)

        print(items)
        print(output_dirs)
        #print(self.data_directory)
        #print(acquire_dates)
        #print(self.product)
        #print(aoi)

    
    def acquire_parser(self):
        # Format inputs.
        PF = parameter_formatting(self.aoi, self.product, self.start_date, self.end_date)
        aoi = None
        if all("MYD" in p for p in self.product):aoi = PF.modis_aoi_formater()
        elif all("SIC" in p for p in self.product):
            raise RuntimeError("The download setup for AMSR-2 sea-ice concentration is not setup, please head over to the European Space Angency's Climate Change Initiative (CCI) toolbox.")
        acquire_dates = PF.acquire_date_formater()
        adjust_dates = PF.adjust_date_formater()

        # Acquire data links & generate output directories.
        AM = acq_tk.acquire_modis(self.data_directory, acquire_dates, self.product[0], aoi)
        LPDAAC = acq_tk.lpdaac()
        product_url = AM.product_url()
        items, output_dirs = AM.url_extractor(product_url)
        # Login & setup workspace for LPDAAC data download.
        netrcdir, urs = LPDAAC.earthdata_authentication()
        # Execute download.
        download = LPDAAC.execute_download(self.data_directory, items, output_dirs, adjust=False)
        '''print("inputs, acq")
        print(self.process)
        print(self.data_directory)
        print(self.product)
        print(self.aoi)
        print(self.start_date)
        print(self.end_date)
        print(self.resample)'''

    def adjust_parser(self):
        print("inputs, adj")
        print(self.process)
        print(self.data_directory)
        print(self.product)
        print(self.aoi)
        print(self.start_date)
        print(self.end_date)
        print(self.resample)

    def fusion_parser(self):
        print("inputs, fusion")
        print(self.process)
        print(self.data_directory)
        print(self.product)
        print(self.aoi)
        print(self.start_date)
        print(self.end_date)
        print(self.resample)

class parameter_formatting:
    def __init__(self, aoi, product, start_date, end_date):
        self.aoi = aoi
        self.product = product
        self.start_date = start_date
        self.end_date = end_date

    def modis_aoi_formater(self):

        """ Turn the AOI in to readable format for the products. """

        # Identify product to set aoi to correct formatting.
        if all("MYD" in p for p in self.product):
            formatting = "MODIS"
        elif any("MYD" in p for p in self.product) and any("SIC" in p for p in self.product):
            formatting = "MODIS, AMSR2"

        # Based on sepecified aoi, assign information.
        if formatting == "MODIS":
            if self.aoi == "antarctica":
                hmin, hmax, vmin, vmax = 14, 24, 15, 17
                return(hmin, hmax, vmin, vmax)
            elif self.aoi == "arctic":
                hmin, hmax, vmin, vmax = 13, 23, 0, 2
                return(hmin, hmax, vmin, vmax)

    def acquire_date_formater(self):

        """ Identify all dates between start and end entries for download. """
        
        dates=[]
        # Pull start date information.
        syear=int(self.start_date[0:4])
        smonth=int(self.start_date[5:7])
        sday=int(self.start_date[8:10])
        d0 = date(syear, smonth, sday)
        # Pull end date information.
        eyear=int(self.end_date[0:4])
        emonth=int(self.end_date[5:7])
        eday=int(self.end_date[8:10])
        d1 = date(eyear, emonth, eday)
        delta=d1-d0
        # Convert to list of datetime formats.
        for i in range(delta.days+1):
            day=d0+timedelta(days=i)
            dates.append(str(day).replace("-", ".")+"/")

        return(dates)

    def adjust_date_formater(self):

        """ Identify all dates between start and end entries for pre-processing. """

        sdate=datetime.strptime(os.path.join(self.start_date), "%Y/%m/%d").date()
        edate=datetime.strptime(os.path.join(self.end_date), "%Y/%m/%d").date()
        dates=[sdate+timedelta(days=x) for x in range((edate-sdate).days+1)]
        dates = ["/".join((str(d.year), str("%02d" %d.month), str("%02d" %d.day))) for d in dates]

        return(dates)