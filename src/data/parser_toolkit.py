#!/usr/bin/env/python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Argument parser for the Automated Polynya Identification Tool.
"""

# Package loader # 
import os, sys
import glob
from datetime import date, datetime, timedelta

# Toolkit loader #
from data.acquire import acquire_toolkit as acq_tk
from data.adjust import adjust_toolkit as adj_tk
from data.fusion import fusion_toolkit as fus_tk
from run import apit_toolkit as apit_tk


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
        # Execute download and pre-processing.#
        # MYD09GA - Reprojection & optional: resampling & mosaicing.
        # MYDTBGA - Extraction, normalising, reprojection & optional: resampling & mosaicing.
        # Check lists are the same length so that incorrect data is not selected.
        if len(output_dirs) == len(items):
            # Turn in to dictionary for efficient downloading for the mosaicing process.
            dir_dict, url_dict = adj_tk.myd_adjust.sort_inputs_by_date(adjust_dates, items, output_dirs)
        else:
            raise RuntimeError(f"The length of the download directories to download product to not match. Directories={len(output_dirs)} & products={len(items)}")

        for url, dir in zip(url_dict.keys(), dir_dict.keys()):
            if not url == dir:
                raise RuntimeError("The key order of the input dictionary is not the same.")
            else:
                download = LPDAAC.execute_download(self.data_directory, url_dict[url], dir_dict[dir], resample=self.resample, adjust=True)
                mosaic = adj_tk.myd_adjust.mosaic(self.data_directory, download, self.product, url, self.aoi) #url=date
                print("\nMosaic created for:", url, "\n")

    
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
        download = LPDAAC.execute_download(self.data_directory, items, output_dirs, resample=self.resample, adjust=False)

    def adjust_parser(self):
        # Format inputs.
        PF = parameter_formatting(self.aoi, self.product, self.start_date, self.end_date)
        aoi = None
        if all("MYD" in p for p in self.product):
            aoi = PF.modis_aoi_formater()
            adjust_dates = PF.adjust_date_formater()
            if self.product[0] == "MYD09GA":product = "MYD09GA_061"
            elif self.product[0] == "MYDTBGA":product = "MYDTBGA_006"
            tiles = adj_tk.myd_adjust.aoi_tile_identifier(aoi[0], aoi[1], aoi[2], aoi[3])
            for d in adjust_dates:
                # Pull all files from those tiles+date folders begining with "02" and ending with "tif"
                files_list = [item for sublist in [glob.glob(self.data_directory+"02_data/MODIS/"+product+"/01_tiles/"+t+"/"+d+"/02*.tif") for t in tiles] for item in sublist]
                mosaic = adj_tk.myd_adjust.mosaic(self.data_directory, files_list, self.product, d, self.aoi)
        elif all("SIC" in p for p in self.product):
            # SIC data extraction, projection, resampled (if specified) and then masked.
            adjust_dates = PF.adjust_date_formater()
            aoi = PF.amsr2_aoi_formater()
            for d in adjust_dates:
                dir = (self.data_directory+"02_data/AMSR2/01_raw/"+d[:4]+"/"+d[5:-3]+"/")
                for file in glob.glob(dir+f"ESACCI*{aoi}-{d.replace('/', '')}*.nc") + glob.glob(dir+f"ice_conc_{aoi.lower()}*{d.replace('/', '')}*.nc"):
                    print(f"Processing {d} AMSR-2 Sea ice concentration data...\n")
                    AA = adj_tk.amsr2_adjust(file, d, self.resample)
                    extract = AA.extract_from_netcdf()
                    rpjct = AA.reproject(extract)
                    if bool(self.resample) == True:
                        # Resample to a given resolution.
                        resample = adj_tk.amsr2_adjust.resample(rpjct, self.resample[0], self.resample[1])
                        adj_tk.amsr2_adjust.mask(resample)
                    else:
                        adj_tk.amsr2_adjust.mask(rpjct)


    def fusion_parser(self):

        PF = parameter_formatting(self.aoi, self.product, self.start_date, self.end_date)
        dates = PF.adjust_date_formater()
        DF = fus_tk.data_fusion(self.data_directory, dates, self.product, self.aoi)
        DF.fusion()
        
    
    def run_parser(self):

        PF = parameter_formatting(self.aoi, self.product, self.start_date, self.end_date)
        # Organise dates.
        dates = PF.adjust_date_formater()
        AP = apit_tk.apit_preparation(self.data_directory, dates, self.product, self.aoi)
        # Pull data that fits all criteria.
        data = AP.data_collector()
        AE = apit_tk.apit_execute(self.data_directory, dates, self.product, self.aoi, data)
        APIT = AE.apit()
        sys.exit()

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

    def amsr2_aoi_formater(self):

        """ Turn the AOI in to readable format for the products. """

        # Identify product to set aoi to correct formatting.
        if all("SIC" in p for p in self.product):
            if self.aoi == "antarctica":
                aoi = "SH"
                return(aoi)
            elif self.aoi == "arctic":
                aoi = "NH"
                return(aoi)

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