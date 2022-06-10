#!/usr/bin/env python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Toolkit containing the Automated Polynya Identification Tool data-fusion steps. 
"""

# Package loader # 
import os, sys, glob, functools
import gdal, gc


class data_fusion:
    def __init__(self, data_directory, dates, products, aoi):
        self.data_directory = data_directory
        self.dates = dates
        self.products = products
        self.aoi = aoi

    def fusion(self):

        """ Fuse together data of given products """

        for d in self.dates:
            print(f"Processing and fusing imagery from: {d}\n")
            imagery = []
            for p in self.products:
                dir = self.data_directory+"02_data/"
                if p == "MYD09GA":
                    imagery.append(glob.glob(dir+f"MODIS/{p}_061/02_mosaic/{d}/02_*{self.aoi}.tif"))
                elif p == "MYDTBGA":
                    imagery.append(glob.glob(dir+f"MODIS/{p}_006/02_mosaic/{d}/02_*{self.aoi}.tif"))
                elif p == "SIC":
                    if self.aoi=="arctic": hem = "NH"
                    elif self.aoi=="antarctica": hem = "SH"
                    imagery.append(glob.glob(dir+f"AMSR2/02_sic/{d}/02_*{hem}*.tif"))
            imgs=functools.reduce(lambda x,y:x+y,(imagery))


            ####################################################################################
            ### Check resolution of all images and make sure they are the same resolution! ###
            # function, takes all images, gdal.Open, rasterXsize, rasterYsize, make sure they're all the same?
            ####################################################################################

            if len(imgs) != len(self.products):
                raise RuntimeError(f"Based on entries, there is missing data from products for {d}. Please acquire and adjust missing data prior to re-trying fusion.")
            else:
                dir = dir+f"fusion/{'_'.join(self.products)}/{d}/"
                if not os.path.isdir(dir): os.makedirs(dir)
                outfile = dir + f"02_{'_'.join(self.products)}_{''.join(d.split('/'))}_{self.aoi}.tif"
                os.system("gdal_merge.py -o -q -of GTIFF -seperate -ot Float32 -o %s %s"%(outfile, " ".join(imgs)))
                self.set_band_description(outfile, self.description_builder(self.products))

    @staticmethod
    def set_band_description(file, band_names):

        """ Assign band names to corresponding bands """

        ds = gdal.Open(file, gdal.GA_Update)
        for b in band_names:
            flat = [item for sublist in b for item in sublist]
            band, desc = flat[0], flat[1]
            rb = ds.GetRasterBand(band)
            rb.SetDescription(desc)
        
        del ds
        gc.collect()
    
    @staticmethod
    def description_builder(products):

        """ Builds description to be passed for the band namer. """

        prods = data_fusion.band_selector(products)
        calculation=[]
        for p in products:
            if p == "MYD09GA":
                bands = ["Red", "Green", "Blue"]
                MYD09GA_numbers = [x for x in range(int(prods[p][0]+1), int(prods[p][1]+1))]
                #MYD09GA_calc = " ".join([f'{n} "Band {n}: {b}"' for n, b in zip(MYD09GA_numbers, bands)])
                MYD09GA_calc = [[[int(f"{n}")], [f"Band {n}: {b}"]] for n, b in zip(MYD09GA_numbers, bands)]
                calculation.append(MYD09GA_calc)
            elif p == "MYDTBGA":
                MYDTBGA_numbers = [x for x in range(int(prods[p][0]+1), int(prods[p][1]+1))]
                #MYDTBGA_calc = " ".join([f'{n} "Band {n}: Thermal"' for n in MYDTBGA_numbers])
                MYDTBGA_calc = [[[int(f"{n}")], [f"Band {n}: Thermal"]] for n in MYDTBGA_numbers]
                calculation.append(MYDTBGA_calc)
            elif p == "SIC":
                SIC_numbers = [x for x in range(int(prods[p][0]+1), int(prods[p][1]+1))]
                #SIC_calc = " ".join([f'{n} "Band {n}: SIC"' for n in SIC_numbers])
                SIC_calc = [[[int(f"{n}")], [f"Band {n}: SIC"]] for n in SIC_numbers]
                calculation.append(SIC_calc)
        calculation = [item for sublist in calculation for item in sublist]

        return(calculation)

    @staticmethod
    def band_selector(products):

        """ Based on order entry, find out which bands correspond to which data. """

        products_dict = dict.fromkeys(products)
        count=0
        for p in products:
            if p == "MYD09GA":
                start=count
                count=count+3
                products_dict[p]=start,count
            elif p == "MYDTBGA" or "SIC":
                start=count
                count=count+1
                products_dict[p]=(start,count)
        
        return(products_dict)
