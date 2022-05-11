#!/usr/bin/env python3
# @James Hickson | Argans UK | jhickson@argans.co.uk
"""
Toolkit for data fusion application for the Automated Polynya Identification Tool
"""
# packages
import os, sys
import gdal, glob, itertools, osr, functools
import numpy as np
import pandas as pd
import scipy
import skimage
from skimage import morphology
from datetime import datetime, timedelta

def fusion(data_folder, sdate, edate, products):

    """ Fuse together data of given products """
    
    sdate=datetime.strptime(os.path.join(sdate), "%Y/%m/%d").date()
    edate=datetime.strptime(os.path.join(edate), "%Y/%m/%d").date()
    date_dt=[sdate+timedelta(days=x) for x in range((edate-sdate).days+1)]
    dates = ["/".join((str(d.year), str("%02d" %d.month), str("%02d" %d.day))) for d in date_dt]
    if any("MYD09GA" in p for p in products):
        version = input("What MYD09GA version would you like, 006 or 061?:\n")
    for d in dates:
        print(f"Processing: {d}...")
        imgs_4_fusion=[]
        for p in products:
            if p == "MYD09GA": 
                fp=glob.glob((data_folder+"MODIS/"+p+"_"+version+"/02_mosaic/"+d+"/02_*.tif"))
            elif p == "MYDTBGA":fp=glob.glob((data_folder+"MODIS/"+p+"_006"+"/02_mosaic/"+d+"/02_*.tif"))
            elif p == "SIC": fp =glob.glob((data_folder+"AMSR2/sic_extracted/"+d+"/02_*tif"))
            imgs_4_fusion.append(fp)
        imgs=functools.reduce(lambda x,y:x+y,(imgs_4_fusion))
        if len(imgs) == len(products):
            if not os.path.isdir(data_folder+"fusion/"+"_".join(products)+"/"+d): os.makedirs(data_folder+"fusion/"+"_".join(products)+"/"+d)
            outfile=data_folder+"fusion/"+"_".join(products)+"/"+d+"/02_"+"_".join(products)+"_"+"".join(d.split("/"))+"_ANTARCTICA.tif"        
            os.system("gdal_merge.py -o -q -of GTIFF -seperate -ot Float32 -o %s %s"%(outfile, " ".join(imgs)))
            #os.system("gdal_set_band_desc.py %s %s"%(outfile, description_builder(products)))
            set_band_descriptions(outfile, description_builder(products))

def description_builder(products):

    """ Builds description to be passed for the band namer. """

    prods = band_selector(products)
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

def set_band_descriptions(filepath, bands):

    """ Assigns band names to corresponding bands """

    ds = gdal.Open(filepath, gdal.GA_Update)
    for b in bands:
        flat = [item for sublist in b for item in sublist]
        band, desc = flat[0], flat[1]
        rb = ds.GetRasterBand(band)
        rb.SetDescription(desc)
    del ds

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