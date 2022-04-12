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
