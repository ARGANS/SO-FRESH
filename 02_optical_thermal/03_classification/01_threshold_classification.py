#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

import argcomplete, argparse
from argcomplete.completers import ChoicesCompleter, FilesCompleter
import gdal
import sys, os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
from PIL import Image
import osr
from tqdm import tqdm
from skimage.filters import threshold_otsu, threshold_minimum
from sklearn.neighbors import KernelDensity
#--------------------------------------------------------------------------------
# Script description:
#--------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="""
# Undertake a classification of fused data. 
**************************************************************************
##Tasks:
-  Threshold classification where the upper limit value is set.
**************************************************************************""",
formatter_class=argparse.RawDescriptionHelpFormatter)

#----------------------------------------------------------------------------------------------------
# Functions:
#----------------------------------------------------------------------------------------------------
class thresholding:
    def __init__(self, data_array, default_threshold=180, max_thresh=185, min_thresh=170):
        self.data_array = data_array
        self.default_threshold = default_threshold
        self.min_thresh = min_thresh
        self.max_thresh = max_thresh
        self.x = None
        self.y = None

    def otsu_threshold(self, bins):
        # Thresholding using Otsu's method
        try:
            t = threshold_otsu(image=array, nbins=bins)
            if self.min_thresh < t < self.max_thresh:
                return t
            else:
                return self.default_threshold
        except:
            return self.default_threshold
    def min_threshold(self, bins=256):
        try:
            t = threshold_minimum(image=array, nbins=bins)
            if self.min_thresh < t < self.max_thresh:
                return t
            else:
                return self.default_thresh
        except:
            return self.default_thresh

def greyscale(array):
    # Returns greyscale of input image (all array bands are summed and divided by total number)
    return np.sum(array, axis=0)/(array.shape[0]) 

def threshold_classify(input, opt_thresh=None, opt_bands=None, therm_thresh=None, therm_bands=None, n_bands=None):
    out_img = os.path.split(input)[0]+"/03a"+os.path.basename(os.path.splitext(input)[0])[2:]+"_classification_o"+str(opt_thresh)+"_t"+str(therm_thresh)+".tif"
    tmpout = out_img[:-4] + "_tmpfile.tif"
    # Assign letters based on band numbers.
    if opt_bands or therm_bands != None:
        opt_bands_let = [chr(ord("@")+x) for x in opt_bands]
        therm_bands_let = [chr(ord("@")+x) for x in therm_bands]
    # All band numbers and corresponding letters.
    band_numbs = sorted(list(opt_bands+therm_bands))
    band_lets = sorted(list(opt_bands_let+therm_bands_let))
    # Command inputs.
    band_expression = ' '.join([f"-{l} {input} --{l}_band={n}" for n, l in zip(band_numbs, band_lets)]) 

    opt_calc = "*".join([f"({b}>1)*({b}<={opt_thresh})" for b in opt_bands_let])
    therm_calc = "*".join([f"({b}>={therm_thresh})" for b in therm_bands_let])
    # Calculation expression.
    calc_expression = opt_calc+" + "+therm_calc

    os.system("gdal_calc.py --quiet %s --outfile=%s --calc='%s' --overwrite"%(band_expression, tmpout, calc_expression))
    os.system("gdal_translate -q -b 1 %s %s"%(tmpout, out_img))
    os.remove(tmpout)

def img_viewer(input):
    img = mpimg.imread(input)
    imgplot = plt.imshow(img)
    plt.show()
#==========================================================
# main:
#----------------------------------------------------------
if __name__ == "__main__":
    try:
        #----------------------------------------------------------------------------------------------------
        # Arguments:
        #----------------------------------------------------------------------------------------------------
        parser.add_argument("-i", "--input-img", nargs="+", required=True, help="Input image - image to be classified. (Will begin with '02'.)").completer = FilesCompleter(allowednames=(".tif"))
        #parser.add_argument("-t", "--threshold", required=False, type=int, help="The upper threshold value as an interger.")
        argcomplete.autocomplete(parser)
        args = parser.parse_args()

        #----------------------------------------------------------------------------------------------------
        # Check for Errors:
        #----------------------------------------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------
        # Code:
        #----------------------------------------------------------------------------------------------------
        print(f"Undertaking threshold classification for {len(args.input_img)} images")
        for img in tqdm(args.input_img):
            products = [string for string in (os.path.basename(img).rsplit("_")) if "MYD" in string]
            if len(products) == 2:
                if gdal.Open(img).ReadAsArray().ndim == 3:
                    img_array = gdal.Open(img).ReadAsArray()
                    # Number of bands.
                    bands_n = img_array.shape[0]
                for p in enumerate(products):
                    if p[1] == "MYD09GA":
                        opt_arr = (img_array)[p[0]:int(p[0]+3)]
                        #opt_arr_slice = (int(p[0]), int(p[0]+3))
                        opt_arr_slice = [x for x in range(int(p[0]+1), int(p[0]+4))]
                        opt_greyscale = greyscale(opt_arr)
                    elif p[1] == "MYDTBGA":
                        if opt_arr.any():
                            therm_arr = (img_array)[int(p[0]+2):int(p[0]+3)]
                            #therm_arr_slice = (int(p[0]+2), int(p[0]+3))
                            therm_arr_slice = [x for x in range(int(p[0]+3), int(p[0]+4))]
                        else:
                            therm_arr = (img_array)[p[0]:int(p[0]+1)]
                            #therm_arr_slice = (int(p[0]), int(p[0]+1))
                            therm_arr_slice = [x for x in range(int(p[0]+1), int(p[0]+2))]


            init = thresholding(opt_greyscale)
            otsu_val = init.otsu_threshold(bins=128)

            threshold_classify(img, opt_thresh=otsu_val, opt_bands=opt_arr_slice, therm_thresh=265, therm_bands=therm_arr_slice, n_bands=bands_n)
        print("Process complete - threshold classifications produced")
        #----------------------------------------------------------------------------------------------------
        # Run and errors:
        #----------------------------------------------------------------------------------------------------
    except RuntimeError as msg:
        print("\033[31mERROR:\033[0m", msg)
